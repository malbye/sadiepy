"""Wave damping across a vegetation patch (theory.md §3).

:func:`propagate` integrates the energy balance (theory.md Eq. 3.1) along a patch of
layered vegetation. With constant drag coefficients it uses the analytic hyperbolic
decay (Eq. 3.7); otherwise it integrates Eq. (3.9) with ``scipy.integrate.solve_ivp``.
"""

import math
import warnings
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Literal

import numpy as np
import numpy.typing as npt
from scipy.integrate import quad, solve_ivp
from scipy.optimize import least_squares

from sadiepy.core.drag import DragClosure, constant
from sadiepy.core.vegetation import Layer, Vegetation
from sadiepy.core.waves import RegularWave, keulegan_carpenter, reynolds, ursell
from sadiepy.units import GRAVITY, Water
from sadiepy.warnings import ExtrapolationWarning

__all__ = ["PropagationResult", "propagate", "s_function"]

FloatArray = npt.NDArray[np.float64]
Method = Literal["auto", "analytic", "ode"]

_N_FIT = 201  # points along the patch for the effective-beta fit
_ODE_RTOL = 1e-11
_QUAD_EPSREL = 1e-12


def s_function(x: float) -> float:
    """S(x) = sinh³x + 3 sinh x, with S'(x) = 3 cosh³x (theory.md Eq. 3.5)."""
    sh = math.sinh(x)
    return sh**3 + 3.0 * sh


@dataclass(frozen=True)
class PropagationResult:
    """Result of :func:`propagate`.

    Attributes
    ----------
    x : numpy.ndarray
        Requested positions along the patch, in m from the leading edge.
    H : numpy.ndarray
        Wave height at ``x``, in m.
    beta : float
        Damping factor, in 1/m: exact (theory.md Eq. 3.7) for the analytic method,
        otherwise a least-squares fit of Eq. (3.7) to H(x).
    beta_fit_rms : float
        RMS of (H_fit - H) / H0 over the patch; 0 for the analytic method.
    layer_fractions : dict of str to float
        Share of the total energy dissipated over the patch in each layer.
    cumulative_dissipation : numpy.ndarray
        Energy dissipated per unit crest width between 0 and ``x``, in W/m.
    energy_flux : numpy.ndarray
        Energy flux E c_g at ``x``, in W/m (theory.md Eq. 2.12).
    method : {"analytic", "ode"}
        Solution method used.
    drag_mode : {"bulk", "per-layer"}
        Whether one closure was applied to every layer or each layer had its own.
    drag_used : dict of str to str
        Closure name used in each layer.
    extrapolated : bool
        True if any closure was evaluated outside its valid range.
    dimensionless : dict
        Range (min, max) of Ur over the patch, and of Re and KC for each closure that
        uses them, with the velocity and length scales used.
    warnings : tuple of str
        Validity and extrapolation messages (also emitted as warnings).
    """

    x: FloatArray
    H: FloatArray
    beta: float
    beta_fit_rms: float
    layer_fractions: dict[str, float]
    cumulative_dissipation: FloatArray
    energy_flux: FloatArray
    method: Literal["analytic", "ode"]
    drag_mode: Literal["bulk", "per-layer"]
    drag_used: dict[str, str]
    extrapolated: bool
    dimensionless: dict[str, object]
    warnings: tuple[str, ...]


def _layer_integral(layer: Layer, k: float, depth: float) -> float:
    """∫ b(s) cosh³(ks) ds over the submerged part of a layer, in m² (theory.md Eq. 3.4).

    Analytic (Eq. 3.5) for constant width, ``scipy.integrate.quad`` otherwise.
    """
    span = layer.submerged_range(depth)
    if span is None:
        return 0.0
    s_b, s_t = span
    if layer.has_constant_width:
        width = layer.width_at(s_b)
        return width * (s_function(k * s_t) - s_function(k * s_b)) / (3.0 * k)
    value, _ = quad(
        lambda s: layer.width_at(s) * math.cosh(k * s) ** 3,
        s_b,
        s_t,
        epsabs=0.0,
        epsrel=_QUAD_EPSREL,
        limit=200,
    )
    return float(value)


def _closure_variable(closure: DragClosure, wave: RegularWave, height: float, nu: float) -> float:
    """Value of the closure's variable for the local wave height (theory.md Eq. 2.13)."""
    if closure.variable == "constant":
        return math.nan
    if closure.variable == "Ur":
        return ursell(height, wave.L, wave.depth)
    u_bar = wave.omega * 0.5 * height / wave.kh  # theory.md Eq. (2.9)
    length = closure.length_scale
    assert length is not None  # guaranteed by DragClosure.__post_init__
    if closure.variable == "Re":
        return reynolds(u_bar, length, nu)
    return keulegan_carpenter(u_bar, wave.T, length)


def _resolve_closures(
    vegetation: Vegetation, drag: DragClosure | float | None
) -> tuple[list[DragClosure], Literal["bulk", "per-layer"]]:
    """One closure per layer, and whether they come from a bulk closure."""
    if drag is not None:
        if any(layer.drag is not None for layer in vegetation.layers):
            raise ValueError(
                "pass either a bulk `drag` to propagate() or per-layer `Layer.drag`, not both"
            )
        closure = drag if isinstance(drag, DragClosure) else constant(float(drag))
        return [closure] * len(vegetation.layers), "bulk"
    closures = []
    for layer in vegetation.layers:
        if layer.drag is None:
            raise ValueError(
                f"layer {layer.name!r} has no drag closure; set Layer.drag or pass `drag`"
            )
        closures.append(layer.drag)
    return closures, "per-layer"


def propagate(
    wave: RegularWave,
    vegetation: Vegetation,
    patch_length: float,
    water: Water | None = None,
    drag: DragClosure | float | None = None,
    x_out: Sequence[float] | npt.ArrayLike | None = None,
    method: Method = "auto",
) -> PropagationResult:
    """Propagate a regular wave across a patch of layered vegetation.

    Implements theory.md §3: dissipation per layer (Eq. 3.4), hyperbolic decay for
    constant C_D (Eq. 3.7), and the ODE (Eq. 3.9) when C_D depends on the wave.

    Parameters
    ----------
    wave : RegularWave
        Incident wave at the leading edge (x = 0); its height is H0, in m.
    vegetation : Vegetation
        Layered vegetation; heights in m above the bed.
    patch_length : float
        Length of the patch in the direction of propagation, in m.
    water : Water, optional
        Water properties; default ``Water()`` (fresh water at 20 °C).
    drag : DragClosure or float, optional
        Bulk drag closure applied to every layer (a float means a constant C_D).
        If omitted, every layer must have its own ``Layer.drag``.
    x_out : array_like, optional
        Positions in m, within ``[0, patch_length]``, at which to report results.
        Default ``[0, patch_length]``.
    method : {"auto", "analytic", "ode"}, optional
        ``"auto"`` uses the analytic solution when every closure is constant and the
        ODE otherwise. ``"ode"`` forces numerical integration (for testing).

    Returns
    -------
    PropagationResult
        Heights, damping factor, layer fractions and diagnostics.

    Raises
    ------
    ValueError
        If inputs are inconsistent, or ``method="analytic"`` with a wave-dependent
        closure.
    """
    water = water or Water()
    if not (math.isfinite(patch_length) and patch_length > 0.0):
        raise ValueError(f"patch_length must be positive, got {patch_length}")
    x = np.atleast_1d(np.asarray([0.0, patch_length] if x_out is None else x_out, dtype=np.float64))
    if np.any(x < 0.0) or np.any(x > patch_length):
        raise ValueError("x_out must lie within [0, patch_length]")

    closures, drag_mode = _resolve_closures(vegetation, drag)
    all_constant = all(c.is_constant for c in closures)
    if method == "auto":
        method = "analytic" if all_constant else "ode"
    if method == "analytic" and not all_constant:
        raise ValueError("method='analytic' requires constant drag coefficients")

    rho, nu, k, h0 = water.rho, water.nu, wave.k, wave.H
    flux_coeff = rho * GRAVITY * wave.cg / 8.0  # E c_g = flux_coeff * H² (Eqs. 2.11-2.12)
    # epsilon_i = C_D,i * g_i * H³ (theory.md Eq. 3.4 with u_m from Eq. 2.8)
    velocity_factor = (wave.omega / (2.0 * math.sinh(wave.kh))) ** 3
    g_layer = np.array(
        [
            (2.0 / (3.0 * math.pi))
            * rho
            * layer.n_per_m2
            * velocity_factor
            * _layer_integral(layer, k, wave.depth)
            for layer in vegetation.layers
        ]
    )

    messages = list(wave.validity_issues())
    x_fit = np.linspace(0.0, patch_length, _N_FIT)

    if method == "analytic":
        cd = np.array([c() for c in closures])
        a_layers = cd * g_layer  # A_i in theory.md Eq. (3.6)
        a_total = float(a_layers.sum())
        beta = 4.0 * a_total * h0 / (rho * GRAVITY * wave.cg)  # theory.md Eq. (3.7)
        heights = h0 / (1.0 + beta * x)
        h_fit_grid = h0 / (1.0 + beta * x_fit)
        if beta > 0.0:
            # ∫_0^x A H0³ / (1 + beta x')³ dx' = A H0³ [1 - (1 + beta x)^-2] / (2 beta)
            cumulative = a_total * h0**3 * (1.0 - (1.0 + beta * x) ** -2) / (2.0 * beta)
        else:
            cumulative = np.zeros_like(x)
        fractions_raw = a_layers
        beta_rms = 0.0
        variable_values: list[list[float]] = [[] for _ in closures]
    else:
        n_layers = len(closures)

        def layer_dissipation(height: float) -> FloatArray:
            cds = np.array([c(_closure_variable(c, wave, height, nu)) for c in closures])
            result: FloatArray = cds * g_layer * height**3
            return result

        def rhs(_x: float, y: FloatArray) -> FloatArray:
            height = float(y[0])
            eps = layer_dissipation(height)
            dh_dx = -4.0 * float(eps.sum()) / (rho * GRAVITY * wave.cg * height)  # Eq. (3.9)
            return np.concatenate(([dh_dx], eps))

        sol = solve_ivp(
            rhs,
            (0.0, patch_length),
            np.concatenate(([h0], np.zeros(n_layers))),
            method="DOP853",
            rtol=_ODE_RTOL,
            atol=1e-14 * h0,
            dense_output=True,
        )
        if not sol.success:
            raise RuntimeError(f"ODE integration failed: {sol.message}")
        states = sol.sol(x)
        heights = states[0]
        cumulative = states[1:].sum(axis=0)
        grid_states = sol.sol(x_fit)
        h_fit_grid = grid_states[0]
        fractions_raw = grid_states[1:, -1]
        beta = _fit_beta(x_fit, h_fit_grid, h0)
        beta_rms = float(np.sqrt(np.mean((h0 / (1.0 + beta * x_fit) - h_fit_grid) ** 2)) / h0)
        variable_values = [
            [_closure_variable(c, wave, float(hh), nu) for hh in h_fit_grid] for c in closures
        ]

    total = float(np.sum(fractions_raw))
    layer_fractions = {
        layer.name: (float(f) / total if total > 0.0 else 0.0)
        for layer, f in zip(vegetation.layers, fractions_raw, strict=True)
    }

    # Ranges of the dimensionless groups over the patch (theory.md §2.9).
    dimensionless: dict[str, object] = {
        "Ur": (
            ursell(float(np.min(h_fit_grid)), wave.L, wave.depth),
            ursell(float(np.max(h_fit_grid)), wave.L, wave.depth),
        ),
        "velocity_scale": "u_bar",
    }
    extrapolated = False
    reported: set[str] = set()
    for closure, values in zip(closures, variable_values, strict=True):
        if closure.is_constant or closure.name in reported:
            continue
        reported.add(closure.name)
        lo, hi = min(values), max(values)
        dimensionless[f"{closure.variable}[{closure.name}]"] = {
            "range": (lo, hi),
            "length_scale": closure.length_scale,
        }
        if not (closure.in_range(lo) and closure.in_range(hi)):
            extrapolated = True
            msg = (
                f"closure {closure.name!r}: {closure.variable} spans [{lo:.4g}, {hi:.4g}], "
                f"outside its valid range {closure.valid_range}"
            )
            messages.append(msg)
            warnings.warn(msg, ExtrapolationWarning, stacklevel=2)

    return PropagationResult(
        x=x,
        H=np.asarray(heights, dtype=np.float64),
        beta=float(beta),
        beta_fit_rms=beta_rms,
        layer_fractions=layer_fractions,
        cumulative_dissipation=np.asarray(cumulative, dtype=np.float64),
        energy_flux=flux_coeff * np.asarray(heights, dtype=np.float64) ** 2,
        method=method,
        drag_mode=drag_mode,
        drag_used={
            layer.name: c.name for layer, c in zip(vegetation.layers, closures, strict=True)
        },
        extrapolated=extrapolated,
        dimensionless=dimensionless,
        warnings=tuple(messages),
    )


def _fit_beta(x: FloatArray, heights: FloatArray, h0: float) -> float:
    """Least-squares fit of H = H0 / (1 + beta x) (theory.md Eq. 3.7) for beta, in 1/m."""
    # Initial guess from the end point: 1/H - 1/H0 = beta x / H0.
    guess = max((h0 / heights[-1] - 1.0) / x[-1], 0.0)
    fit = least_squares(
        lambda b: h0 / (1.0 + b[0] * x) - heights,
        x0=[guess],
        bounds=([0.0], [np.inf]),
        xtol=1e-15,
        ftol=1e-15,
        gtol=1e-15,
    )
    return float(fit.x[0])

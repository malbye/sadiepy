"""Tests for sadiepy.core.solver (theory.md §3): spec M1 tests 2-7."""

import math
import warnings
from collections.abc import Callable

import numpy as np
import pytest
from scipy.integrate import quad, trapezoid

from sadiepy import ExtrapolationWarning, Layer, RegularWave, Vegetation, Water, propagate
from sadiepy.core import drag
from sadiepy.units import GRAVITY


def stems(top: float = 0.4, d: float = 0.01, n: float = 36.0, **kw: object) -> Layer:
    return Layer(name="stem", s_bottom=0.0, s_top=top, n_per_m2=n, frontal_width=d, **kw)  # type: ignore[arg-type]


def zhang_eq4(cd: float, d: float, n: float, wave: RegularWave, h_vs: float) -> float:
    """Zhang et al. (2023) Eq. (4), p. 3 = theory.md Eq. (3.8)."""
    k, h = wave.k, wave.depth
    return (
        4.0
        / (9.0 * math.pi)
        * cd
        * d
        * n
        * wave.H
        * k
        * (math.sinh(k * h_vs) ** 3 + 3 * math.sinh(k * h_vs))
        / ((math.sinh(2 * k * h) + 2 * k * h) * math.sinh(k * h))
    )


# --- Test 2: no vegetation ------------------------------------------------------


def test_no_layers_gives_constant_height() -> None:
    wave = RegularWave(H=0.08, T=1.4, depth=0.4)
    res = propagate(wave, Vegetation(layers=()), 2.0, drag=1.0, x_out=[0.0, 1.0, 2.0])
    assert res.beta == 0.0
    np.testing.assert_array_equal(res.H, [0.08, 0.08, 0.08])


def test_layer_above_water_dissipates_nothing() -> None:
    wave = RegularWave(H=0.08, T=1.4, depth=0.3)
    dry = Layer("canopy", s_bottom=0.35, s_top=0.5, n_per_m2=36.0, frontal_width=0.1)
    res = propagate(wave, Vegetation(layers=(dry,)), 2.0, drag=1.0)
    assert res.beta == 0.0


# --- Test 3: Zhang Eq. (4) ------------------------------------------------------


@pytest.mark.parametrize(
    ("depth", "top", "h_vs"),
    [(0.3, 0.4, 0.3), (0.4, 0.4, 0.4), (0.5, 0.4, 0.4), (0.5, 0.25, 0.25)],
)
def test_single_uniform_layer_matches_zhang_eq4(depth: float, top: float, h_vs: float) -> None:
    wave = RegularWave(H=0.08, T=1.4, depth=depth)
    res = propagate(wave, Vegetation(layers=(stems(top=top),)), 2.0, drag=1.3)
    assert res.method == "analytic"
    assert res.beta == pytest.approx(zhang_eq4(1.3, 0.01, 36.0, wave, h_vs), rel=1e-10)


# --- Test 4: splitting a layer -------------------------------------------------


def test_splitting_layer_leaves_beta_unchanged() -> None:
    wave = RegularWave(H=0.08, T=1.4, depth=0.5)
    whole = Vegetation(layers=(stems(top=0.4),))
    split = Vegetation(
        layers=(
            Layer("lower", 0.0, 0.13, 36.0, 0.01),
            Layer("upper", 0.13, 0.4, 36.0, 0.01),
        )
    )
    b_whole = propagate(wave, whole, 2.0, drag=1.1).beta
    b_split = propagate(wave, split, 2.0, drag=1.1).beta
    assert b_split == pytest.approx(b_whole, rel=1e-12)


def test_beta_is_sum_of_layer_betas() -> None:
    wave = RegularWave(H=0.08, T=1.4, depth=0.5)
    roots = Layer("root", 0.0, 0.1, 1008.0, 0.005, drag=drag.constant(2.0))
    stem = Layer("stem", 0.0, 0.4, 36.0, 0.01, drag=drag.constant(1.0))
    both = propagate(wave, Vegetation(layers=(roots, stem)), 2.0)
    each = [propagate(wave, Vegetation(layers=(layer,)), 2.0).beta for layer in (roots, stem)]
    assert both.beta == pytest.approx(sum(each), rel=1e-12)
    assert both.drag_mode == "per-layer"
    assert both.layer_fractions["root"] == pytest.approx(each[0] / sum(each), rel=1e-12)


# --- Test 5: numeric (callable) vs analytic width -------------------------------


def test_callable_width_matches_constant_width() -> None:
    wave = RegularWave(H=0.08, T=1.4, depth=0.5)
    analytic = Layer("stem", 0.0, 0.4, 36.0, 0.01)
    numeric = Layer("stem", 0.0, 0.4, 36.0, lambda s: 0.01)
    b_a = propagate(wave, Vegetation(layers=(analytic,)), 2.0, drag=1.0).beta
    b_n = propagate(wave, Vegetation(layers=(numeric,)), 2.0, drag=1.0).beta
    assert b_n == pytest.approx(b_a, rel=1e-8)


def test_triangular_canopy_against_hand_integral() -> None:
    # Inverted triangle, width 0 at s=0.23 rising to W=0.145 at s=0.40 (fully submerged).
    wave = RegularWave(H=0.08, T=1.4, depth=0.5)
    s0, s1, w = 0.23, 0.40, 0.145
    canopy = Layer("canopy", s0, s1, 36.0, lambda s: w * (s - s0) / (s1 - s0))
    res = propagate(wave, Vegetation(layers=(canopy,)), 2.0, drag=1.0)
    # Same integral by a fine trapezoidal rule, then theory.md Eqs. (3.4), (3.7).
    s = np.linspace(s0, s1, 200_001)
    integral = float(trapezoid(w * (s - s0) / (s1 - s0) * np.cosh(wave.k * s) ** 3, s))
    rho = Water().rho
    a = (2 / (3 * math.pi)) * rho * 36.0 * (wave.omega / (2 * math.sinh(wave.kh))) ** 3 * integral
    beta = 4 * a * wave.H / (rho * GRAVITY * wave.cg)
    assert res.beta == pytest.approx(beta, rel=1e-8)


# --- Test 6: ODE vs analytic for constant C_D -----------------------------------


def test_ode_matches_hyperbolic_decay() -> None:
    wave = RegularWave(H=0.08, T=1.4, depth=0.4)
    veg = Vegetation(layers=(Layer("root", 0.0, 0.1, 1008.0, 0.005), stems()))
    x = np.linspace(0.0, 2.0, 9)
    exact = propagate(wave, veg, 2.0, drag=1.5, x_out=x, method="analytic")
    ode = propagate(wave, veg, 2.0, drag=1.5, x_out=x, method="ode")
    assert ode.method == "ode"
    np.testing.assert_allclose(ode.H, wave.H / (1.0 + exact.beta * x), rtol=1e-6)
    assert ode.beta == pytest.approx(exact.beta, rel=1e-6)
    assert ode.beta_fit_rms < 1e-8


# --- Test 7: energy conservation ------------------------------------------------


def kc_closure() -> drag.DragClosure:
    return drag.custom("test_kc", lambda kc: 0.6 + (100.0 / kc) ** 1.2, "KC", length_scale=0.01)


def test_energy_flux_loss_equals_dissipation_constant_cd_ode() -> None:
    # Only the ODE branch is a real check: in the analytic branch the identity holds
    # by construction once beta = 4 A H0 / (rho g c_g).
    wave = RegularWave(H=0.08, T=1.4, depth=0.4)
    res = propagate(
        wave, Vegetation(layers=(stems(),)), 2.0, drag=2.0, x_out=[0.0, 2.0], method="ode"
    )
    loss = res.energy_flux[0] - res.energy_flux[-1]
    assert res.cumulative_dissipation[-1] == pytest.approx(loss, rel=1e-6)


def test_energy_flux_loss_equals_dissipation_variable_cd() -> None:
    wave = RegularWave(H=0.12, T=1.4, depth=0.4)
    res = propagate(wave, Vegetation(layers=(stems(),)), 2.0, drag=kc_closure(), x_out=[0.0, 2.0])
    assert res.method == "ode"
    loss = res.energy_flux[0] - res.energy_flux[-1]
    assert res.cumulative_dissipation[-1] == pytest.approx(loss, rel=1e-6)


# --- Variable C_D and diagnostics ---------------------------------------------


def test_variable_cd_reports_ranges_and_fit() -> None:
    wave = RegularWave(H=0.12, T=1.4, depth=0.4)
    res = propagate(wave, Vegetation(layers=(stems(),)), 2.0, drag=kc_closure())
    kc_info = res.dimensionless["KC[test_kc]"]
    assert isinstance(kc_info, dict)
    lo, hi = kc_info["range"]
    # KC at the leading edge, from theory.md Eqs. (2.9), (2.13).
    assert hi == pytest.approx(wave.u_bar * wave.T / 0.01, rel=1e-12)
    assert lo < hi
    assert res.H[-1] < res.H[0]
    assert res.drag_used == {"stem": "test_kc"}
    assert res.drag_mode == "bulk"
    assert not res.extrapolated


def test_extrapolation_warns() -> None:
    wave = RegularWave(H=0.12, T=1.4, depth=0.4)
    narrow = drag.custom("narrow", lambda kc: 1.0 + 0 * kc, "KC", 0.01, valid_range=(1.0, 2.0))
    with pytest.warns(ExtrapolationWarning, match="narrow"):
        res = propagate(wave, Vegetation(layers=(stems(),)), 2.0, drag=narrow)
    assert res.extrapolated
    assert any("narrow" in m for m in res.warnings)


def test_ode_layer_fractions_match_analytic_beta_shares() -> None:
    # The ODE integrates each layer's dissipation separately; with constant C_D its
    # shares must equal the ratio of single-layer betas (theory.md §3.7).
    wave = RegularWave(H=0.08, T=1.4, depth=0.5)
    root = Layer("root", 0.0, 0.1, 1008.0, 0.005)
    stem = stems()
    ode = propagate(wave, Vegetation(layers=(root, stem)), 2.0, drag=1.0, method="ode")
    b_root = propagate(wave, Vegetation(layers=(root,)), 2.0, drag=1.0).beta
    b_stem = propagate(wave, Vegetation(layers=(stem,)), 2.0, drag=1.0).beta
    assert ode.layer_fractions["root"] == pytest.approx(b_root / (b_root + b_stem), rel=1e-8)


def test_more_vegetation_more_damping() -> None:
    wave = RegularWave(H=0.08, T=1.4, depth=0.4)
    betas = [
        propagate(wave, Vegetation(layers=(stems(n=n),)), 2.0, drag=1.0).beta
        for n in (10.0, 36.0, 100.0)
    ]
    assert betas == sorted(betas)


# --- Input checking ---------------------------------------------------------------


def test_bulk_and_per_layer_drag_together_raise() -> None:
    wave = RegularWave(H=0.08, T=1.4, depth=0.4)
    veg = Vegetation(layers=(stems(drag=drag.constant(1.0)),))
    with pytest.raises(ValueError, match="not both"):
        propagate(wave, veg, 2.0, drag=1.0)


def test_missing_drag_raises() -> None:
    wave = RegularWave(H=0.08, T=1.4, depth=0.4)
    with pytest.raises(ValueError, match="no drag closure"):
        propagate(wave, Vegetation(layers=(stems(),)), 2.0)


def test_analytic_with_variable_cd_raises() -> None:
    wave = RegularWave(H=0.08, T=1.4, depth=0.4)
    with pytest.raises(ValueError, match="analytic"):
        propagate(wave, Vegetation(layers=(stems(),)), 2.0, drag=kc_closure(), method="analytic")


@pytest.mark.parametrize("x_out", [[-0.1], [2.1]])
def test_x_out_outside_patch_raises(x_out: list[float]) -> None:
    wave = RegularWave(H=0.08, T=1.4, depth=0.4)
    with pytest.raises(ValueError, match="x_out"):
        propagate(wave, Vegetation(layers=(stems(),)), 2.0, drag=1.0, x_out=x_out)


def test_no_warnings_for_normal_case() -> None:
    wave = RegularWave(H=0.08, T=1.4, depth=0.4)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        res = propagate(wave, Vegetation(layers=(stems(),)), 2.0, drag=1.0)
    assert res.warnings == ()


# --- First-principles dissipation (independent of theory.md Eqs. 3.3-3.6) ----------


def beta_from_first_principles(
    wave: RegularWave,
    cd: float,
    n: float,
    width: Callable[[float], float],
    s_b: float,
    s_t: float,
) -> float:
    """beta = 4 eps(H0) / (rho g c_g H0²), with eps = <∫ ½ rho C_D N b |u| u u ds>.

    u is taken from phi_x of the potential (theory.md Eq. 2.1), the period average and
    the depth integral are done by quadrature, and c_g by a numerical d(omega)/dk.
    """
    rho = Water().rho
    k, h, a, om = wave.k, wave.depth, wave.a, wave.omega
    top = min(s_t, h)

    def u(s: float, theta: float) -> float:
        return G_ * a * k / om * math.cosh(k * s) / math.cosh(k * h) * math.cos(theta)

    def per_height(s: float) -> float:
        avg = quad(lambda th: abs(u(s, th)) ** 3, 0.0, 2 * math.pi, epsrel=1e-11)[0]
        return float(0.5 * rho * cd * n * width(s) * avg / (2 * math.pi))

    eps = quad(per_height, s_b, top, epsrel=1e-11)[0]

    def omega_of_k(kk: float) -> float:
        return math.sqrt(G_ * kk * math.tanh(kk * h))

    dk = 1e-6 * k
    cg = (omega_of_k(k + dk) - omega_of_k(k - dk)) / (2 * dk)
    return float(4 * eps / (rho * G_ * cg * wave.H**2))


G_ = GRAVITY


@pytest.mark.parametrize(
    ("depth", "s_b", "s_t", "width"),
    [
        (0.4, 0.0, 0.4, lambda s: 0.01),  # emergent stem
        (0.5, 0.23, 0.40, lambda s: 0.145 * (s - 0.23) / 0.17),  # submerged canopy
        (0.3, 0.23, 0.40, lambda s: 0.145 * (s - 0.23) / 0.17),  # canopy cut by surface
        (0.5, 0.1, 0.3, lambda s: 0.02),  # layer with s_b > 0
    ],
)
def test_beta_matches_first_principles(
    depth: float, s_b: float, s_t: float, width: Callable[[float], float]
) -> None:
    wave = RegularWave(H=0.08, T=1.4, depth=depth)
    layer = Layer("layer", s_b, s_t, 36.0, width)
    res = propagate(wave, Vegetation(layers=(layer,)), 2.0, drag=1.3)
    expected = beta_from_first_principles(wave, 1.3, 36.0, width, s_b, s_t)
    assert res.beta == pytest.approx(expected, rel=1e-7)


def test_constant_layer_above_surface_matches_clipped_callable() -> None:
    wave = RegularWave(H=0.08, T=1.4, depth=0.4)
    const = Layer("c", 0.1, 0.5, 36.0, 0.02)
    func = Layer("c", 0.1, 0.5, 36.0, lambda s: 0.02)
    b_c = propagate(wave, Vegetation(layers=(const,)), 2.0, drag=1.0).beta
    b_f = propagate(wave, Vegetation(layers=(func,)), 2.0, drag=1.0).beta
    assert b_f == pytest.approx(b_c, rel=1e-8)


def test_numpy_scalar_drag_is_constant() -> None:
    wave = RegularWave(H=0.08, T=1.4, depth=0.4)
    res = propagate(wave, Vegetation(layers=(stems(),)), 2.0, drag=np.int64(1))
    assert res.beta == pytest.approx(
        propagate(wave, Vegetation(layers=(stems(),)), 2.0, drag=1.0).beta, rel=1e-15
    )

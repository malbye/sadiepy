"""Linear (Airy) wave theory for regular waves on a flat bed.

Derivations, notation and validity checks are in ``docs/theory.md`` §2. Heights are
measured upward from the bed, ``s = z + h`` (theory.md §2.1).
"""

import math
import warnings
from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt

from sadiepy.units import GRAVITY, Water
from sadiepy.warnings import ValidityWarning

__all__ = [
    "BREAKING_DEPTH_RATIO",
    "BREAKING_STEEPNESS_COEFF",
    "RegularWave",
    "keulegan_carpenter",
    "reynolds",
    "ursell",
    "wavenumber",
]

FloatArray = npt.NDArray[np.float64]

BREAKING_DEPTH_RATIO = 0.78
"""Depth-limited breaking ratio H/h (theory.md §2.10).

Dean & Dalrymple (1991), Eq. (4.122a), p. 113 (McCowan, 1894);
``refs/Dean_Dalrymple_1991.pdf``.
"""

BREAKING_STEEPNESS_COEFF = 0.142
"""Coefficient of the Miche steepness limit H/L = 0.142 tanh(kh) (theory.md §2.10).

TODO(verify): citable source for the tanh (Miche, 1944) form. Dean & Dalrymple (1991),
Eq. (12.12), p. 335 (``refs/Dean_Dalrymple_1991.pdf``) gives only the deep-water limit
with a nonlinear wavelength; see theory.md §2.10 and DECISIONS D011.
"""

_NEWTON_RTOL = 1e-14  # relative residual |f|/y at which Newton stops (theory.md §2.5)
_NEWTON_MAX_ITER = 50


def wavenumber(period: npt.ArrayLike, depth: npt.ArrayLike, g: float = GRAVITY) -> FloatArray:
    """Solve the linear dispersion relation for the wavenumber.

    Solves ``omega**2 = g k tanh(k h)`` (theory.md Eq. 2.2) by Newton iteration on
    ``x tanh x = y`` with ``x = k h``, ``y = omega**2 h / g``, starting from
    ``x0 = y / sqrt(tanh y)`` (theory.md Eqs. 2.3-2.5). Vectorised over broadcastable
    ``period`` and ``depth``.

    Parameters
    ----------
    period : array_like
        Wave period T, in s. Must be positive.
    depth : array_like
        Still-water depth h, in m. Must be positive.
    g : float, optional
        Acceleration of gravity, in m/s². Default :data:`sadiepy.units.GRAVITY`.

    Returns
    -------
    numpy.ndarray
        Wavenumber k, in rad/m, with the broadcast shape of the inputs (0-d for
        scalar inputs).

    Raises
    ------
    ValueError
        If any period or depth is not positive and finite.
    RuntimeError
        If Newton iteration does not converge (never observed for
        1e-6 <= y <= 1e4; theory.md §2.5).

    Examples
    --------
    >>> k = wavenumber(1.4, 0.4)
    >>> round(float(k), 4)
    2.6264
    """
    t_arr, h_arr = np.broadcast_arrays(
        np.asarray(period, dtype=np.float64), np.asarray(depth, dtype=np.float64)
    )
    if not (np.all(np.isfinite(t_arr)) and np.all(t_arr > 0.0)):
        raise ValueError("period must be positive and finite")
    if not (np.all(np.isfinite(h_arr)) and np.all(h_arr > 0.0)):
        raise ValueError("depth must be positive and finite")

    omega = 2.0 * np.pi / t_arr
    y = omega**2 * h_arr / g
    x = y / np.sqrt(np.tanh(y))  # theory.md Eq. (2.4)
    for _ in range(_NEWTON_MAX_ITER):
        tanh_x = np.tanh(x)
        f = x * tanh_x - y  # theory.md Eq. (2.3)
        if np.all(np.abs(f) <= _NEWTON_RTOL * y):
            break
        # sech²x = 1 - tanh²x avoids overflow in cosh for very deep entries.
        f_prime = tanh_x + x * (1.0 - tanh_x**2)
        x = x - f / f_prime  # theory.md Eq. (2.5)
    else:
        raise RuntimeError("dispersion relation did not converge")
    result: FloatArray = x / h_arr
    return result


def reynolds(velocity: float, length: float, nu: float) -> float:
    """Reynolds number Re = u l / nu (theory.md Eq. 2.13).

    Parameters
    ----------
    velocity : float
        Chosen velocity scale u, in m/s (e.g. :attr:`RegularWave.u_bar`).
    length : float
        Chosen length scale l, in m (e.g. a stem diameter).
    nu : float
        Kinematic viscosity, in m²/s (e.g. :attr:`sadiepy.Water.nu`).

    Returns
    -------
    float
        Reynolds number (dimensionless).

    References
    ----------
    Zhang et al. (2023), Eq. (5), p. 3. ``refs/Zhang_et_al_2023.pdf``.
    """
    return velocity * length / nu


def keulegan_carpenter(velocity: float, period: float, length: float) -> float:
    """Keulegan-Carpenter number KC = u T / l (theory.md Eq. 2.13).

    Parameters
    ----------
    velocity : float
        Chosen velocity scale u, in m/s.
    period : float
        Wave period T, in s.
    length : float
        Chosen length scale l, in m.

    Returns
    -------
    float
        Keulegan-Carpenter number (dimensionless).

    References
    ----------
    Zhang et al. (2023), Eq. (6), p. 3. ``refs/Zhang_et_al_2023.pdf``.
    """
    return velocity * period / length


def ursell(height: float, wavelength: float, depth: float) -> float:
    """Ursell number Ur = H L² / h³ (theory.md Eq. 2.13).

    Equal to ``8 pi**3 (H/L) / (k h)**3``, the form in Zhang et al. (2023).

    Parameters
    ----------
    height : float
        Wave height H, in m.
    wavelength : float
        Linear wavelength L, in m.
    depth : float
        Still-water depth h, in m.

    Returns
    -------
    float
        Ursell number (dimensionless).

    References
    ----------
    Zhang et al. (2023), Eq. (7), p. 3. ``refs/Zhang_et_al_2023.pdf``.
    """
    return height * wavelength**2 / depth**3


@dataclass(frozen=True)
class RegularWave:
    """A regular wave on a flat bed, described by linear wave theory.

    On construction, the wave is checked against the validity limits of theory.md
    §2.10; each violated limit emits :class:`~sadiepy.ValidityWarning`.

    Parameters
    ----------
    H : float
        Wave height, in m. Must be positive.
    T : float
        Wave period, in s. Must be positive.
    depth : float
        Still-water depth h, in m. Must be positive.
    ursell_max : float or None, optional
        Optional extra threshold: warn if the Ursell number exceeds it. ``None``
        (default) disables this check. Nonlinearity is always checked with
        theory.md Eq. (2.15) regardless.

    Examples
    --------
    >>> wave = RegularWave(H=0.08, T=1.4, depth=0.4)
    >>> round(wave.L, 3)
    2.392
    """

    H: float
    T: float
    depth: float
    ursell_max: float | None = None
    k: float = field(init=False)
    """Wavenumber, in rad/m (theory.md Eq. 2.2)."""

    def __post_init__(self) -> None:
        """Solve for k and run the validity checks of theory.md §2.10."""
        if not (math.isfinite(self.H) and self.H > 0.0):
            raise ValueError(f"H must be positive and finite, got {self.H}")
        object.__setattr__(self, "k", float(wavenumber(self.T, self.depth)))
        for message in self.validity_issues():
            warnings.warn(message, ValidityWarning, stacklevel=3)

    # --- kinematics (theory.md §2.4-2.6) -------------------------------------
    @property
    def a(self) -> float:
        """Amplitude H/2, in m."""
        return 0.5 * self.H

    @property
    def omega(self) -> float:
        """Angular frequency 2 pi / T, in rad/s."""
        return 2.0 * math.pi / self.T

    @property
    def kh(self) -> float:
        """Relative depth k h (dimensionless)."""
        return self.k * self.depth

    @property
    def L(self) -> float:
        """Wavelength 2 pi / k, in m (theory.md Eq. 2.6)."""
        return 2.0 * math.pi / self.k

    @property
    def c(self) -> float:
        """Phase speed omega / k, in m/s (theory.md Eq. 2.6)."""
        return self.omega / self.k

    @property
    def n(self) -> float:
        """Ratio of group to phase speed, (1 + 2kh / sinh 2kh) / 2 (theory.md Eq. 2.7)."""
        two_kh = 2.0 * self.kh
        return 0.5 * (1.0 + two_kh / math.sinh(two_kh))

    @property
    def cg(self) -> float:
        """Group speed n c, in m/s (theory.md Eq. 2.7)."""
        return self.n * self.c

    # --- velocities (theory.md §2.7) -----------------------------------------
    def u_m(self, s: npt.ArrayLike) -> FloatArray:
        """Horizontal orbital velocity amplitude at height s above the bed.

        Implements theory.md Eq. (2.8): ``u_m(s) = a omega cosh(k s) / sinh(k h)``.

        Parameters
        ----------
        s : array_like
            Height above the bed, in m, with ``0 <= s <= depth``.

        Returns
        -------
        numpy.ndarray
            Velocity amplitude, in m/s, with the shape of ``s``.

        Raises
        ------
        ValueError
            If any ``s`` lies outside ``[0, depth]``.
        """
        s_arr = self._check_heights(s)
        result: FloatArray = self.a * self.omega * np.cosh(self.k * s_arr) / math.sinh(self.kh)
        return result

    def orbital_excursion(self, s: npt.ArrayLike) -> FloatArray:
        """Horizontal orbital excursion amplitude A_w = u_m / omega, in m.

        Implements theory.md Eq. (2.10).

        Parameters
        ----------
        s : array_like
            Height above the bed, in m, with ``0 <= s <= depth``.

        Returns
        -------
        numpy.ndarray
            Excursion amplitude, in m, with the shape of ``s``.
        """
        result: FloatArray = self.u_m(s) / self.omega
        return result

    @property
    def u_bar(self) -> float:
        """Depth-averaged velocity amplitude omega a / (k h), in m/s.

        Implements theory.md Eq. (2.9); Zhang et al. (2023), Eq. (8), p. 3
        (``refs/Zhang_et_al_2023.pdf``).
        """
        return self.omega * self.a / self.kh

    # --- energy (theory.md §2.8) ---------------------------------------------
    def energy(self, water: Water | None = None) -> float:
        """Wave energy per unit surface area E = rho g H² / 8, in J/m².

        Implements theory.md Eq. (2.11).

        Parameters
        ----------
        water : Water, optional
            Water properties; default ``Water()`` (fresh water at 20 °C).

        Returns
        -------
        float
            Energy density, in J/m².
        """
        rho = (water or Water()).rho
        return rho * GRAVITY * self.H**2 / 8.0

    def energy_flux(self, water: Water | None = None) -> float:
        """Energy flux per unit crest width F = E c_g, in W/m.

        Implements theory.md Eq. (2.12).

        Parameters
        ----------
        water : Water, optional
            Water properties; default ``Water()`` (fresh water at 20 °C).

        Returns
        -------
        float
            Energy flux, in W/m.
        """
        return self.energy(water) * self.cg

    # --- dimensionless groups and validity (theory.md §2.9-2.10) ------------
    @property
    def steepness(self) -> float:
        """Wave steepness H/L (dimensionless)."""
        return self.H / self.L

    @property
    def relative_height(self) -> float:
        """Relative wave height H/h (dimensionless)."""
        return self.H / self.depth

    @property
    def ursell(self) -> float:
        """Ursell number H L² / h³ (theory.md Eq. 2.13)."""
        return ursell(self.H, self.L, self.depth)

    @property
    def second_order_ratio(self) -> float:
        """Ratio R of second- to first-order velocity potential (theory.md Eq. 2.14).

        Dean & Dalrymple (1991), Eq. (11.34), p. 303 (``refs/Dean_Dalrymple_1991.pdf``).
        Linear theory requires R << 1.
        """
        kh = self.kh
        return 0.375 * self.k * self.a * math.cosh(2.0 * kh) / (math.cosh(kh) * math.sinh(kh) ** 3)

    @property
    def trough_bump_limit(self) -> float:
        """Largest ka with no secondary trough crest at second order (theory.md Eq. 2.15).

        ``sinh³(kh) / (cosh(kh) (2 + cosh(2kh)))``; Dean & Dalrymple (1991),
        Eq. (11.41), p. 304 (``refs/Dean_Dalrymple_1991.pdf``).
        """
        kh = self.kh
        return math.sinh(kh) ** 3 / (math.cosh(kh) * (2.0 + math.cosh(2.0 * kh)))

    def validity_issues(self) -> list[str]:
        """List the validity limits of theory.md §2.10 that this wave violates.

        Returns
        -------
        list of str
            One message per violated limit; empty if the wave is within all limits.
        """
        issues = []
        if self.relative_height > BREAKING_DEPTH_RATIO:
            issues.append(
                f"H/h = {self.relative_height:.3g} exceeds the depth-limited breaking "
                f"ratio {BREAKING_DEPTH_RATIO} (Dean & Dalrymple 1991, Eq. 4.122a)"
            )
        steepness_limit = BREAKING_STEEPNESS_COEFF * math.tanh(self.kh)
        if self.steepness > steepness_limit:
            issues.append(
                f"H/L = {self.steepness:.3g} exceeds the steepness-limited breaking "
                f"limit 0.142 tanh(kh) = {steepness_limit:.3g} (Miche form)"
            )
        ka = self.k * self.a
        if ka > self.trough_bump_limit:
            issues.append(
                f"ka = {ka:.3g} exceeds {self.trough_bump_limit:.3g}, where second-order "
                "Stokes theory predicts a secondary crest in the trough (Dean & Dalrymple "
                f"1991, Eq. 11.41; Ur = {self.ursell:.3g}, R = {self.second_order_ratio:.3g}): "
                "nonlinearity is strong and linear wave theory is unreliable"
            )
        if self.ursell_max is not None and self.ursell > self.ursell_max:
            issues.append(
                f"Ursell number {self.ursell:.3g} exceeds ursell_max = {self.ursell_max}; "
                "linear wave theory is questionable"
            )
        return issues

    def _check_heights(self, s: npt.ArrayLike) -> FloatArray:
        """Convert heights to an array and check ``0 <= s <= depth``."""
        s_arr = np.asarray(s, dtype=np.float64)
        if np.any(s_arr < 0.0) or np.any(s_arr > self.depth):
            raise ValueError(f"heights above the bed must lie in [0, depth = {self.depth}] m")
        return s_arr

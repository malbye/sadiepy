"""Physical constants and water properties (density, viscosity).

All correlations are from Sharqawy, Lienhard & Zubair (2010),
``refs/Sharqawy_et_al_2010.pdf``; derivations and ranges in ``docs/theory.md`` §1.
"""

import warnings
from dataclasses import dataclass

from sadiepy.warnings import ValidityWarning

__all__ = [
    "GRAVITY",
    "Water",
    "pure_water_dynamic_viscosity",
    "seawater_density",
    "seawater_dynamic_viscosity",
]

GRAVITY = 9.80665
"""Standard acceleration of gravity, in m/s² (3rd CGPM, 1901; theory.md §2)."""

# Validity ranges (theory.md §1.2-1.3). Bounds are treated as inclusive (theory.md §1.6).
_T_MIN_C = 0.0  # °C, Sharqawy et al. (2010) Table 2 Eq. (8) p. 357; Table 5 Eq. (22) p. 364
_T_MAX_C = 180.0  # °C, same sources
_S_MAX_DENSITY = 0.16  # kg/kg, Sharqawy et al. (2010) Table 2, Eq. (8), p. 357
_S_MAX_VISCOSITY = 0.15  # kg/kg, Sharqawy et al. (2010) Table 5, Eq. (22), p. 364

# Density coefficients: Sharqawy et al. (2010), Table 2, Eq. (8), p. 357
# (refs/Sharqawy_et_al_2010.pdf); theory.md Eq. (1.1).
_A1, _A2, _A3, _A4, _A5 = 9.999e2, 2.034e-2, -6.162e-3, 2.261e-5, -4.657e-8
_B1, _B2, _B3, _B4, _B5 = 8.020e2, -2.001, 1.677e-2, -3.060e-5, -1.613e-5


def seawater_density(t: float, s: float) -> float:
    """Density of seawater at atmospheric pressure.

    Implements theory.md Eq. (1.1). No range checking; use :class:`Water` for that.

    Parameters
    ----------
    t : float
        Temperature on the ITS-90 scale, in °C.
    s : float
        Absolute (reference-composition) salinity, in **kg/kg** (35 g/kg = 0.035).

    Returns
    -------
    float
        Density, in kg/m³.

    References
    ----------
    Sharqawy, M.H., Lienhard V, J.H., Zubair, S.M. (2010). Desalination and Water
    Treatment 16, 354-380. Table 2, Eq. (8), p. 357. ``refs/Sharqawy_et_al_2010.pdf``.
    """
    pure = _A1 + _A2 * t + _A3 * t**2 + _A4 * t**3 + _A5 * t**4
    saline = _B1 * s + _B2 * s * t + _B3 * s * t**2 + _B4 * s * t**3 + _B5 * s**2 * t**2
    return pure + saline


def pure_water_dynamic_viscosity(t: float) -> float:
    """Dynamic viscosity of pure water at atmospheric pressure.

    Implements theory.md Eq. (1.4). No range checking.

    Parameters
    ----------
    t : float
        Temperature on the ITS-90 scale, in °C.

    Returns
    -------
    float
        Dynamic viscosity, in Pa·s (kg/(m·s)).

    References
    ----------
    Sharqawy et al. (2010), Table 5, Eq. (23), p. 364. ``refs/Sharqawy_et_al_2010.pdf``.
    """
    return 4.2844e-5 + 1.0 / (0.157 * (t + 64.993) ** 2 - 91.296)


def seawater_dynamic_viscosity(t: float, s: float) -> float:
    """Dynamic viscosity of seawater at atmospheric pressure.

    Implements theory.md Eqs. (1.2)-(1.4). No range checking; use :class:`Water`.

    Parameters
    ----------
    t : float
        Temperature on the ITS-90 scale, in °C.
    s : float
        Absolute (reference-composition) salinity, in **kg/kg** (35 g/kg = 0.035).

    Returns
    -------
    float
        Dynamic viscosity, in Pa·s (kg/(m·s)).

    References
    ----------
    Sharqawy et al. (2010), Table 5, Eqs. (22)-(23), p. 364.
    ``refs/Sharqawy_et_al_2010.pdf``.
    """
    a = 1.541 + 1.998e-2 * t - 9.52e-5 * t**2
    b = 7.974 - 7.561e-2 * t + 4.724e-4 * t**2
    return pure_water_dynamic_viscosity(t) * (1.0 + a * s + b * s**2)


@dataclass(frozen=True)
class Water:
    """Water at a given temperature and salinity, at atmospheric pressure.

    Density and viscosity follow Sharqawy et al. (2010); see theory.md §1. Inputs
    outside the correlations' validity ranges (0-180 °C; salinity up to 160 g/kg for
    density and 150 g/kg for viscosity) emit :class:`~sadiepy.ValidityWarning` and
    still return a value.

    Parameters
    ----------
    temperature : float
        Temperature on the ITS-90 scale, in °C. Default 20.
    salinity : float
        Absolute salinity, in **g/kg** (fresh water 0, typical seawater 35).
        Default 0. Converted to kg/kg internally.

    Raises
    ------
    ValueError
        If ``salinity`` is negative.

    References
    ----------
    Sharqawy et al. (2010), ``refs/Sharqawy_et_al_2010.pdf``: Table 2, Eq. (8),
    p. 357 (density); Table 5, Eqs. (22)-(23), p. 364 (viscosity).

    Examples
    --------
    >>> w = Water(temperature=20.0, salinity=0.0)
    >>> round(w.rho, 2)
    998.02
    """

    temperature: float = 20.0
    salinity: float = 0.0

    def __post_init__(self) -> None:
        """Check inputs against the correlations' validity ranges."""
        if self.salinity < 0.0:
            raise ValueError(f"salinity must be >= 0 g/kg, got {self.salinity}")
        t = self.temperature
        if not _T_MIN_C <= t <= _T_MAX_C:
            warnings.warn(
                f"temperature {t} °C is outside the 0-180 °C validity range of the "
                "Sharqawy et al. (2010) density and viscosity correlations",
                ValidityWarning,
                stacklevel=3,
            )
        s = self._s_kg_per_kg
        if s > _S_MAX_DENSITY:
            warnings.warn(
                f"salinity {self.salinity} g/kg exceeds the 160 g/kg density and "
                "150 g/kg viscosity validity limits of Sharqawy et al. (2010)",
                ValidityWarning,
                stacklevel=3,
            )
        elif s > _S_MAX_VISCOSITY:
            warnings.warn(
                f"salinity {self.salinity} g/kg exceeds the 150 g/kg validity limit of "
                "the Sharqawy et al. (2010) viscosity correlation",
                ValidityWarning,
                stacklevel=3,
            )

    @property
    def _s_kg_per_kg(self) -> float:
        """Salinity in kg/kg, the unit used by the correlations."""
        return self.salinity / 1000.0

    @property
    def rho(self) -> float:
        """Density, in kg/m³ (theory.md Eq. 1.1)."""
        return seawater_density(self.temperature, self._s_kg_per_kg)

    @property
    def mu(self) -> float:
        """Dynamic viscosity, in Pa·s (theory.md Eqs. 1.2-1.4)."""
        return seawater_dynamic_viscosity(self.temperature, self._s_kg_per_kg)

    @property
    def nu(self) -> float:
        """Kinematic viscosity, in m²/s (theory.md Eq. 1.5)."""
        return self.mu / self.rho

"""Drag-coefficient closures for vegetation (theory.md §3.8-3.9 and §5.4-5.5).

A closure gives C_D either as a constant or as a function of one dimensionless
variable (Re, KC or Ur). Besides :func:`constant` and the user-defined hook
:func:`custom`, the module holds a registry of published closures, transcribed
from Zhang et al. (2023), Table 3 (refs/Zhang_et_al_2023.pdf, p. 12); see
:data:`REGISTRY` and :func:`get`.

All published closures here are **bulk** coefficients for the whole plant
(theory.md §3.9). Their C_D values were back-calculated with Zhang et al.'s
Eq. (4) and a depth-averaged frontal area (theory.md §5.3), and their Re and KC use
the depth-averaged velocity of the incident wave (theory.md §5.4).
"""

import math
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import Literal

__all__ = [
    "REGISTRY",
    "DragClosure",
    "LengthScale",
    "Variable",
    "constant",
    "custom",
    "get",
    "names",
]

Variable = Literal["constant", "Re", "KC", "Ur"]

LengthScale = Literal["d", "l_e", "l_ew", "r_v"]
"""Named characteristic lengths, resolved from the vegetation at the wave's depth:
reference diameter, effective vegetation length, effective water length, and
vegetation hydraulic radius (theory.md Eq. 5.3)."""

_NAMED_LENGTHS = ("d", "l_e", "l_ew", "r_v")


@dataclass(frozen=True)
class DragClosure:
    """A drag-coefficient closure C_D = formula(variable).

    Parameters
    ----------
    name : str
        Identifier recorded in results.
    formula : callable
        Maps the variable's value (dimensionless) to C_D (dimensionless). For
        ``variable="constant"`` it is called with ``nan`` and must ignore it.
    variable : {"constant", "Re", "KC", "Ur"}
        Dimensionless group C_D depends on (theory.md Eq. 2.13).
    length_scale : float or {"d", "l_e", "l_ew", "r_v"} or None
        Length scale for Re and KC: a value in m, or the name of a characteristic
        length of the vegetation (theory.md Eq. 5.3). Required for Re and KC.
    velocity_scale : {"u_bar"}
        Velocity scale for Re and KC; the depth-averaged amplitude (theory.md Eq. 2.9).
    valid_range : tuple of float
        (min, max) of the variable over which the closure was fitted, at the
        closure's source scale. Evaluating outside it emits
        :class:`~sadiepy.ExtrapolationWarning` in the solver.
    scale : float or None
        Geometric scale lambda of the source experiment (1:lambda), used to convert
        Re between scales (theory.md Eq. 5.6). None if unknown or not applicable.
    vegetation_type : str
        What the closure was fitted to.
    components : str
        Plant components present in the source experiment, e.g. ``"stem+root"``.
    citation : str
        Source of the closure.
    notes : str
        Free-text notes.
    """

    name: str
    formula: Callable[[float], float]
    variable: Variable = "constant"
    length_scale: float | LengthScale | None = None
    velocity_scale: Literal["u_bar"] = "u_bar"
    valid_range: tuple[float, float] = (-math.inf, math.inf)
    scale: float | None = None
    vegetation_type: str = ""
    components: str = ""
    citation: str = ""
    notes: str = ""

    def __post_init__(self) -> None:
        """Check the length scale for Re/KC closures and the source scale."""
        if self.variable in ("Re", "KC"):
            ls = self.length_scale
            if isinstance(ls, str):
                if ls not in _NAMED_LENGTHS:
                    raise ValueError(
                        f"closure {self.name!r}: unknown length scale {ls!r}; "
                        f"use a value in m or one of {_NAMED_LENGTHS}"
                    )
            elif not (ls is not None and ls > 0.0):
                raise ValueError(
                    f"closure {self.name!r} on {self.variable} needs a length_scale > 0"
                )
        if self.scale is not None and not self.scale > 0.0:
            raise ValueError(f"closure {self.name!r}: scale must be positive")

    @property
    def is_constant(self) -> bool:
        """True if C_D does not depend on the wave."""
        return self.variable == "constant"

    def __call__(self, value: float = math.nan) -> float:
        """Evaluate C_D at the variable's value (dimensionless)."""
        return self.formula(value)

    def in_range(self, value: float) -> bool:
        """Return True if ``value`` lies inside :attr:`valid_range` (inclusive)."""
        return self.valid_range[0] <= value <= self.valid_range[1]


def constant(cd: float) -> DragClosure:
    """Create a constant drag coefficient closure.

    Parameters
    ----------
    cd : float
        Drag coefficient C_D (dimensionless). Must be positive.

    Returns
    -------
    DragClosure
        Closure named ``"constant(<cd>)"``.

    Examples
    --------
    >>> constant(1.2)()
    1.2
    """
    if not (math.isfinite(cd) and cd > 0.0):
        raise ValueError(f"cd must be positive and finite, got {cd}")
    return DragClosure(name=f"constant({cd:g})", formula=lambda _value: cd)


def custom(
    name: str,
    formula: Callable[[float], float],
    variable: Variable,
    length_scale: float | LengthScale | None = None,
    valid_range: tuple[float, float] = (-math.inf, math.inf),
    citation: str = "user-defined",
    scale: float | None = None,
) -> DragClosure:
    """Create a user-defined closure C_D = formula(variable).

    Parameters
    ----------
    name : str
        Identifier recorded in results.
    formula : callable
        Maps the variable's value to C_D (dimensionless).
    variable : {"constant", "Re", "KC", "Ur"}
        Dimensionless group C_D depends on.
    length_scale : float or str, optional
        Length scale for Re and KC: a value in m or a named characteristic length.
    valid_range : tuple of float, optional
        Fitted range of the variable; default unbounded.
    citation : str, optional
        Source of the closure.
    scale : float, optional
        Geometric scale lambda of the experiment the closure was fitted to.

    Returns
    -------
    DragClosure
    """
    return DragClosure(
        name=name,
        formula=formula,
        variable=variable,
        length_scale=length_scale,
        valid_range=valid_range,
        citation=citation,
        scale=scale,
    )


# --- Published closures: Zhang et al. (2023), Table 3 ---------------------------------
# Transcribed from refs/Zhang_et_al_2023.pdf, Table 3, p. 12 (formula, scope, scale,
# vegetation). Subscripts: Re_d and KC_d use the stem diameter d (DBH); Re_rv and
# KC_rv the vegetation hydraulic radius r_v (Zhang Eq. 11). Original sources are in
# the reference list of Zhang et al. (2023).

_TABLE3 = "Zhang et al. (2023), Table 3 (refs/Zhang_et_al_2023.pdf, p. 12)"
_STEM_ROOT = "artificial mangrove with stem, root"
_STEM_ROOT_CANOPY = "artificial mangrove with stem, root, canopy"


def _power_law(a1: float, a2: float, a3: float) -> Callable[[float], float]:
    """C_D = a1 + (a2 / phi)**a3 (Zhang et al. 2023, p. 12, form of Table 3 fits)."""

    def formula(phi: float) -> float:
        return float(a1 + (a2 / phi) ** a3)

    return formula


def _maza2017(re: float) -> float:
    return float(9.49 * re**-0.31)


def _he2019(kc: float) -> float:
    return 18.025 * math.exp(-0.043 * kc)


_CLOSURES = (
    DragClosure(
        "maza2017_re",
        _maza2017,
        "Re",
        "d",
        valid_range=(500.0, 1700.0),
        scale=12.0,
        vegetation_type=_STEM_ROOT,
        components="stem+root",
        citation=f"Maza et al. (2017), via {_TABLE3}",
    ),
    DragClosure(
        "he2019_kc",
        _he2019,
        "KC",
        "d",
        valid_range=(10.0, 37.0),
        scale=10.0,
        vegetation_type=_STEM_ROOT_CANOPY,
        components="stem+root+canopy",
        citation=(
            f"He, Chen & Jiang (2019), Coastal Eng. 152, 103509, Eq. (15), p. 6 "
            f"(refs/He_et_al_2019.pdf); via {_TABLE3}"
        ),
        notes=(
            "He et al. define KC = u_max T / b_v (Eqs. 6-7, p. 2), with u_max the orbital "
            "velocity at the still-water surface and b_v the plant area per unit height; "
            "their roots are an 8 cm polyethylene dome, not rods. Zhang et al. Table 3 "
            "lists the fit as KC_d (u_bar and stem diameter), which SADIE follows; with "
            "He's own definitions KC would differ."
        ),
    ),
    DragClosure(
        "wang2022b_kc",
        _power_law(0.42, 0.77, 0.41),
        "KC",
        "r_v",
        valid_range=(0.01, 0.28),
        scale=10.0,
        vegetation_type=_STEM_ROOT,
        components="stem+root",
        citation=(
            "Wang, Yin & Liu (2022b), Appl. Ocean Res. 118, 102974, Eq. (7), p. 11 "
            f"(refs/Wang_et_al_2022b.pdf); via {_TABLE3}"
        ),
        notes=(
            "Wang et al. define KC_rv = u_p T / r_v (p. 11) with the pore velocity "
            "u_p = u_m / (1 - phi), u_m the maximum orbital velocity at the vegetation top "
            "(p. 2), r_v = (pi/4)(1 - phi)/phi b_v and phi = pi N b_v^2 / 4, where b_v = d* = "
            "A_f / (N h) for mangrove models (p. 6). SADIE follows Zhang et al. (u_bar and "
            "r_v from DBH and SVF, Eq. 11), so KC_rv differs from Wang's. Wang states the "
            "range 0.02 <= KC_rv <= 0.28 (p. 12); Table 3 gives 0.01-0.28."
        ),
    ),
    DragClosure(
        "chang2022_kc",
        _power_law(1.15, 14.70, 1.065),
        "KC",
        "d",
        valid_range=(2.5, 20.0),
        scale=7.0,
        vegetation_type=_STEM_ROOT,
        components="stem+root",
        citation=f"Chang et al. (2022), via {_TABLE3}",
    ),
    DragClosure(
        "chang2022_re",
        _power_law(1.15, 1015.21, 1.044),
        "Re",
        "d",
        valid_range=(180.0, 1300.0),
        scale=7.0,
        vegetation_type=_STEM_ROOT,
        components="stem+root",
        citation=f"Chang et al. (2022), via {_TABLE3}",
    ),
    DragClosure(
        "bryant2022_kc",
        _power_law(0.86, 1.24, 0.78),
        "KC",
        "d",
        valid_range=(0.0, 12.0),
        scale=2.1,
        vegetation_type=_STEM_ROOT,
        components="stem+root",
        citation=f"Bryant et al. (2022), via {_TABLE3}",
    ),
    DragClosure(
        "kelty2022_re",
        _power_law(0.60, 300000.0, 1.0),
        "Re",
        "d",
        valid_range=(4900.0, 119000.0),
        scale=1.0,
        vegetation_type=_STEM_ROOT,
        components="stem+root",
        citation=f"Kelty et al. (2022), via {_TABLE3}",
    ),
    DragClosure(
        "zhang2023_kc",
        _power_law(0.60, 101.99, 1.22),
        "KC",
        "d",
        valid_range=(6.05, 64.90),
        scale=10.0,
        vegetation_type=_STEM_ROOT_CANOPY,
        components="stem+root+canopy",
        citation=_TABLE3,
        notes="Fitted to M1-M3 together. The Conclusions give KC_d 12.97-77.39 (erratum E5).",
    ),
    DragClosure(
        "zhang2023_rerv",
        _power_law(0.60, 420000.0, 1.41),
        "Re",
        "r_v",
        valid_range=(34900.0, 407000.0),
        scale=10.0,
        vegetation_type=_STEM_ROOT_CANOPY,
        components="stem+root+canopy",
        citation=_TABLE3,
        notes="Fitted to M1-M3 together.",
    ),
    DragClosure(
        "zhang2023_ur",
        _power_law(0.60, 80.00, 0.75),
        "Ur",
        valid_range=(0.73, 21.33),
        scale=10.0,
        vegetation_type=_STEM_ROOT_CANOPY,
        components="stem+root+canopy",
        citation=_TABLE3,
        notes="Fitted to M1-M3 together.",
    ),
)

REGISTRY: Mapping[str, DragClosure] = MappingProxyType({c.name: c for c in _CLOSURES})
"""Published closures by name (read-only)."""


def names() -> list[str]:
    """Names of the registered closures, in the order of Zhang et al. (2023), Table 3."""
    return list(REGISTRY)


def get(name: str) -> DragClosure:
    """Look up a registered closure by name.

    Parameters
    ----------
    name : str
        A key of :data:`REGISTRY`, e.g. ``"zhang2023_kc"``.

    Returns
    -------
    DragClosure

    Examples
    --------
    >>> round(get("zhang2023_kc")(20.0), 3)
    7.898
    """
    try:
        return REGISTRY[name]
    except KeyError:
        raise KeyError(f"unknown drag closure {name!r}; available: {names()}") from None

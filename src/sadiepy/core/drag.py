"""Drag-coefficient closures for vegetation (theory.md §3.8-3.9).

A closure gives C_D either as a constant or as a function of one dimensionless
variable (Re, KC or Ur). M1 provides :func:`constant` and a user-defined hook
(:func:`custom`); the published closures of Zhang et al. (2023), Table 3, are added
in M2.
"""

import math
from collections.abc import Callable
from dataclasses import dataclass
from typing import Literal

__all__ = ["DragClosure", "Variable", "constant", "custom"]

Variable = Literal["constant", "Re", "KC", "Ur"]


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
    length_scale : float or None
        Length scale l for Re and KC, in m. Required for those variables.
    velocity_scale : {"u_bar"}
        Velocity scale for Re and KC; the depth-averaged amplitude (theory.md Eq. 2.9).
    valid_range : tuple of float
        (min, max) of the variable over which the closure was fitted. Evaluating
        outside it emits :class:`~sadiepy.ExtrapolationWarning` in the solver.
    citation : str
        Source of the closure.
    notes : str
        Free-text notes.
    """

    name: str
    formula: Callable[[float], float]
    variable: Variable = "constant"
    length_scale: float | None = None
    velocity_scale: Literal["u_bar"] = "u_bar"
    valid_range: tuple[float, float] = (-math.inf, math.inf)
    citation: str = ""
    notes: str = ""

    def __post_init__(self) -> None:
        """Check that Re/KC closures have a positive length scale."""
        if self.variable in ("Re", "KC") and not (
            self.length_scale is not None and self.length_scale > 0.0
        ):
            raise ValueError(f"closure {self.name!r} on {self.variable} needs a length_scale > 0")

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
    length_scale: float | None = None,
    valid_range: tuple[float, float] = (-math.inf, math.inf),
    citation: str = "user-defined",
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
    length_scale : float, optional
        Length scale for Re and KC, in m.
    valid_range : tuple of float, optional
        Fitted range of the variable; default unbounded.
    citation : str, optional
        Source of the closure.

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
    )

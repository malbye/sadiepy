"""Vegetation described as a vertical stack of layers (theory.md §3.2).

Heights are measured upward from the bed (theory.md §2.1). M1 covers the geometry
needed for drag dissipation (frontal width, density, vertical extent); volumes and
characteristic lengths (l_e, l_ew, r_v) are added in M2.
"""

import math
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from sadiepy.core.drag import DragClosure

__all__ = ["Layer", "Vegetation", "Width"]

Width = float | Callable[[float], float]
"""Frontal width of one element: a constant in m, or a function b(s) of height s (m)."""


@dataclass(frozen=True)
class Layer:
    """One layer of vegetation elements, e.g. roots, stems or canopy.

    Parameters
    ----------
    name : str
        Layer name (unique within a :class:`Vegetation`).
    s_bottom, s_top : float
        Bottom and top of the layer, in m above the bed; ``0 <= s_bottom < s_top``.
        Only the part below the still-water level dissipates energy.
    n_per_m2 : float
        Number of elements per m² of bed. Must be positive.
    frontal_width : float or callable
        Frontal (projected) width of one element, in m: a constant, or a function
        b(s) of the height above the bed, in m, defined on ``[s_bottom, s_top]``.
    drag : DragClosure, optional
        Per-layer drag closure. Leave ``None`` when a bulk closure is passed to the
        solver instead (theory.md §3.9).
    """

    name: str
    s_bottom: float
    s_top: float
    n_per_m2: float
    frontal_width: Width
    drag: DragClosure | None = None

    def __post_init__(self) -> None:
        """Validate the layer geometry."""
        if not (0.0 <= self.s_bottom < self.s_top and math.isfinite(self.s_top)):
            raise ValueError(
                f"layer {self.name!r}: need 0 <= s_bottom < s_top, "
                f"got [{self.s_bottom}, {self.s_top}]"
            )
        if not (math.isfinite(self.n_per_m2) and self.n_per_m2 > 0.0):
            raise ValueError(f"layer {self.name!r}: n_per_m2 must be positive")
        if not callable(self.frontal_width) and not (
            math.isfinite(self.frontal_width) and self.frontal_width > 0.0
        ):
            raise ValueError(f"layer {self.name!r}: frontal_width must be positive")

    @property
    def has_constant_width(self) -> bool:
        """True if the frontal width is a constant (analytic integral, theory.md Eq. 3.5)."""
        return not callable(self.frontal_width)

    def width_at(self, s: float) -> float:
        """Frontal width of one element at height s, in m (0 outside the layer)."""
        if not self.s_bottom <= s <= self.s_top:
            return 0.0
        if callable(self.frontal_width):
            return float(self.frontal_width(s))
        return float(self.frontal_width)

    def submerged_range(self, depth: float) -> tuple[float, float] | None:
        """Submerged extent ``(s_bottom, min(s_top, depth))``, or None if dry (m)."""
        top = min(self.s_top, depth)
        if top <= self.s_bottom:
            return None
        return (self.s_bottom, top)


@dataclass(frozen=True)
class Vegetation:
    """An ordered stack of :class:`Layer` objects describing one vegetation type.

    Parameters
    ----------
    layers : tuple of Layer
        Layers, normally ordered from the bed up. Names must be unique.
    name : str, optional
        Name of the vegetation type.
    source : str, optional
        Citation for the parameter values.
    notes : str, optional
        Free-text notes.
    """

    layers: tuple[Layer, ...]
    name: str = ""
    source: str = ""
    notes: str = ""

    def __post_init__(self) -> None:
        """Check that layer names are unique and convert lists to tuples."""
        object.__setattr__(self, "layers", tuple(self.layers))
        names = [layer.name for layer in self.layers]
        if len(set(names)) != len(names):
            raise ValueError(f"layer names must be unique, got {names}")

    def frontal_area_profile(self, s: npt.ArrayLike) -> npt.NDArray[np.float64]:
        """Frontal area per unit volume a(s) = sum_i N_i b_i(s), in m²/m³ (1/m).

        Parameters
        ----------
        s : array_like
            Heights above the bed, in m.

        Returns
        -------
        numpy.ndarray
            Frontal area per unit volume at each height, in 1/m.
        """
        s_arr = np.asarray(s, dtype=np.float64)
        flat = [
            sum(layer.n_per_m2 * layer.width_at(float(si)) for layer in self.layers)
            for si in s_arr.ravel()
        ]
        return np.asarray(flat, dtype=np.float64).reshape(s_arr.shape)

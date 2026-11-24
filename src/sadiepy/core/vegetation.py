"""Vegetation described as a vertical stack of layers (theory.md §3.2 and §5).

Heights are measured upward from the bed (theory.md §2.1). Besides the geometry
needed for drag dissipation (frontal width, density, vertical extent), a layer may
carry the solid cross-sectional area of its elements, from which the submerged
volume, volume fraction and characteristic lengths of theory.md §5.1-5.2 follow.
"""

import math
from collections.abc import Callable
from dataclasses import dataclass, field, replace

import numpy as np
import numpy.typing as npt
from scipy.integrate import quad

from sadiepy.core.drag import DragClosure, LengthScale

__all__ = ["LengthScale", "Layer", "Vegetation", "Width"]

Width = float | Callable[[float], float]
"""A constant in SI units, or a function of the height s above the bed (m)."""


_QUAD_EPSREL = 1e-12


def _integrate(profile: Callable[[float], float], s_b: float, s_t: float) -> float:
    """∫ profile(s) ds over [s_b, s_t] by ``scipy.integrate.quad``."""
    value, _ = quad(profile, s_b, s_t, epsabs=0.0, epsrel=_QUAD_EPSREL, limit=200)
    return float(value)


def _check_positive_width(name: str, label: str, value: Width | None) -> None:
    if value is None or callable(value):
        return
    if not (math.isfinite(value) and value > 0.0):
        raise ValueError(f"layer {name!r}: {label} must be positive")


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
    solid_area : float or callable, optional
        Solid (plan) cross-sectional area of one element, in m², constant or a
        function a(s) of height (theory.md §5.1), e.g. ``pi * d**2 / 4`` for a
        cylinder. Needed only for volumes, SVF and characteristic lengths.
    """

    name: str
    s_bottom: float
    s_top: float
    n_per_m2: float
    frontal_width: Width
    drag: DragClosure | None = None
    solid_area: Width | None = None

    def __post_init__(self) -> None:
        """Validate the layer geometry."""
        if not (0.0 <= self.s_bottom < self.s_top and math.isfinite(self.s_top)):
            raise ValueError(
                f"layer {self.name!r}: need 0 <= s_bottom < s_top, "
                f"got [{self.s_bottom}, {self.s_top}]"
            )
        if not (math.isfinite(self.n_per_m2) and self.n_per_m2 > 0.0):
            raise ValueError(f"layer {self.name!r}: n_per_m2 must be positive")
        _check_positive_width(self.name, "frontal_width", self.frontal_width)
        _check_positive_width(self.name, "solid_area", self.solid_area)

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

    def frontal_area(self, depth: float) -> float:
        """Submerged frontal area of one element, ∫ b ds, in m² (theory.md Eq. 5.1)."""
        span = self.submerged_range(depth)
        if span is None:
            return 0.0
        if callable(self.frontal_width):
            return _integrate(self.frontal_width, span[0], span[1])
        return float(self.frontal_width) * (span[1] - span[0])

    def solid_volume(self, depth: float) -> float:
        """Submerged solid volume of one element, ∫ a ds, in m³ (theory.md Eq. 5.1).

        Raises
        ------
        ValueError
            If the layer has no ``solid_area``.
        """
        if self.solid_area is None:
            raise ValueError(f"layer {self.name!r} has no solid_area; volumes are undefined")
        span = self.submerged_range(depth)
        if span is None:
            return 0.0
        if callable(self.solid_area):
            return _integrate(self.solid_area, span[0], span[1])
        return float(self.solid_area) * (span[1] - span[0])


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
    reference_diameter : float, optional
        The diameter d, in m, used for Re_d, KC_d and r_v (theory.md Eq. 5.3), e.g.
        the stem diameter at breast height (DBH) of a mangrove.
    derived_from : Vegetation, optional
        Set by :meth:`layer_averaged` and :meth:`depth_averaged`: the original
        vegetation, from which volumes and characteristic lengths are computed.
    derived_depth : float, optional
        Still-water depth, in m, for which an averaged geometry was derived.
    """

    layers: tuple[Layer, ...]
    name: str = ""
    source: str = ""
    notes: str = ""
    reference_diameter: float | None = None
    derived_from: "Vegetation | None" = field(default=None, repr=False)
    derived_depth: float | None = None

    def __post_init__(self) -> None:
        """Check that layer names are unique and convert lists to tuples."""
        object.__setattr__(self, "layers", tuple(self.layers))
        names = [layer.name for layer in self.layers]
        if len(set(names)) != len(names):
            raise ValueError(f"layer names must be unique, got {names}")
        d = self.reference_diameter
        if d is not None and not (math.isfinite(d) and d > 0.0):
            raise ValueError(f"reference_diameter must be positive, got {d}")

    @property
    def height(self) -> float:
        """Height of the tallest layer above the bed, h_v, in m."""
        return max(layer.s_top for layer in self.layers)

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

    # --- Volumes and characteristic lengths (theory.md §5.1-5.2) ----------------------

    @property
    def _base(self) -> "Vegetation":
        """The original vegetation (for averaged geometries), else self."""
        return self.derived_from._base if self.derived_from is not None else self

    def frontal_area(self, depth: float) -> float:
        """Submerged frontal area per m² of bed, A_f, in m²/m² (theory.md Eq. 5.1).

        For an averaged geometry this is the area of the original vegetation.
        """
        return sum(layer.n_per_m2 * layer.frontal_area(depth) for layer in self._base.layers)

    def submerged_volume(self, depth: float) -> float:
        """Submerged solid volume per m² of bed, V_s, in m³/m² (theory.md Eq. 5.1)."""
        return sum(layer.n_per_m2 * layer.solid_volume(depth) for layer in self._base.layers)

    def svf(self, depth: float) -> float:
        """Submerged volume fraction V_s / h, dimensionless (theory.md Eq. 5.2)."""
        return self.submerged_volume(depth) / depth

    def l_e(self, depth: float) -> float:
        """Effective vegetation length V_s / A_f, in m (theory.md Eq. 5.3)."""
        return self.submerged_volume(depth) / self.frontal_area(depth)

    def l_ew(self, depth: float) -> float:
        """Effective water length (h - V_s) / A_f, in m (theory.md Eq. 5.3)."""
        return (depth - self.submerged_volume(depth)) / self.frontal_area(depth)

    def hydraulic_radius(self, depth: float) -> float:
        """Vegetation hydraulic radius r_v = (pi/4) d (1 - SVF)/SVF, in m (Eq. 5.3)."""
        phi = self.svf(depth)
        return math.pi / 4.0 * self._diameter() * (1.0 - phi) / phi

    def length_scale(self, name: LengthScale, depth: float) -> float:
        """Value of a named characteristic length at a given depth, in m.

        Parameters
        ----------
        name : {"d", "l_e", "l_ew", "r_v"}
            Reference diameter, effective vegetation length, effective water
            length, or vegetation hydraulic radius (theory.md Eq. 5.3).
        depth : float
            Still-water depth, in m.
        """
        if name == "d":
            return self._diameter()
        if name == "l_e":
            return self.l_e(depth)
        if name == "l_ew":
            return self.l_ew(depth)
        if name == "r_v":
            return self.hydraulic_radius(depth)
        raise ValueError(f"unknown length scale {name!r}")

    def _diameter(self) -> float:
        d = self._base.reference_diameter
        if d is None:
            raise ValueError(f"vegetation {self.name!r} has no reference_diameter")
        return d

    # --- Averaged geometries (theory.md §5.3) -----------------------------------------

    def layer_averaged(self, depth: float) -> "Vegetation":
        """Replace each layer's width by its mean over the submerged part.

        This is the geometry of Zhang et al. (2023), Eqs. (19)-(20)
        (refs/Zhang_et_al_2023.pdf, p. 9); theory.md Eq. (5.4).

        Parameters
        ----------
        depth : float
            Still-water depth, in m.

        Returns
        -------
        Vegetation
            Constant-width layers, valid at this depth only.
        """
        layers = []
        for layer in self.layers:
            span = layer.submerged_range(depth)
            if span is None or layer.has_constant_width:
                layers.append(layer)
                continue
            mean_width = layer.frontal_area(depth) / (span[1] - span[0])
            layers.append(replace(layer, frontal_width=mean_width))
        return Vegetation(
            layers=tuple(layers),
            name=f"{self.name} [layer-averaged, h={depth:g} m]",
            source=self.source,
            notes=self.notes,
            reference_diameter=self.reference_diameter,
            derived_from=self,
            derived_depth=depth,
        )

    def depth_averaged(self, depth: float) -> "Vegetation":
        """One uniform layer on [0, min(h_v, h)] with frontal area A_f(h)/h per unit volume.

        This is the geometry of Zhang et al. (2023), Eq. (4) with
        D̄ = (1/h)∫A dz (refs/Zhang_et_al_2023.pdf, p. 3); theory.md Eq. (5.5).
        The layer has ``n_per_m2 = 1`` and a frontal width A_f(h)/h, so only their
        product is meaningful.

        Parameters
        ----------
        depth : float
            Still-water depth, in m.

        Returns
        -------
        Vegetation
            A single layer named ``"depth-averaged"``, valid at this depth only.
        """
        top = min(self.height, depth)
        layer = Layer(
            name="depth-averaged",
            s_bottom=0.0,
            s_top=top,
            n_per_m2=1.0,
            frontal_width=self.frontal_area(depth) / depth,
        )
        return Vegetation(
            layers=(layer,),
            name=f"{self.name} [depth-averaged, h={depth:g} m]",
            source=self.source,
            notes=self.notes,
            reference_diameter=self.reference_diameter,
            derived_from=self,
            derived_depth=depth,
        )

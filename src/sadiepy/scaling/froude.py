"""Froude similarity between a model at scale 1:lambda and its prototype (theory.md §4, §5.5).

``scale`` is the geometric ratio lambda = prototype length / model length (e.g. 10
for a 1:10 model). The same fluid density and viscosity are assumed in both.
"""

import math
from collections.abc import Callable
from dataclasses import replace
from typing import Literal

from sadiepy.core.vegetation import Layer, Vegetation, Width
from sadiepy.core.waves import RegularWave

__all__ = [
    "EXPONENTS",
    "Quantity",
    "rescale_reynolds",
    "scale_factor",
    "scale_vegetation",
    "to_model",
    "to_prototype",
    "vegetation_to_prototype",
    "wave_to_prototype",
]

Quantity = Literal[
    "length",
    "velocity",
    "time",
    "wavenumber",
    "beta",
    "force",
    "reynolds",
    "kc",
    "ursell",
    "elastic_modulus",
]

EXPONENTS: dict[Quantity, float] = {
    "length": 1.0,
    "velocity": 0.5,
    "time": 0.5,
    "wavenumber": -1.0,
    "beta": -1.0,
    "force": 3.0,
    "reynolds": 1.5,
    "kc": 0.0,
    "ursell": 0.0,
    "elastic_modulus": 1.0,
}
"""Exponent n such that prototype / model = lambda**n (theory.md §4 table)."""


def scale_factor(quantity: Quantity, scale: float) -> float:
    """Prototype/model ratio of a quantity, lambda**n (theory.md §4).

    Parameters
    ----------
    quantity : str
        One of the keys of :data:`EXPONENTS`.
    scale : float
        Geometric scale lambda = prototype length / model length. Must be positive.

    Returns
    -------
    float
        Ratio prototype / model (dimensionless).
    """
    if not scale > 0.0:
        raise ValueError(f"scale must be positive, got {scale}")
    return float(scale ** EXPONENTS[quantity])


def to_prototype(value: float, quantity: Quantity, scale: float) -> float:
    """Convert a model-scale value to prototype scale (theory.md §4).

    Parameters
    ----------
    value : float
        Model-scale value, in SI units.
    quantity : str
        Kind of quantity (a key of :data:`EXPONENTS`).
    scale : float
        Geometric scale lambda (prototype / model).

    Returns
    -------
    float
        Prototype-scale value, in SI units.

    Examples
    --------
    >>> to_prototype(1.4, "time", 100.0)
    14.0
    """
    return value * scale_factor(quantity, scale)


def to_model(value: float, quantity: Quantity, scale: float) -> float:
    """Convert a prototype-scale value to model scale (inverse of :func:`to_prototype`)."""
    return value / scale_factor(quantity, scale)


def wave_to_prototype(wave: RegularWave, scale: float) -> RegularWave:
    """Scale a model wave to the prototype: H and h by lambda, T by sqrt(lambda).

    Parameters
    ----------
    wave : RegularWave
        Model-scale wave.
    scale : float
        Geometric scale lambda (prototype / model).

    Returns
    -------
    RegularWave
        Prototype-scale wave (``ursell_max`` is carried over unchanged).
    """
    return RegularWave(
        H=to_prototype(wave.H, "length", scale),
        T=to_prototype(wave.T, "time", scale),
        depth=to_prototype(wave.depth, "length", scale),
        ursell_max=wave.ursell_max,
    )


def rescale_reynolds(re: float, from_scale: float, to_scale: float) -> float:
    """Convert a Reynolds number between Froude-similar scales (theory.md Eq. 5.6).

    Re_to = Re_from * (lambda_from / lambda_to)**1.5, as used by Zhang et al. (2023),
    §4.4 and Conclusions (refs/Zhang_et_al_2023.pdf, pp. 12-13), after Kelty et al.
    (2022).

    Parameters
    ----------
    re : float
        Reynolds number at scale 1:``from_scale`` (dimensionless).
    from_scale, to_scale : float
        Geometric scales lambda (prototype / model) of the two cases; 1 is prototype.

    Returns
    -------
    float
        Reynolds number of the Froude-similar flow at scale 1:``to_scale``.

    Examples
    --------
    >>> rescale_reynolds(1000.0, 10.0, 1.0)  # 1:10 model to prototype
    31622.776601683792
    """
    if not (from_scale > 0.0 and to_scale > 0.0):
        raise ValueError("scales must be positive")
    return float(re * (from_scale / to_scale) ** 1.5)


def _scaled_profile(profile: Width | None, factor: float, power: int) -> Width | None:
    """Profile of a geometrically scaled element: factor**power * p(s / factor)."""
    if profile is None:
        return None
    if callable(profile):
        f: Callable[[float], float] = profile
        return lambda s: factor**power * f(s / factor)
    return factor**power * float(profile)


def scale_vegetation(vegetation: Vegetation, factor: float) -> Vegetation:
    """Scale all lengths of a vegetation by ``factor`` (geometric similarity).

    Heights and widths scale by ``factor``, element densities (per m²) by
    ``factor**-2`` and solid cross-sections by ``factor**2``. Per-layer drag closures
    are kept unchanged.

    Parameters
    ----------
    vegetation : Vegetation
        Vegetation to scale (not an averaged geometry).
    factor : float
        Length ratio new / old, e.g. lambda to go from model to prototype.

    Returns
    -------
    Vegetation
    """
    if not (math.isfinite(factor) and factor > 0.0):
        raise ValueError(f"factor must be positive, got {factor}")
    if vegetation.derived_from is not None:
        raise ValueError("scale the original vegetation, not an averaged geometry")
    layers = tuple(
        Layer(
            name=layer.name,
            s_bottom=layer.s_bottom * factor,
            s_top=layer.s_top * factor,
            n_per_m2=layer.n_per_m2 / factor**2,
            frontal_width=_scaled_profile(layer.frontal_width, factor, 1),  # type: ignore[arg-type]
            drag=layer.drag,
            solid_area=_scaled_profile(layer.solid_area, factor, 2),
        )
        for layer in vegetation.layers
    )
    d = vegetation.reference_diameter
    return replace(vegetation, layers=layers, reference_diameter=None if d is None else d * factor)


def vegetation_to_prototype(vegetation: Vegetation, scale: float) -> Vegetation:
    """Scale a model-scale vegetation to the prototype (lengths times lambda)."""
    return scale_vegetation(vegetation, to_prototype(1.0, "length", scale))

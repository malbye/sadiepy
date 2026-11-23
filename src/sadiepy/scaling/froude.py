"""Froude similarity between a model at scale 1:lambda and its prototype (theory.md §4).

``scale`` is the geometric ratio lambda = prototype length / model length (e.g. 10
for a 1:10 model). The same fluid density and viscosity are assumed in both.
"""

from typing import Literal

from sadiepy.core.waves import RegularWave

__all__ = ["EXPONENTS", "Quantity", "scale_factor", "to_model", "to_prototype", "wave_to_prototype"]

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

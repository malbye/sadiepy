"""Generic vegetation presets."""

import math

from sadiepy.core.vegetation import Layer, Vegetation

__all__ = ["rigid_cylinders"]


def rigid_cylinders(d: float, n_per_m2: float, height: float) -> Vegetation:
    """Uniform array of rigid vertical cylinders standing on the bed.

    Parameters
    ----------
    d : float
        Cylinder diameter, in m (also the reference diameter for Re_d and KC_d).
    n_per_m2 : float
        Number of cylinders per m² of bed.
    height : float
        Cylinder height above the bed, in m.

    Returns
    -------
    Vegetation
        One layer ``"cylinders"`` on [0, height], frontal width d and solid
        cross-section pi d²/4.

    Examples
    --------
    >>> veg = rigid_cylinders(d=0.01, n_per_m2=400.0, height=0.3)
    >>> round(veg.svf(0.3), 5)
    0.03142
    """
    layer = Layer("cylinders", 0.0, height, n_per_m2, d, solid_area=math.pi * d**2 / 4.0)
    return Vegetation(
        layers=(layer,),
        name=f"rigid cylinders (d={d:g} m, N={n_per_m2:g}/m², height={height:g} m)",
        reference_diameter=d,
    )

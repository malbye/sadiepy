"""Artificial mangrove models of Zhang et al. (2023), Coastal Engineering 186, 104403.

doi:10.1016/j.coastaleng.2023.104403 (refs/Zhang_et_al_2023.pdf). Page numbers
below are the journal page numbers printed on the PDF.

Model scale 1:10. Each tree is 0.40 m high: 28 prop roots (0.10 m high), a stem
(DBH 0.01 m) and, in M2 and M3, a canopy 0.145 m wide and 0.17 m high on top
(Table 1, p. 4; Fig. 2, p. 5). The patch is 2.0 m long and 0.5 m wide with 36 trees,
i.e. 36 trees/m² (§3.2, pp. 5-6).

=====================  ==================  =====================================
Parameter              Model (1:10)        Source
=====================  ==================  =====================================
Tree density N         36 /m²              Table 1 (p. 4); Table 2 (p. 5)
Tree height h_v        0.40 m              Table 1; Fig. 2a
Stem diameter (DBH)    0.01 m, 0-0.40 m    Table 1; Fig. 3 caption (p. 6)
Roots                  28/tree, d 0.005 m, Table 1; §3.2 (p. 5); Table 2
                       0-0.10 m, 1008 /m²
Canopy                 W_c 0.145 m,        Table 1; Fig. 2a (10 + 13 + 17 cm)
                       h_c 0.17 m, on
                       0.23-0.40 m
Canopy volume V_c      M2 5e-5 m³,         Table 1; §3.2 (p. 6)
                       M3 1e-4 m³
=====================  ==================  =====================================

Geometry conventions (Fig. 3 caption, p. 6, and decisions of 2026-11-23):

* Canopy frontal area per tree is an inverted triangle, apex at the canopy base:
  0.5 W_cs h_cs with W_cs = W_c h_cs / h_c, so b(s) = f W_c (s - 0.23) / h_c.
  The factor f is 1 for M2 and 2 for M3, which reproduces Fig. 3a and 3b (the
  caption formula alone gives M2 and M3 the same area; spec erratum E7).
* Canopy solid volume is an inverted cone, (h_cs/h_c)³ V_c, so the solid
  cross-section is a(s) = 3 V_c (s - 0.23)² / h_c³.
* The stem spans the full tree height, as in the caption's DBH h_vs term.
"""

import math
from dataclasses import dataclass
from typing import Literal

from sadiepy.core.vegetation import Layer, Vegetation
from sadiepy.core.waves import RegularWave
from sadiepy.scaling.froude import vegetation_to_prototype

__all__ = [
    "CANOPY_AREA_FACTOR",
    "CANOPY_VOLUME",
    "FIG4_BETA",
    "HSVI_COEFF",
    "HSVI_WEIGHTS",
    "PSI",
    "MODEL_SCALE",
    "PATCH_LENGTH",
    "PATCH_WIDTH",
    "REPORTED_RANGES",
    "TEST_MATRIX",
    "Config",
    "ZhangTest",
    "beta_hsvf",
    "beta_hsvi",
    "hsvf",
    "hsvi",
    "svf_parts",
    "zhang2023",
]

Config = Literal["M1", "M2", "M3"]

MODEL_SCALE = 10.0
"""Geometric scale 1:10 (Table 1, p. 4)."""
PATCH_LENGTH = 2.0
"""Patch length along the flume, m, model scale (§3.1, p. 4)."""
PATCH_WIDTH = 0.5
"""Patch width, m, model scale (§3.1, p. 4)."""

N_TREES = 36.0  # trees per m² (Table 1, p. 4)
TREE_HEIGHT = 0.40  # h_v, m (Table 1)
DBH = 0.01  # stem diameter, m (Table 1)
ROOTS_PER_TREE = 28  # §3.2, p. 5
ROOT_DIAMETER = 0.005  # m (Table 1)
ROOT_HEIGHT = 0.10  # m (Table 1)
CANOPY_WIDTH = 0.145  # W_c, m (Table 1)
CANOPY_HEIGHT = 0.17  # h_c, m (Table 1)
CANOPY_BOTTOM = TREE_HEIGHT - CANOPY_HEIGHT  # 0.23 m (Fig. 2a: 10 + 13 cm)

CANOPY_VOLUME: dict[str, float] = {"M2": 5e-5, "M3": 1e-4}
"""Canopy volume per tree V_c, m³, model scale (Table 1, p. 4; §3.2, p. 6)."""
CANOPY_AREA_FACTOR: dict[str, float] = {"M2": 1.0, "M3": 2.0}
"""Multiplier on the triangular canopy frontal area (reproduces Fig. 3a; erratum E7)."""


def zhang2023(
    config: Config,
    scale: Literal["model", "prototype"] = "model",
    root_frontal_factor: float = 1.0,
) -> Vegetation:
    """Artificial mangrove model M1, M2 or M3 of Zhang et al. (2023).

    Parameters
    ----------
    config : {"M1", "M2", "M3"}
        M1 stem + root; M2 + sparse canopy; M3 + dense canopy (Table 2, p. 5).
    scale : {"model", "prototype"}, optional
        Model scale (1:10, as built) or prototype scale (lengths times 10).
    root_frontal_factor : float, optional
        Multiplier in (0, 1] on the summed frontal width of the 28 roots per tree;
        default 1.0, the root area as built and as plotted in Fig. 3a (p. 6). Zhang
        et al.'s beta and C_D analysis (Fig. 8, Eq. 20, p. 9; M1 mean C_D 2.63, p. 10)
        is reproduced with about 0.5; the paper does not state this, so it is an
        inferred analysis convention, not a measured property. Root solid volume is
        unaffected. See the reproduction notebook, section 7b.

    Returns
    -------
    Vegetation
        Layers ``"root"``, ``"stem"`` and, for M2 and M3, ``"canopy"``, with
        ``reference_diameter`` equal to the DBH.

    Examples
    --------
    >>> veg = zhang2023("M3")
    >>> round(veg.frontal_area(0.5) / 36, 5)  # A_p per tree at h = 0.5 m, m²
    0.04265
    """
    if config not in ("M1", "M2", "M3"):
        raise ValueError(f"config must be 'M1', 'M2' or 'M3', got {config!r}")
    if not 0.0 < root_frontal_factor <= 1.0:
        raise ValueError(f"root_frontal_factor must be in (0, 1], got {root_frontal_factor}")
    if scale not in ("model", "prototype"):
        raise ValueError(f"scale must be 'model' or 'prototype', got {scale!r}")
    layers = [
        Layer(
            "root",
            0.0,
            ROOT_HEIGHT,
            ROOTS_PER_TREE * N_TREES,
            ROOT_DIAMETER * root_frontal_factor,
            solid_area=math.pi * ROOT_DIAMETER**2 / 4.0,
        ),
        Layer("stem", 0.0, TREE_HEIGHT, N_TREES, DBH, solid_area=math.pi * DBH**2 / 4.0),
    ]
    if config in ("M2", "M3"):
        factor = CANOPY_AREA_FACTOR[config]
        volume = CANOPY_VOLUME[config]
        layers.append(
            Layer(
                "canopy",
                CANOPY_BOTTOM,
                TREE_HEIGHT,
                N_TREES,
                lambda s: factor * CANOPY_WIDTH * (s - CANOPY_BOTTOM) / CANOPY_HEIGHT,
                solid_area=lambda s: 3.0 * volume * (s - CANOPY_BOTTOM) ** 2 / CANOPY_HEIGHT**3,
            )
        )
    veg = Vegetation(
        layers=tuple(layers),
        name=f"Zhang et al. (2023) {config}, 1:10 model",
        source="Zhang et al. (2023), Coastal Eng. 186, 104403, Table 1 and Fig. 3",
        notes="Canopy frontal area doubled for M3 to match Fig. 3a (erratum E7)."
        + (
            f" Root frontal area x{root_frontal_factor:g} (inferred analysis convention)."
            if root_frontal_factor != 1.0
            else ""
        ),
        reference_diameter=DBH,
    )
    if scale == "prototype":
        proto = vegetation_to_prototype(veg, MODEL_SCALE)
        return Vegetation(
            layers=proto.layers,
            name=f"Zhang et al. (2023) {config}, prototype",
            source=veg.source,
            notes=veg.notes,
            reference_diameter=proto.reference_diameter,
        )
    return veg


@dataclass(frozen=True)
class ZhangTest:
    """One wave condition of the test matrix (model scale).

    Attributes
    ----------
    H : float
        Target incident wave height, m.
    T : float
        Wave period, s.
    h : float
        Still-water depth, m.
    configs : tuple of str
        Configurations tested, including the bare bed M0.
    """

    H: float
    T: float
    h: float
    configs: tuple[str, ...]


def _rows(
    h: float, configs: tuple[str, ...], spec: dict[float, tuple[float, ...]]
) -> list[ZhangTest]:
    return [ZhangTest(H, T, h, configs) for H, periods in spec.items() for T in periods]


_EMERGENT = ("M0", "M1", "M3")
_ALL = ("M0", "M1", "M2", "M3")

TEST_MATRIX: tuple[ZhangTest, ...] = tuple(
    _rows(0.3, _EMERGENT, {0.04: (1.0, 1.2, 1.4), 0.06: (1.2, 1.4), 0.08: (1.2, 1.4, 1.6)})
    + _rows(0.4, _ALL, {0.06: (1.0, 1.2, 1.4), 0.08: (1.2, 1.4, 1.6), 0.10: (1.4, 1.6, 2.0)})
    + _rows(
        0.5,
        _ALL,
        {
            0.04: (1.0,),
            0.06: (1.2,),
            0.08: (1.4, 1.6, 1.8),
            0.12: (1.4, 1.6, 1.8),
            0.16: (1.4, 1.8, 2.0),
        },
    )
)
"""Wave conditions of Zhang et al. (2023), Table 2 (p. 5): 28 conditions, 104 tests."""

FIG4_BETA: dict[tuple[float, float, float], dict[str, float]] = {
    (0.06, 1.4, 0.3): {"M1": 0.042, "M3": 0.411},
    (0.08, 1.4, 0.3): {"M1": 0.055, "M3": 0.429},
    (0.08, 1.6, 0.3): {"M1": 0.047, "M3": 0.290},
    (0.08, 1.4, 0.4): {"M1": 0.066, "M2": 0.310, "M3": 0.709},
    (0.08, 1.6, 0.4): {"M1": 0.050, "M2": 0.236, "M3": 0.533},
    (0.10, 1.6, 0.4): {"M1": 0.070, "M2": 0.306, "M3": 0.697},
    (0.08, 1.4, 0.5): {"M1": 0.019, "M2": 0.122, "M3": 0.300},
    (0.08, 1.6, 0.5): {"M1": 0.011, "M2": 0.096, "M3": 0.226},
    (0.16, 1.8, 0.5): {"M1": 0.042, "M2": 0.185, "M3": 0.286},
}
"""Measured damping factors beta, 1/m, keyed by (H, T, h) in m, s, m.

As printed in the panels of Zhang et al. (2023), Fig. 4 (p. 7): Eq. (3) fitted to the
gauge data after subtracting the bare-bed (M0) reduction (§4.1, p. 7).
"""

REPORTED_RANGES: dict[str, tuple[float, float]] = {
    "kh": (0.68, 2.08),
    "H/L": (0.0186, 0.0622),
    "H/h": (0.08, 0.32),
    "Ur": (0.73, 21.33),
    "Re_d": (1.28e3, 4.27e3),
    "KC_d (Table 3)": (6.05, 64.90),
    "KC_d (Conclusions)": (12.97, 77.39),
    "Re_rv (Table 3)": (34900.0, 407000.0),
}
"""Hydraulic ranges stated by Zhang et al. (2023): Conclusions (pp. 13-14) and
Table 3 (p. 12). The two KC_d ranges disagree (spec erratum E5)."""


# --- Empirical damping model, Zhang et al. (2023) Eqs. (12)-(17) ----------------------
# refs/Zhang_et_al_2023.pdf, pp. 8-9; theory.md §5.7. Model scale (1:10), B = 2 m only.

PSI: dict[str, tuple[float, float]] = {
    "M1": (199.94, 14.48),
    "M2": (833.82, 34.07),
    "M3": (1319.90, 54.23),
}
"""beta = psi * HSVF: (psi, standard error) in 1/m for each model (Eqs. 13-15, p. 9)."""
HSVI_WEIGHTS = (1.6, 1.0)
"""(c1, c2) weighting canopy and cylindrical SVF in HSVI (Eq. 16, p. 9)."""
HSVI_COEFF = (6.67e4, 2.20e3)
"""beta = coeff * HSVI: (coefficient, standard error) in 1/m (Eq. 17, p. 9; erratum E2)."""


def svf_parts(config: Config, depth: float) -> tuple[float, float]:
    """Submerged volume fractions of the cylindrical part and the canopy (theory.md §5.7).

    Parameters
    ----------
    config : {"M1", "M2", "M3"}
        Model configuration.
    depth : float
        Still-water depth, m (model scale).

    Returns
    -------
    tuple of float
        (SVF_cy, SVF_ca): stem + root and canopy solid volume per unit water volume,
        dimensionless; they sum to the SVF.
    """
    veg = zhang2023(config)
    canopy = sum(
        layer.n_per_m2 * layer.solid_volume(depth) for layer in veg.layers if layer.name == "canopy"
    )
    svf_ca = canopy / depth
    return veg.svf(depth) - svf_ca, svf_ca


def hsvf(config: Config, wave: RegularWave, patch_length: float = PATCH_LENGTH) -> float:
    """Hydraulic submerged volume fraction SVF (H/L)(B/L)(h_vs/h) (Zhang Eq. 12; theory Eq. 5.7).

    Parameters
    ----------
    config : {"M1", "M2", "M3"}
        Model configuration.
    wave : RegularWave
        Incident wave at model scale.
    patch_length : float, optional
        Patch length B, m; default 2.0.

    Returns
    -------
    float
        HSVF, dimensionless.
    """
    h = wave.depth
    h_vs = min(TREE_HEIGHT, h)
    return zhang2023(config).svf(h) * wave.steepness * (patch_length / wave.L) * (h_vs / h)


def hsvi(config: Config, wave: RegularWave, patch_length: float = PATCH_LENGTH) -> float:
    """Hydraulic submerged volume index HSVF (c1 SVF_ca + c2 SVF_cy) (Zhang Eq. 16; Eq. 5.8)."""
    svf_cy, svf_ca = svf_parts(config, wave.depth)
    c1, c2 = HSVI_WEIGHTS
    return hsvf(config, wave, patch_length) * (c1 * svf_ca + c2 * svf_cy)


def beta_hsvf(config: Config, wave: RegularWave) -> float:
    """Damping factor from the per-model fits beta = psi HSVF (Zhang Eqs. 13-15), 1/m.

    Empirical, for the 1:10 model with B = 2 m only (theory.md §5.7).
    """
    return PSI[config][0] * hsvf(config, wave)


def beta_hsvi(config: Config, wave: RegularWave) -> float:
    """Damping factor from the unified fit beta = 6.67e4 HSVI (Zhang Eq. 17), 1/m.

    Empirical, for the 1:10 model with B = 2 m only (theory.md §5.7).
    """
    return HSVI_COEFF[0] * hsvi(config, wave)

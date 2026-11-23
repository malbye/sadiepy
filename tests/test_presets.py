"""Presets: Zhang et al. (2023) mangrove models and generic cylinders."""

import math

import pytest

import sadiepy as sd
from sadiepy.presets.zhang2023 import FIG4_BETA, TEST_MATRIX


def a_p_caption(config: str, h: float) -> float:
    """A_p per tree from the Fig. 3 caption (refs/Zhang_et_al_2023.pdf, p. 6).

    A_p = N_root/N D_root h_rs + DBH h_vs + f 0.5 W_cs h_cs, with W_cs = W_c h_cs/h_c
    and f = 0, 1, 2 for M1, M2, M3 (decision on erratum E7).
    """
    f = {"M1": 0.0, "M2": 1.0, "M3": 2.0}[config]
    h_rs, h_vs = min(0.1, h), min(0.4, h)
    h_cs = max(0.0, min(h, 0.4) - 0.23)
    return 28 * 0.005 * h_rs + 0.01 * h_vs + f * 0.5 * (0.145 * h_cs / 0.17) * h_cs


def sv_caption(config: str, h: float) -> float:
    """Submerged volume per tree from the Fig. 3 caption (p. 6)."""
    v_c = {"M1": 0.0, "M2": 5e-5, "M3": 1e-4}[config]
    h_rs, h_vs = min(0.1, h), min(0.4, h)
    h_cs = max(0.0, min(h, 0.4) - 0.23)
    return (
        28 * math.pi * 0.005**2 / 4 * h_rs + math.pi * 0.01**2 / 4 * h_vs + (h_cs / 0.17) ** 3 * v_c
    )


@pytest.mark.parametrize("config", ["M1", "M2", "M3"])
@pytest.mark.parametrize("h", [0.2, 0.3, 0.35, 0.4, 0.45, 0.5])
def test_fig3_area_and_volume(config: str, h: float) -> None:
    veg = sd.presets.zhang2023(config)  # type: ignore[arg-type]
    assert veg.frontal_area(h) / 36 == pytest.approx(a_p_caption(config, h), rel=1e-12)
    assert veg.submerged_volume(h) / 36 == pytest.approx(sv_caption(config, h), rel=1e-12)
    # SVF = N SV / V with V = B W h (Fig. 3 caption): 36 trees on 1 m².
    assert veg.svf(h) == pytest.approx(36 * sv_caption(config, h) / (2.0 * 0.5 * h), rel=1e-12)


def test_svf_range_matches_fig10() -> None:
    # Zhang et al. (2023), Fig. 10 (p. 8): 0.0062 <= SVF <= 0.0168 over the tests.
    svfs = [
        sd.presets.zhang2023(c).svf(t.h)  # type: ignore[arg-type]
        for t in TEST_MATRIX
        for c in t.configs
        if c != "M0"
    ]
    assert round(min(svfs), 4) == 0.0062
    assert round(max(svfs), 4) == 0.0168


def test_test_matrix_counts() -> None:
    # Table 2 (p. 5) and §3.2 (p. 6): 8, 9 and 11 conditions; 104 tests in total.
    assert [sum(t.h == h for t in TEST_MATRIX) for h in (0.3, 0.4, 0.5)] == [8, 9, 11]
    assert sum(len(t.configs) for t in TEST_MATRIX) == 104
    assert all("M2" not in t.configs for t in TEST_MATRIX if t.h == 0.3)
    conditions = {(t.H, t.T, t.h) for t in TEST_MATRIX}
    assert set(FIG4_BETA) <= conditions


def test_prototype_scale() -> None:
    m, p = sd.presets.zhang2023("M3"), sd.presets.zhang2023("M3", scale="prototype")
    # Table 1 (p. 4), prototype column: N 0.36/m², root density 10.08/m², DBH 0.1 m.
    assert p.layers[1].n_per_m2 == pytest.approx(0.36)
    assert p.layers[0].n_per_m2 == pytest.approx(10.08)
    assert p.reference_diameter == pytest.approx(0.1)
    assert p.height == pytest.approx(4.0)
    # Prototype canopy volume 1e-1 m³ (Table 1): submerged volume of the canopy layer.
    assert p.layers[2].solid_volume(10.0) == pytest.approx(0.1, rel=1e-10)
    assert p.svf(5.0) == pytest.approx(m.svf(0.5), rel=1e-12)


def test_bad_arguments() -> None:
    with pytest.raises(ValueError, match="config"):
        sd.presets.zhang2023("M4")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="scale"):
        sd.presets.zhang2023("M1", scale="full")  # type: ignore[arg-type]


def test_rigid_cylinders() -> None:
    veg = sd.presets.rigid_cylinders(d=0.02, n_per_m2=100.0, height=0.5)
    assert veg.l_e(0.6) == pytest.approx(math.pi * 0.02 / 4)
    assert veg.frontal_area(0.6) == pytest.approx(100 * 0.02 * 0.5)

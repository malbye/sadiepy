"""Published drag closures, Re rescaling and closure evaluation in the solver (M2)."""

import math
from collections.abc import Callable

import numpy as np
import pytest

from sadiepy import Layer, RegularWave, Vegetation, Water, propagate
from sadiepy.core import drag
from sadiepy.scaling.froude import (
    rescale_reynolds,
    scale_vegetation,
    vegetation_to_prototype,
    wave_to_prototype,
)

# Independent transcription of Zhang et al. (2023), Table 3 (refs/Zhang_et_al_2023.pdf,
# p. 12): name -> (C_D formula, variable, length scale, (min, max), source scale).
TABLE3: dict[str, tuple[Callable[[float], float], str, str | None, tuple[float, float], float]] = {
    "maza2017_re": (lambda x: 9.49 * x**-0.31, "Re", "d", (500, 1700), 12),
    "he2019_kc": (lambda x: 18.025 * math.exp(-0.043 * x), "KC", "d", (10, 37), 10),
    "wang2022b_kc": (lambda x: 0.42 + (0.77 / x) ** 0.41, "KC", "r_v", (0.01, 0.28), 10),
    "chang2022_kc": (lambda x: 1.15 + (14.70 / x) ** 1.065, "KC", "d", (2.5, 20), 7),
    "chang2022_re": (lambda x: 1.15 + (1015.21 / x) ** 1.044, "Re", "d", (180, 1300), 7),
    "bryant2022_kc": (lambda x: 0.86 + (1.24 / x) ** 0.78, "KC", "d", (0, 12), 2.1),
    "kelty2022_re": (lambda x: 0.60 + (300000 / x) ** 1.0, "Re", "d", (4900, 119000), 1),
    "zhang2023_kc": (lambda x: 0.60 + (101.99 / x) ** 1.22, "KC", "d", (6.05, 64.90), 10),
    "zhang2023_rerv": (
        lambda x: 0.60 + (420000 / x) ** 1.41,
        "Re",
        "r_v",
        (34900, 407000),
        10,
    ),
    "zhang2023_ur": (lambda x: 0.60 + (80.00 / x) ** 0.75, "Ur", None, (0.73, 21.33), 10),
}


def test_registry_has_table3_in_order() -> None:
    assert drag.names() == list(TABLE3)


@pytest.mark.parametrize("name", list(TABLE3))
def test_closure_matches_table3(name: str) -> None:
    formula, variable, length, valid, scale = TABLE3[name]
    c = drag.get(name)
    assert c.variable == variable
    assert c.length_scale == length
    assert c.valid_range == valid
    assert c.scale == scale
    assert "Table 3" in c.citation
    lo, hi = valid
    for x in np.linspace(max(lo, 1e-3 * hi), hi, 7):
        assert c(float(x)) == pytest.approx(formula(float(x)), rel=1e-14)
    # Every published fit decreases with its variable over its range.
    xs = np.linspace(max(lo, 1e-3 * hi), hi, 50)
    cds = [c(float(x)) for x in xs]
    assert all(a > b for a, b in zip(cds, cds[1:], strict=False))


def test_unknown_names_raise() -> None:
    with pytest.raises(KeyError, match="available"):
        drag.get("nope")
    with pytest.raises(ValueError, match="unknown length scale"):
        drag.custom("x", lambda v: 1.0, "KC", length_scale="diameter")  # type: ignore[arg-type]


# --- Re rescaling (theory.md Eq. 5.6) -------------------------------------------------


def test_rescale_reynolds() -> None:
    assert rescale_reynolds(1000.0, 10.0, 1.0) == pytest.approx(1000.0 * 10**1.5)
    assert rescale_reynolds(1000.0, 7.0, 7.0) == 1000.0
    re = rescale_reynolds(rescale_reynolds(1234.0, 10.0, 2.1), 2.1, 10.0)
    assert re == pytest.approx(1234.0, rel=1e-14)


# --- Vegetation scaling ---------------------------------------------------------------


def _tree() -> Vegetation:
    return Vegetation(
        layers=(
            Layer("root", 0.0, 0.1, 1008.0, 0.005, solid_area=math.pi * 0.005**2 / 4),
            Layer("stem", 0.0, 0.4, 36.0, 0.01, solid_area=math.pi * 0.01**2 / 4),
            Layer(
                "canopy",
                0.23,
                0.40,
                36.0,
                lambda s: 0.145 * (s - 0.23) / 0.17,
                solid_area=lambda s: 3e-4 * (s - 0.23) ** 2 / 0.17**3,
            ),
        ),
        reference_diameter=0.01,
    )


def test_scaled_vegetation_geometry() -> None:
    lam, h = 10.0, 0.5
    m = _tree()
    p = vegetation_to_prototype(m, lam)
    assert p.layers[2].s_bottom == pytest.approx(2.3)
    assert p.layers[0].n_per_m2 == pytest.approx(10.08)
    assert p.reference_diameter == pytest.approx(0.1)
    assert p.svf(h * lam) == pytest.approx(m.svf(h), rel=1e-12)
    assert p.l_e(h * lam) == pytest.approx(lam * m.l_e(h), rel=1e-12)
    assert p.hydraulic_radius(h * lam) == pytest.approx(lam * m.hydraulic_radius(h), rel=1e-12)
    assert p.frontal_area(h * lam) == pytest.approx(m.frontal_area(h), rel=1e-12)  # m²/m²
    with pytest.raises(ValueError, match="averaged"):
        scale_vegetation(m.depth_averaged(h), lam)


_IGNORE_EXTRAP = pytest.mark.filterwarnings("ignore::sadiepy.ExtrapolationWarning")


@_IGNORE_EXTRAP
@pytest.mark.parametrize("name", ["zhang2023_kc", "zhang2023_ur", "zhang2023_rerv"])
def test_closures_are_froude_consistent_with_scale_conversion(name: str) -> None:
    # KC and Ur are Froude-invariant; Re is too once converted to the source scale
    # (theory.md §5.5). With equal C_D at both scales, beta scales as 1/lambda.
    lam = 10.0
    wave_m = RegularWave(0.08, 1.4, 0.5)
    wave_p = wave_to_prototype(wave_m, lam)
    veg_m = _tree()
    veg_p = vegetation_to_prototype(veg_m, lam)
    kw = {"drag": name, "cd_evaluation": "incident"}
    r_m = propagate(wave_m, veg_m, 2.0, scale=lam, **kw)  # type: ignore[arg-type]
    r_p = propagate(wave_p, veg_p, 2.0 * lam, scale=1.0, **kw)  # type: ignore[arg-type]
    assert r_p.drag_coefficients["stem"][0] == pytest.approx(
        r_m.drag_coefficients["stem"][0], rel=1e-10
    )
    assert r_p.beta == pytest.approx(r_m.beta / lam, rel=1e-10)


@_IGNORE_EXTRAP
def test_without_scale_re_closure_sees_prototype_re() -> None:
    lam = 10.0
    wave_p = wave_to_prototype(RegularWave(0.08, 1.4, 0.5), lam)
    veg_p = vegetation_to_prototype(_tree(), lam)
    kw = {"drag": "zhang2023_rerv", "cd_evaluation": "incident"}
    conv = propagate(wave_p, veg_p, 20.0, scale=1.0, **kw)  # type: ignore[arg-type]
    raw = propagate(wave_p, veg_p, 20.0, **kw)  # type: ignore[arg-type]
    key = "Re[zhang2023_rerv]"
    re_conv = conv.dimensionless[key]["range"][0]  # type: ignore[index]
    re_raw = raw.dimensionless[key]["range"][0]  # type: ignore[index]
    assert re_raw == pytest.approx(re_conv * lam**1.5, rel=1e-12)
    assert conv.dimensionless[key]["re_scale_factor"] == pytest.approx(lam**-1.5)  # type: ignore[index]
    assert raw.dimensionless[key]["re_scale_factor"] == 1.0  # type: ignore[index]


# --- Solver options -------------------------------------------------------------------


def test_incident_evaluation_equals_constant_cd_at_h0() -> None:
    wave, veg, water = RegularWave(0.08, 1.4, 0.4), _tree(), Water()
    res = propagate(wave, veg, 2.0, drag="zhang2023_kc", cd_evaluation="incident")
    u_bar = wave.omega * wave.a / wave.kh
    kc = u_bar * wave.T / 0.01
    cd = 0.60 + (101.99 / kc) ** 1.22
    ref = propagate(wave, veg, 2.0, water=water, drag=cd)
    assert res.method == "analytic"
    assert res.cd_evaluation == "incident"
    assert res.beta == pytest.approx(ref.beta, rel=1e-13)
    assert res.drag_coefficients["canopy"] == pytest.approx((cd, cd), rel=1e-13)
    assert res.drag_used["root"] == "zhang2023_kc"
    assert res.dimensionless["KC[zhang2023_kc]"]["range"] == pytest.approx((kc, kc))  # type: ignore[index]


def test_local_evaluation_damps_more_for_decreasing_closure() -> None:
    # C_D falls as KC rises; KC falls as H decays, so C_D(H(x)) >= C_D(H0) and
    # local evaluation damps at least as much as incident evaluation.
    wave, veg = RegularWave(0.08, 1.4, 0.4), _tree()
    inc = propagate(wave, veg, 2.0, drag="zhang2023_kc", cd_evaluation="incident")
    loc = propagate(wave, veg, 2.0, drag="zhang2023_kc")
    assert loc.method == "ode"
    assert loc.beta > inc.beta
    lo, hi = loc.drag_coefficients["stem"]
    assert lo == pytest.approx(inc.drag_coefficients["stem"][0], rel=1e-9)
    assert hi > lo


def test_named_length_resolved_from_vegetation() -> None:
    wave, veg = RegularWave(0.08, 1.4, 0.5), _tree()
    res = propagate(wave, veg, 2.0, drag="zhang2023_rerv", cd_evaluation="incident")
    entry = res.dimensionless["Re[zhang2023_rerv]"]
    assert entry["length_scale"] == "r_v"  # type: ignore[index]
    assert entry["length_m"] == pytest.approx(veg.hydraulic_radius(0.5), rel=1e-14)  # type: ignore[index]


def test_averaged_geometry_depth_must_match_wave() -> None:
    veg = _tree().depth_averaged(0.4)
    with pytest.raises(ValueError, match="averaged for h"):
        propagate(RegularWave(0.08, 1.4, 0.5), veg, 2.0, drag=1.0)
    with pytest.raises(ValueError, match="cd_evaluation"):
        propagate(RegularWave(0.08, 1.4, 0.4), veg, 2.0, drag=1.0, cd_evaluation="mid")  # type: ignore[arg-type]

"""Tests for core.vegetation, core.drag and scaling.froude (spec M1 test 8)."""

import math

import numpy as np
import pytest

from sadiepy import Layer, RegularWave, Vegetation, Water, propagate
from sadiepy.core import drag
from sadiepy.core.waves import keulegan_carpenter, reynolds
from sadiepy.scaling.froude import scale_factor, to_model, to_prototype, wave_to_prototype

# --- Vegetation -------------------------------------------------------------------


def test_frontal_area_profile_sums_layers() -> None:
    veg = Vegetation(
        layers=(
            Layer("root", 0.0, 0.1, 1008.0, 0.005),
            Layer("stem", 0.0, 0.4, 36.0, 0.01),
        )
    )
    a = veg.frontal_area_profile([0.05, 0.2, 0.5])
    np.testing.assert_allclose(a, [1008 * 0.005 + 36 * 0.01, 36 * 0.01, 0.0], rtol=1e-15)


def test_submerged_range() -> None:
    layer = Layer("canopy", 0.23, 0.40, 36.0, 0.1)
    assert layer.submerged_range(0.5) == (0.23, 0.40)
    assert layer.submerged_range(0.3) == (0.23, 0.3)
    assert layer.submerged_range(0.2) is None


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"s_bottom": 0.4, "s_top": 0.1}, "s_bottom < s_top"),
        ({"s_bottom": -0.1}, "s_bottom < s_top"),
        ({"n_per_m2": 0.0}, "n_per_m2"),
        ({"frontal_width": -0.01}, "frontal_width"),
    ],
)
def test_layer_validation(kwargs: dict[str, float], match: str) -> None:
    args: dict[str, object] = {
        "name": "x",
        "s_bottom": 0.0,
        "s_top": 0.4,
        "n_per_m2": 36.0,
        "frontal_width": 0.01,
    }
    args.update(kwargs)
    with pytest.raises(ValueError, match=match):
        Layer(**args)  # type: ignore[arg-type]


def test_duplicate_layer_names_raise() -> None:
    layer = Layer("stem", 0.0, 0.4, 36.0, 0.01)
    with pytest.raises(ValueError, match="unique"):
        Vegetation(layers=(layer, layer))


# --- Drag -------------------------------------------------------------------------


def test_constant_closure() -> None:
    c = drag.constant(1.2)
    assert c() == 1.2
    assert c.is_constant
    assert c.name == "constant(1.2)"
    with pytest.raises(ValueError):
        drag.constant(0.0)


def test_re_closure_needs_length_scale() -> None:
    with pytest.raises(ValueError, match="length_scale"):
        drag.custom("bad", lambda re: 1.0, "Re")


# --- Test 8: Froude invariance ----------------------------------------------------


@pytest.mark.parametrize("lam", [2.0, 10.0, 25.0])
def test_froude_invariance(lam: float) -> None:
    model = RegularWave(H=0.08, T=1.4, depth=0.4)
    proto = wave_to_prototype(model, lam)
    d_model, d_proto = 0.01, to_prototype(0.01, "length", lam)
    nu = Water().nu  # same fluid at both scales

    assert proto.kh == pytest.approx(model.kh, rel=1e-12)
    assert proto.ursell == pytest.approx(model.ursell, rel=1e-12)
    kc_m = keulegan_carpenter(model.u_bar, model.T, d_model)
    kc_p = keulegan_carpenter(proto.u_bar, proto.T, d_proto)
    assert kc_p == pytest.approx(kc_m, rel=1e-12)
    re_m = reynolds(model.u_bar, d_model, nu)
    re_p = reynolds(proto.u_bar, d_proto, nu)
    assert re_p / re_m == pytest.approx(lam**1.5, rel=1e-12)
    assert re_p / re_m == pytest.approx(scale_factor("reynolds", lam), rel=1e-12)


def test_beta_scales_inversely_with_length() -> None:
    # Same C_D at both scales: beta (1/m) should scale as 1/lambda (theory.md §4).
    lam = 10.0
    model = RegularWave(H=0.08, T=1.4, depth=0.4)
    veg_m = Vegetation(layers=(Layer("stem", 0.0, 0.3, 36.0, 0.01),))
    veg_p = Vegetation(layers=(Layer("stem", 0.0, 0.3 * lam, 36.0 / lam**2, 0.01 * lam),))
    b_m = propagate(model, veg_m, 2.0, drag=1.0).beta
    b_p = propagate(wave_to_prototype(model, lam), veg_p, 2.0 * lam, drag=1.0).beta
    assert b_p == pytest.approx(b_m / lam, rel=1e-12)


def test_to_model_inverts_to_prototype() -> None:
    assert to_model(to_prototype(1.4, "time", 10.0), "time", 10.0) == pytest.approx(1.4)
    assert to_prototype(2.0, "kc", 10.0) == 2.0
    assert math.isclose(to_prototype(1.0, "force", 10.0), 1000.0)

"""Volumes, characteristic lengths and averaged geometries (theory.md §5.1-5.3)."""

import math

import pytest

from sadiepy import Layer, RegularWave, Vegetation, propagate

D, N, HV = 0.01, 36.0, 0.4
W_C, H_C, S0 = 0.145, 0.17, 0.23  # triangle canopy, apex at S0 (Zhang-like numbers)
V_C = 1e-4


def _sinh_s(x: float) -> float:
    return math.sinh(x) ** 3 + 3.0 * math.sinh(x)


def cylinders(top: float = HV) -> Vegetation:
    return Vegetation(
        layers=(Layer("stem", 0.0, top, N, D, solid_area=math.pi * D**2 / 4),),
        reference_diameter=D,
    )


def tree() -> Vegetation:
    """Roots + stem + inverted-triangle canopy with a cone-shaped solid volume."""
    return Vegetation(
        layers=(
            Layer("root", 0.0, 0.1, 1008.0, 0.005, solid_area=math.pi * 0.005**2 / 4),
            Layer("stem", 0.0, HV, N, D, solid_area=math.pi * D**2 / 4),
            Layer(
                "canopy",
                S0,
                S0 + H_C,
                N,
                lambda s: W_C * (s - S0) / H_C,
                solid_area=lambda s: 3.0 * V_C * (s - S0) ** 2 / H_C**3,
            ),
        ),
        reference_diameter=D,
    )


# --- Volumes and lengths (Eqs. 5.1-5.3) -----------------------------------------------


def test_cylinder_volumes_and_lengths() -> None:
    veg, h = cylinders(), 0.5
    a_f = N * D * HV
    v_s = N * math.pi * D**2 / 4 * HV
    assert veg.frontal_area(h) == pytest.approx(a_f, rel=1e-14)
    assert veg.submerged_volume(h) == pytest.approx(v_s, rel=1e-14)
    assert veg.svf(h) == pytest.approx(v_s / h, rel=1e-14)
    # For cylinders l_e = (pi d²/4) / d = pi d / 4 at any depth.
    assert veg.l_e(h) == pytest.approx(math.pi * D / 4, rel=1e-14)
    assert veg.l_e(0.3) == pytest.approx(math.pi * D / 4, rel=1e-14)
    assert veg.l_ew(h) == pytest.approx((h - v_s) / a_f, rel=1e-14)
    phi = v_s / h
    assert veg.hydraulic_radius(h) == pytest.approx(math.pi / 4 * D * (1 - phi) / phi, rel=1e-14)
    assert veg.length_scale("d", h) == D
    assert veg.length_scale("r_v", h) == veg.hydraulic_radius(h)


@pytest.mark.parametrize("h", [0.3, 0.35, 0.4, 0.5])
def test_triangle_canopy_area_and_cone_volume(h: float) -> None:
    canopy = tree().layers[2]
    h_cs = max(0.0, min(h, S0 + H_C) - S0)
    # Triangle with apex at the bottom: area 0.5 * W_cs * h_cs, W_cs = W_c h_cs / h_c
    # (Zhang et al. 2023, Fig. 3 caption, refs/Zhang_et_al_2023.pdf p. 6).
    assert canopy.frontal_area(h) == pytest.approx(0.5 * (W_C * h_cs / H_C) * h_cs, rel=1e-12)
    # Submerged cone volume (h_cs/h_c)³ V_c (same caption).
    assert canopy.solid_volume(h) == pytest.approx((h_cs / H_C) ** 3 * V_C, rel=1e-12)


def test_dry_layer_contributes_nothing() -> None:
    canopy = tree().layers[2]
    assert canopy.frontal_area(0.2) == 0.0
    assert canopy.solid_volume(0.2) == 0.0


def test_missing_inputs_raise() -> None:
    veg = Vegetation(layers=(Layer("stem", 0.0, 0.4, N, D),))
    with pytest.raises(ValueError, match="solid_area"):
        veg.svf(0.5)
    with pytest.raises(ValueError, match="reference_diameter"):
        veg.length_scale("d", 0.5)
    with pytest.raises(ValueError, match="reference_diameter"):
        Vegetation(layers=veg.layers, reference_diameter=0.0)
    with pytest.raises(ValueError, match="solid_area"):
        Layer("x", 0.0, 1.0, 1.0, 0.01, solid_area=-1.0)


# --- Layer-averaged geometry (Eq. 5.4 = Zhang Eqs. 19-20) -----------------------------


@pytest.mark.parametrize(("h", "period"), [(0.3, 1.4), (0.4, 1.6), (0.5, 2.0)])
def test_layer_averaged_ratios_match_zhang_eqs_19_20(h: float, period: float) -> None:
    wave = RegularWave(0.08, period, h)
    k = wave.k
    veg = tree()
    res = propagate(wave, veg.layer_averaged(h), 2.0, drag=1.0)
    f = res.layer_fractions
    # Zhang et al. (2023), Eqs. (19)-(20), refs/Zhang_et_al_2023.pdf p. 9, written out
    # independently: D̄ are the layer-depth-averaged projected widths.
    h_vs, h_rs = min(HV, h), min(0.1, h)
    h_ct, h_cb = min(S0 + H_C, h), S0
    d_canopy = 0.5 * W_C * (h_ct - h_cb) / H_C  # mean width of the submerged triangle
    beta_rc = d_canopy * (_sinh_s(k * h_ct) - _sinh_s(k * h_cb)) / (D * _sinh_s(k * h_vs))
    beta_rr = 1008.0 * 0.005 * _sinh_s(k * h_rs) / (N * D * _sinh_s(k * h_vs))
    assert f["canopy"] / f["stem"] == pytest.approx(beta_rc, rel=1e-10)
    assert f["root"] / f["stem"] == pytest.approx(beta_rr, rel=1e-10)


@pytest.mark.parametrize("h", [0.3, 0.4, 0.5])
def test_resolved_canopy_dissipates_at_least_layer_averaged(h: float) -> None:
    # Chebyshev's sum inequality (theory.md §5.3): b(s) and cosh³(ks) both increase.
    wave = RegularWave(0.08, 1.4, h)
    veg = tree()
    resolved = propagate(wave, veg, 2.0, drag=1.0).beta
    averaged = propagate(wave, veg.layer_averaged(h), 2.0, drag=1.0).beta
    assert resolved >= averaged


# --- Depth-averaged geometry (Eq. 5.5 = Zhang Eq. 4 with D̄) ---------------------------


@pytest.mark.parametrize(("h", "period"), [(0.3, 1.4), (0.4, 1.2), (0.5, 1.8)])
def test_depth_averaged_matches_zhang_eq4_with_dbar(h: float, period: float) -> None:
    wave = RegularWave(0.08, period, h)
    veg = tree()
    k, kh = wave.k, wave.kh
    # Zhang et al. (2023), Eq. (4), refs/Zhang_et_al_2023.pdf p. 3, with
    # D̄ = (1/h) ∫ A dz per tree and N trees/m²: N D̄ = A_f(h) / h.
    a_p = (
        28 * 0.005 * min(0.1, h)
        + D * min(HV, h)
        + 0.5 * W_C * max(0.0, min(h, S0 + H_C) - S0) ** 2 / H_C
    )
    n_dbar = N * a_p / h
    h_vs = min(S0 + H_C, h)
    beta_eq4 = (
        4.0
        / (9.0 * math.pi)
        * 1.0
        * n_dbar
        * wave.H
        * k
        * _sinh_s(k * h_vs)
        / ((math.sinh(2 * kh) + 2 * kh) * math.sinh(kh))
    )
    beta = propagate(wave, veg.depth_averaged(h), 2.0, drag=1.0).beta
    assert beta == pytest.approx(beta_eq4, rel=1e-10)


def test_averaged_geometry_keeps_characteristic_lengths_of_original() -> None:
    veg, h = tree(), 0.5
    for derived in (veg.depth_averaged(h), veg.layer_averaged(h)):
        assert derived.derived_depth == h
        assert derived.svf(h) == pytest.approx(veg.svf(h), rel=1e-14)
        assert derived.length_scale("r_v", h) == pytest.approx(veg.hydraulic_radius(h))
        assert derived.frontal_area(h) == pytest.approx(veg.frontal_area(h))


def test_depth_averaged_equals_layered_for_emergent_uniform_cylinders() -> None:
    # One uniform emergent layer is already depth-uniform: both geometries agree.
    wave = RegularWave(0.06, 1.2, 0.3)
    veg = cylinders(top=0.3)
    b1 = propagate(wave, veg, 2.0, drag=1.3).beta
    b2 = propagate(wave, veg.depth_averaged(0.3), 2.0, drag=1.3).beta
    assert b2 == pytest.approx(b1, rel=1e-12)

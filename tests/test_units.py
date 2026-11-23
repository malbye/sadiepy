"""Tests for sadiepy.units (theory.md §1)."""

import warnings

import pytest

from sadiepy import ValidityWarning, Water
from sadiepy.units import pure_water_dynamic_viscosity

# Hand-calculated from theory.md Eqs. (1.1), (1.4), (1.5) at t = 20 °C, S = 0
# (theory.md §1.5).
RHO_20 = 999.9 + 0.4068 - 2.4648 + 0.18088 - 0.0074512  # = 998.0154288 kg/m³
MU_20 = 4.2844e-5 + 1.0 / (0.157 * 84.993**2 - 91.296)  # = 1.00176187e-3 Pa·s


def test_fresh_water_20c_matches_hand_calculation() -> None:
    w = Water(temperature=20.0, salinity=0.0)
    assert w.rho == pytest.approx(RHO_20, rel=1e-12)
    assert w.mu == pytest.approx(MU_20, rel=1e-12)
    assert w.nu == pytest.approx(MU_20 / RHO_20, rel=1e-12)


# Seawater at t = 20 °C, 35 g/kg (S = 0.035 kg/kg), term by term from Sharqawy et al.
# (2010) Eq. (8) p. 357 and Eqs. (22)-(23) p. 364, coefficients re-typed from the PDF.
_S, _T = 0.035, 20.0
RHO_SEA_20 = RHO_20 + (
    8.020e2 * _S
    - 2.001 * _S * _T
    + 1.677e-2 * _S * _T**2
    - 3.060e-5 * _S * _T**3
    - 1.613e-5 * _S**2 * _T**2
)  # = 998.0154288 + 26.8955... kg/m³
_A = 1.541 + 1.998e-2 * _T - 9.52e-5 * _T**2  # = 1.90252
_B = 7.974 - 7.561e-2 * _T + 4.724e-4 * _T**2  # = 6.65076
MU_SEA_20 = MU_20 * (1.0 + _A * _S + _B * _S**2)


def test_seawater_20c_35gkg_matches_hand_calculation() -> None:
    w = Water(temperature=20.0, salinity=35.0)
    assert w.rho == pytest.approx(RHO_SEA_20, rel=1e-12)
    assert w.mu == pytest.approx(MU_SEA_20, rel=1e-12)
    # Pinned values, so a change in any coefficient digit is caught.
    assert w.rho == pytest.approx(1024.911, abs=5e-4)
    assert w.mu == pytest.approx(1.07663e-3, rel=1e-5)


def test_fresh_water_20c_within_stated_accuracy_of_reference_values() -> None:
    w = Water(temperature=20.0, salinity=0.0)
    # IAPWS-95 pure water at 20 °C, 0.1 MPa: 998.21 kg/m³ (theory.md §1.5);
    # Eq. (1.1) accuracy ±0.1 %.
    assert w.rho == pytest.approx(998.21, rel=1e-3)
    # 1.002e-3 Pa·s at 20 °C (Swindells et al., quoted in Sharqawy et al. 2010, p. 364);
    # Eq. (1.2) accuracy ±1.5 %.
    assert w.mu == pytest.approx(1.002e-3, rel=1.5e-2)


def test_defaults_are_fresh_water_at_20c() -> None:
    assert Water() == Water(temperature=20.0, salinity=0.0)


def test_salinity_is_in_grams_per_kilogram() -> None:
    # 35 g/kg must enter the correlation as 0.035 kg/kg: the density increase over
    # fresh water is then of order 800 * 0.035 ≈ 28 kg/m³, not ~28 000.
    increase = Water(20.0, 35.0).rho - Water(20.0, 0.0).rho
    assert 20.0 < increase < 35.0


def test_pure_water_limit_of_viscosity() -> None:
    assert Water(15.0, 0.0).mu == pure_water_dynamic_viscosity(15.0)


@pytest.mark.parametrize("t", [5.0, 20.0, 40.0])
def test_density_and_viscosity_increase_with_salinity(t: float) -> None:
    fresh, sea = Water(t, 0.0), Water(t, 35.0)
    assert sea.rho > fresh.rho
    assert sea.mu > fresh.mu


def test_viscosity_decreases_with_temperature() -> None:
    mus = [Water(t, 35.0).mu for t in (5.0, 20.0, 40.0, 80.0)]
    assert mus == sorted(mus, reverse=True)


@pytest.mark.parametrize(("t", "s"), [(0.0, 0.0), (20.0, 35.0), (180.0, 150.0)])
def test_no_warning_inside_validity_range(t: float, s: float) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        Water(t, s)


@pytest.mark.parametrize(
    ("t", "s", "match"),
    [
        (-2.0, 0.0, "temperature"),
        (190.0, 0.0, "temperature"),
        (20.0, 155.0, "150 g/kg"),
        (20.0, 170.0, "160 g/kg"),
    ],
)
def test_warns_outside_validity_range_but_returns_value(t: float, s: float, match: str) -> None:
    with pytest.warns(ValidityWarning, match=match):
        w = Water(t, s)
    assert w.rho > 0.0
    assert w.nu > 0.0


def test_negative_salinity_raises() -> None:
    with pytest.raises(ValueError, match="salinity"):
        Water(20.0, -1.0)


def test_water_is_frozen() -> None:
    w = Water()
    with pytest.raises(AttributeError):
        w.temperature = 25.0  # type: ignore[misc]

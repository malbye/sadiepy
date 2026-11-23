"""Tests for sadiepy.core.waves (theory.md §2)."""

import math
import warnings

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.optimize import brentq

from sadiepy import ValidityWarning, Water
from sadiepy.core.waves import (
    RegularWave,
    keulegan_carpenter,
    reynolds,
    ursell,
    wavenumber,
)
from sadiepy.units import GRAVITY

G = GRAVITY


# --- Spec M1 test 1: dispersion -------------------------------------------------


@pytest.mark.parametrize("y", np.logspace(-6, 4, 201))
def test_dispersion_residual_below_1e_12(y: float) -> None:
    h = 1.0
    omega = math.sqrt(y * G / h)
    k = float(wavenumber(2.0 * math.pi / omega, h))
    residual = abs(omega**2 - G * k * math.tanh(k * h)) / omega**2
    assert residual < 1e-12


def test_deep_water_limit() -> None:
    # kh ≈ 63 here, so tanh(kh) = 1 to machine precision and k = omega²/g.
    T, h = 2.0, 100.0
    omega = 2.0 * math.pi / T
    assert float(wavenumber(T, h)) == pytest.approx(omega**2 / G, rel=1e-12)


def test_shallow_water_limit() -> None:
    # For small kh, k = (omega/sqrt(gh)) * (1 + (kh)²/6 + ...); here kh ≈ 1e-4,
    # so the relative difference from the shallow-water value is ≈ 1.7e-9.
    h = 1.0
    omega = math.sqrt(1e-8 * G / h)
    k = float(wavenumber(2.0 * math.pi / omega, h))
    assert k == pytest.approx(omega / math.sqrt(G * h), rel=1e-8)


def test_wavenumber_matches_independent_root_finder() -> None:
    T, h = 1.4, 0.4
    omega = 2.0 * math.pi / T
    k_ref = brentq(lambda k: omega**2 - G * k * math.tanh(k * h), 1e-6, 100.0, xtol=1e-15)
    assert float(wavenumber(T, h)) == pytest.approx(k_ref, rel=1e-12)


def test_wavenumber_is_vectorised() -> None:
    T = np.array([1.0, 1.4, 2.0])
    h = np.array([[0.3], [0.5]])
    k = wavenumber(T, h)
    assert k.shape == (2, 3)
    assert float(k[1, 2]) == pytest.approx(float(wavenumber(2.0, 0.5)), rel=1e-15)


@pytest.mark.parametrize(("T", "h"), [(0.0, 1.0), (-1.0, 1.0), (1.0, 0.0), (math.nan, 1.0)])
def test_wavenumber_rejects_invalid_inputs(T: float, h: float) -> None:
    with pytest.raises(ValueError):
        wavenumber(T, h)


# --- Wave quantities (theory.md §2.6-2.8) --------------------------------------


@pytest.fixture
def wave() -> RegularWave:
    return RegularWave(H=0.08, T=1.4, depth=0.4)


def test_basic_relations(wave: RegularWave) -> None:
    assert wave.L == pytest.approx(2.0 * math.pi / wave.k, rel=1e-15)
    assert wave.c == pytest.approx(wave.L / wave.T, rel=1e-14)
    assert wave.kh == pytest.approx(wave.k * 0.4, rel=1e-15)


def test_group_speed_equals_numerical_derivative(wave: RegularWave) -> None:
    # c_g = d omega / dk (theory.md §2.6), checked by central differences of Eq. (2.2).
    def omega_of_k(k: float) -> float:
        return math.sqrt(G * k * math.tanh(k * wave.depth))

    dk = 1e-6 * wave.k
    cg_numeric = (omega_of_k(wave.k + dk) - omega_of_k(wave.k - dk)) / (2.0 * dk)
    assert wave.cg == pytest.approx(cg_numeric, rel=1e-8)


def test_n_limits() -> None:
    assert RegularWave(H=0.01, T=1.0, depth=50.0).n == pytest.approx(0.5, rel=1e-12)
    # Shallow water: n = 1 - (kh)²/3 + ...; here kh ≈ 6e-5, so 1 - n ≈ 1.3e-9.
    assert RegularWave(H=0.001, T=1e4, depth=0.1).n == pytest.approx(1.0, rel=1e-8)


def test_u_bar_equals_depth_average_of_u_m(wave: RegularWave) -> None:
    integral, _ = quad(lambda s: float(wave.u_m(s)), 0.0, wave.depth, epsabs=0, epsrel=1e-13)
    assert wave.u_bar == pytest.approx(integral / wave.depth, rel=1e-12)


def test_u_m_values_and_shape(wave: RegularWave) -> None:
    s = np.array([0.0, 0.2, 0.4])
    expected = wave.a * wave.omega * np.cosh(wave.k * s) / math.sinh(wave.kh)
    np.testing.assert_allclose(wave.u_m(s), expected, rtol=1e-15)
    np.testing.assert_allclose(wave.orbital_excursion(s), expected / wave.omega, rtol=1e-15)


@pytest.mark.parametrize("s", [-0.01, 0.41])
def test_u_m_rejects_heights_outside_water_column(wave: RegularWave, s: float) -> None:
    with pytest.raises(ValueError, match="heights"):
        wave.u_m(s)


def test_energy_density(wave: RegularWave) -> None:
    water = Water(20.0, 0.0)
    assert wave.energy(water) == pytest.approx(water.rho * G * 0.08**2 / 8.0, rel=1e-15)
    assert wave.energy() == wave.energy(Water())


def test_energy_flux_equals_pressure_work(wave: RegularWave) -> None:
    # theory.md Eq. (2.12): F = <∫ p u ds>, evaluated by quadrature, equals E c_g.
    water = Water()
    k, h, a, om = wave.k, wave.depth, wave.a, wave.omega
    p_amp = water.rho * G * a / math.cosh(k * h)  # times cosh(ks) cos(theta)
    u_amp = a * om / math.sinh(k * h)  # times cosh(ks) cos(theta)
    integral, _ = quad(lambda s: math.cosh(k * s) ** 2, 0.0, h, epsabs=0, epsrel=1e-13)
    flux = 0.5 * p_amp * u_amp * integral  # <cos²> = 1/2
    assert wave.energy_flux(water) == pytest.approx(flux, rel=1e-12)


# --- Dimensionless groups (theory.md §2.9) -------------------------------------


def test_ursell_forms_agree(wave: RegularWave) -> None:
    zhang_form = 8.0 * math.pi**3 * (wave.H / wave.L) / wave.kh**3
    assert wave.ursell == pytest.approx(zhang_form, rel=1e-13)
    assert ursell(0.08, 2.0, 0.4) == pytest.approx(0.08 * 4.0 / 0.064, rel=1e-15)


def test_reynolds_and_kc() -> None:
    assert reynolds(0.2, 0.01, 1e-6) == pytest.approx(2000.0, rel=1e-15)
    assert keulegan_carpenter(0.2, 1.4, 0.01) == pytest.approx(28.0, rel=1e-15)


# --- Validity checks (theory.md §2.10) -----------------------------------------


@pytest.mark.parametrize(
    ("H", "T", "h"),
    # Corners of the Zhang et al. (2023) test range (H 0.04-0.16 m, T 1.0-2.0 s,
    # h 0.3-0.5 m), all well inside the breaking limits.
    [(0.04, 1.0, 0.5), (0.16, 2.0, 0.5), (0.08, 1.4, 0.4), (0.04, 2.0, 0.3)],
)
def test_no_warning_for_moderate_waves(H: float, T: float, h: float) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        RegularWave(H=H, T=T, depth=h)


def test_depth_limited_breaking_warns() -> None:
    with pytest.warns(ValidityWarning, match="depth-limited"):
        wave = RegularWave(H=0.32, T=3.0, depth=0.4)  # H/h = 0.8
    assert wave.k > 0.0


def test_steepness_limited_breaking_warns() -> None:
    # Deep water (kh ≈ 25): limit ≈ 0.142; L ≈ 1.56 m, so H = 0.25 m gives H/L ≈ 0.16.
    with pytest.warns(ValidityWarning, match="steepness"):
        RegularWave(H=0.25, T=1.0, depth=10.0)


def test_ursell_check_is_off_by_default_and_warns_when_set() -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        RegularWave(H=0.1, T=4.0, depth=0.5)  # Ur ≈ 37, no limit set
    with pytest.warns(ValidityWarning, match="Ursell"):
        RegularWave(H=0.1, T=4.0, depth=0.5, ursell_max=10.0)


def test_invalid_height_raises() -> None:
    with pytest.raises(ValueError, match="H must be positive"):
        RegularWave(H=0.0, T=1.0, depth=1.0)


def test_wave_is_frozen(wave: RegularWave) -> None:
    with pytest.raises(AttributeError):
        wave.H = 0.1  # type: ignore[misc]

# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and
this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Documentation site (Sphinx, PyData theme, MyST-NB; Read the Docs config) with the
  theory, an API reference and example notebooks in `docs/examples`
  (`01_quickstart.ipynb`); `[examples]` and `[docs]` extras; CI docs build.
- `Water`: seawater density and viscosity from temperature and salinity (g/kg),
  after Sharqawy et al. (2010), with validity warnings.
- `RegularWave` and `core.waves`: linear wave theory (dispersion solver, wavelength,
  phase and group speed, orbital velocity, energy and energy flux, Re/KC/Ur) with
  breaking and steepness validity warnings.
- `ValidityWarning`, `ExtrapolationWarning`; `GRAVITY` = 9.80665 m/s².
- `Layer`, `Vegetation`: layered vegetation with constant or height-dependent
  frontal width.
- `drag.constant`, `drag.custom`: drag-coefficient closures (constant or a function of
  Re, KC or Ur).
- `propagate`: wave damping across a vegetation patch (analytic hyperbolic decay for
  constant C_D, ODE otherwise), with effective beta, layer fractions, energy
  bookkeeping and extrapolation warnings.
- `scaling.froude`: Froude scale factors and model/prototype conversion.
- `docs/theory.md` §1-4: water properties, linear wave theory, vegetation damping,
  Froude scaling.
- Project scaffold: `pyproject.toml` (hatchling, `src/` layout), package skeleton
  (`core`, `presets`, `scaling`, `analysis`, `design`), smoke tests.
- Tooling: pytest, ruff (lint and format), mypy (strict).
- GitHub Actions CI: lint/type-check job and test matrix for Python 3.10–3.12.
- `README.md`, `CITATION.cff`, `CHANGELOG.md`.

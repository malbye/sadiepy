# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and
this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Zhang et al. (2023) mangrove presets M1-M3 at model and prototype scale, with the
  Table 2 test matrix and the Fig. 4 damping factors; generic `rigid_cylinders` preset.
- Registry of the ten published drag closures of Zhang et al. (2023), Table 3; closures
  record source scale, vegetation and components, and may use named length scales
  (`d`, `l_e`, `l_ew`, `r_v`).
- Vegetation volumes, SVF and characteristic lengths; layer-averaged and depth-averaged
  geometries (Zhang Eqs. 19-20 and Eq. 4).
- `propagate(..., cd_evaluation="incident" | "local", scale=...)` and Re rescaling
  between Froude scales; vegetation scaling to prototype.
- Theory section 5.
- Release workflow: publishing a GitHub release uploads the same version to PyPI.
- `02_tour` notebook: a plain-English tour of every current capability, with figures.
- `01_quickstart` notebook rewritten as a step-by-step tutorial covering every object, option, result field and error message.
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

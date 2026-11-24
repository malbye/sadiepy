# SADIE: Simulation, Analysis, and Design In Ecohydraulics

<!-- Badges to add once each exists: Python versions (after a final release), DOI (Zenodo). -->
[![PyPI](https://img.shields.io/pypi/v/sadiepy.svg)](https://pypi.org/project/sadiepy/)
[![CI](https://github.com/malbye/sadiepy/actions/workflows/ci.yml/badge.svg)](https://github.com/malbye/sadiepy/actions/workflows/ci.yml)
[![Documentation](https://img.shields.io/readthedocs/sadiepy/latest)](https://sadiepy.readthedocs.io)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/malbye/sadiepy/blob/main/LICENSE)

> **Status:** early development (pre-release).

Physics-based prediction of wave damping by vegetation, with tradespace tools for
designing and analysing flume experiments.

SADIE (`sadiepy`) predicts wave attenuation through patches of coastal vegetation
(mangroves, salt marsh, seagrass, kelp, or user-defined) described as a vertical stack
of layers, using linear wave theory and the Dalrymple et al. (1984) energy-dissipation
framework. It supports flume experiment design before an experiment (tradespace
analysis of candidate designs) and data analysis after it (reflection separation,
damping-factor fitting, and drag-coefficient back-calculation).

Documentation: [https://sadiepy.readthedocs.io](https://sadiepy.readthedocs.io)

Repository: [https://github.com/malbye/sadiepy](https://github.com/malbye/sadiepy)

## Installation

SADIE requires Python 3.10&ndash;3.12. Install the latest release from PyPI:
```shell
pip install sadiepy
```

### Development Installation

To work with the source code, clone the repository and install the `sadiepy` library in "editable" mode:
```shell
git clone https://github.com/malbye/sadiepy.git
cd sadiepy
pip install -e .
```

Note: the following optional dependencies are available with bracket notation:
 * `pip install -e ".[dev]"`: for development (unit testing, linting, and type checking)
 * `pip install -e ".[examples]"`: for running the example notebooks in `docs/examples`
 * `pip install -e ".[docs]"`: for building the documentation site (includes `examples`)

## Getting Started

New to SADIE? Start with the two example notebooks in `docs/examples`, in order:

1. [**01 &ndash; Tutorial: how SADIE works**](https://github.com/malbye/sadiepy/blob/main/docs/examples/01_quickstart.ipynb):
   a slow, step-by-step introduction to every object (water, waves, vegetation layers,
   drag laws), how to run a calculation, how to read the result, and what the warnings
   and error messages mean.
2. [**02 &ndash; A tour of SADIE**](https://github.com/malbye/sadiepy/blob/main/docs/examples/02_tour.ipynb):
   a visual, plain-English overview of everything SADIE can currently do, with larger
   examples and figures.

To run them yourself, install the `examples` extra (`pip install -e ".[examples]"`) and
open them in Jupyter.

A minimal example:
```python
import sadiepy as sd

wave = sd.RegularWave(H=0.08, T=1.4, depth=0.4)  # height (m), period (s), depth (m)
veg = sd.presets.zhang2023("M3")  # 1:10 mangrove model, Zhang et al. (2023)
res = sd.propagate(wave, veg, patch_length=2.0, drag=1.0)
print(res.beta, res.H)  # damping factor (1/m), heights (m)
```

## Development Tools

Development tools are applicable when working with the source code.

### Unit Tests

Run unit tests with pytest (installed with the `dev` optional dependencies):
```shell
pytest
```

Run a subset of tests by file, or by a keyword expression with `-k`:
```shell
pytest tests/test_package.py -k "version"
```

### Documentation

Build the documentation site (Sphinx; this also runs the example notebooks in
`docs/examples`) from the project root with:
```shell
sphinx-build -b html docs docs/_build/html
```
then open `docs/_build/html/index.html`.

### Code Style

This project uses ruff for linting and formatting, applied from the project root:
```shell
ruff check .
ruff format .
```

Type checking uses mypy in strict mode (configured in `pyproject.toml`):
```shell
mypy
```

### Releasing

Releases are made from GitHub; PyPI is updated automatically. Bump `version` in
`pyproject.toml` and add a `CHANGELOG.md` entry (merged through a pull request), then
publish a GitHub release with the tag `v<version>` (e.g. `v0.1.0`). The `Release`
workflow checks that the tag matches the version, runs the tests, builds the package
and uploads it to PyPI (trusted publishing).

## License

This project is licensed under the MIT License &mdash; see
[LICENSE](https://github.com/malbye/sadiepy/blob/main/LICENSE) for details.

## Citation

If you use SADIE in your research, please cite it using the metadata in
[CITATION.cff](https://github.com/malbye/sadiepy/blob/main/CITATION.cff).

## Contact

Mallory Byers <mgbyers@ucsc.edu>

## Acknowledgements

Current Project Team
 * Mallory Byers, Independent Researcher, University of California, Santa Cruz

# SADIE: Simulation, Analysis, Design In Ecohydraulics

<!-- Badges to add once each exists: PyPI version, Python versions, CI (unit tests), documentation, DOI (Zenodo). -->
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

> **Status:** early development (pre-release).

SADIE (`sadiepy`) predicts wave attenuation through patches of coastal vegetation
(mangroves, salt marsh, seagrass, kelp, or user-defined) described as a vertical stack
of layers, using linear wave theory and the Dalrymple et al. (1984) energy-dissipation
framework. It supports flume experiment design before an experiment (tradespace
analysis of candidate designs) and data analysis after it (reflection separation,
damping-factor fitting, and drag-coefficient back-calculation).

<!-- TODO: one-line description to be finalized. -->

Repository: [https://github.com/malbye/sadiepy](https://github.com/malbye/sadiepy)

## Installation

<!-- TODO: supported Python versions to be decided. -->

`sadiepy` is not yet released on PyPI. Once released, it will install with:
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

<!-- TODO: list optional dependency groups (e.g., dev, docs, examples) once pyproject.toml is written. -->

## Development Tools

Development tools are applicable when working with the source code.

<!-- TODO: unit tests, documentation and code style sections, once tooling is decided. -->

## License

This project is licensed under the MIT License &mdash; see
[LICENSE](LICENSE) for details.

## Citation

If you use SADIE in your research, please cite it using the metadata in
[CITATION.cff](CITATION.cff).

## Contact

Mallory Byers <mgbyers@ucsc.edu>

## Acknowledgements

Current Project Team
 * Mallory Byers, Independent Researcher, University of California, Santa Cruz

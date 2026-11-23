"""Vegetation presets with sourced parameters.

``zhang2023`` builds the artificial mangrove models of Zhang et al. (2023); the
test matrix and published values are in the module
:mod:`sadiepy.presets.zhang2023` (``from sadiepy.presets.zhang2023 import TEST_MATRIX``).
"""

from sadiepy.presets.generic import rigid_cylinders
from sadiepy.presets.zhang2023 import zhang2023

__all__ = ["rigid_cylinders", "zhang2023"]

"""SADIE: Simulation, Analysis, and Design In Ecohydraulics.

Physics-based prediction of wave damping by vegetation, with tradespace tools for
designing and analysing flume experiments.
"""

from importlib.metadata import version as _version

from sadiepy.core import drag
from sadiepy.core.solver import PropagationResult, propagate
from sadiepy.core.vegetation import Layer, Vegetation
from sadiepy.core.waves import RegularWave
from sadiepy.units import GRAVITY, Water
from sadiepy.warnings import ExtrapolationWarning, ValidityWarning

__version__ = _version("sadiepy")

__all__ = [
    "GRAVITY",
    "ExtrapolationWarning",
    "Layer",
    "PropagationResult",
    "RegularWave",
    "ValidityWarning",
    "Vegetation",
    "Water",
    "__version__",
    "drag",
    "propagate",
]

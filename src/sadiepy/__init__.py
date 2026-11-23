"""SADIE: Simulation, Analysis, and Design In Ecohydraulics.

Physics-based prediction of wave damping by vegetation, with tradespace tools for
designing and analysing flume experiments.
"""

from importlib.metadata import version as _version

from sadiepy.units import Water
from sadiepy.warnings import ExtrapolationWarning, ValidityWarning

__version__ = _version("sadiepy")

__all__ = ["ExtrapolationWarning", "ValidityWarning", "Water", "__version__"]

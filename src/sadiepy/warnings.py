"""Warning categories emitted by sadiepy.

Calculations that leave the validity range of an equation or correlation never fail
silently: they return a value and emit one of these warnings. Use the standard
:mod:`warnings` filters to silence them or turn them into errors, e.g.
``warnings.simplefilter("error", sadiepy.ValidityWarning)``.
"""

__all__ = ["ExtrapolationWarning", "ValidityWarning"]


class ValidityWarning(UserWarning):
    """An input lies outside the range where an equation or theory is valid.

    Examples: water temperature or salinity outside a property correlation's range,
    or a wave that breaks or is too nonlinear for linear wave theory.
    """


class ExtrapolationWarning(ValidityWarning):
    """An empirical fit (e.g. a drag closure) is evaluated outside its fitted range.

    Subclass of :class:`ValidityWarning`, so filtering ``ValidityWarning`` also
    covers extrapolation.
    """

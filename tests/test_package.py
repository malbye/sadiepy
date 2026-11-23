import importlib

import sadiepy


def test_version_is_set() -> None:
    assert sadiepy.__version__ == "0.0.1.dev0"


def test_subpackages_import() -> None:
    for name in ("core", "presets", "scaling", "analysis", "design"):
        importlib.import_module(f"sadiepy.{name}")

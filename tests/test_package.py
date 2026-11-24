import importlib
import re
from pathlib import Path

import sadiepy


def test_version_matches_pyproject() -> None:
    text = (Path(__file__).parents[1] / "pyproject.toml").read_text(encoding="utf-8")
    match = re.search(r'^version = "([^"]+)"', text, flags=re.MULTILINE)
    assert match is not None
    assert sadiepy.__version__ == match.group(1)


def test_subpackages_import() -> None:
    for name in ("core", "presets", "scaling", "analysis", "design"):
        importlib.import_module(f"sadiepy.{name}")

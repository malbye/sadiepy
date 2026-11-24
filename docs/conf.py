"""Sphinx configuration for the SADIE documentation site."""

from importlib.metadata import version as _version

project = "SADIE"
author = "Mallory Byers"
copyright = "2026, Mallory Byers"  # noqa: A001
release = _version("sadiepy")
version = release

extensions = [
    "myst_nb",  # Markdown (MyST) pages and Jupyter notebooks
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",  # NumPy-style docstrings
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx_copybutton",
]

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]

# --- MyST / notebooks -------------------------------------------------------------
myst_enable_extensions = ["dollarmath", "amsmath", "colon_fence"]
myst_heading_anchors = 3
nb_execution_mode = "cache"  # run notebooks at build time; cache unchanged ones
nb_execution_raise_on_error = True  # a failing notebook fails the build
nb_execution_timeout = 300
nb_render_markdown_format = "myst"  # render Markdown cell outputs (e.g. tables) with MyST

# --- API reference ----------------------------------------------------------------
autosummary_generate = True
autodoc_typehints = "description"
autodoc_member_order = "bysource"
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_rtype = False

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
}

# --- HTML -------------------------------------------------------------------------
html_theme = "pydata_sphinx_theme"
html_title = "SADIE"
html_theme_options = {
    "github_url": "https://github.com/malbye/sadiepy",
    "navbar_align": "left",
    "show_toc_level": 2,
    "use_edit_page_button": False,
}
html_context = {
    "github_user": "malbye",
    "github_repo": "sadiepy",
    "github_version": "main",
    "doc_path": "docs",
}

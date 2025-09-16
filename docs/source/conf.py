# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'CASTOR Exposure Time Calculator'
copyright = '2025, CASTOR Team'
author = 'CASTOR Team'

import sys
from pathlib import Path

sys.path.insert(0, str(Path('..', '').resolve()))

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.napoleon',
              'myst_nb',
              'sphinx.ext.mathjax',  # Print mathematical expressions
              "autodoc2",
              'sphinx_rtd_theme'
              ]

# -- Options for Autodoc2 output -------------------------------------------------
# https://sphinx-autodoc2.readthedocs.io/en/stable/config.html
## TODO: redo this so that it analyzes each "submodule" on its own
autodoc2_packages = [
    "../../castor_etc"
]
autodoc2_docstring_parser_regexes = [
    # this will render all docstrings as Markdown
    (r".*", "myst"),
]

source_suffix = {
    '.rst': 'restructuredtext'
}


templates_path = ['_templates']
html_static_path = ['_static']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'

# Set the CASTOR logo for the top and for favicon
html_logo = "_static/images/castor-logo.png"
html_favicon = "_static/images/castor-favicon.ico"

html_css_files = [
    "css/custom.css"
]

# The full documentation for this is found here: https://sphinx-rtd-theme.readthedocs.io/en/stable/configuring.html
html_theme_options = {
}

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'CASTOR ETC'
copyright = '2025, CASTOR Team'
author = 'CASTOR Team'

import sys
from pathlib import Path

sys.path.insert(0, str(Path('..', '').resolve()))

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.napoleon',
              'myst_parser',
              "autodoc2"
              ]

autodoc2_packages = [
    "../../castor_etc"
]
# source_suffix = {
#     '.rst': 'restructuredtext',
#     '.ipynb': 'myst-nb',
#     '.myst': 'myst-nb',
# }
templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_book_theme'
html_static_path = ['_static']


# Automatically extract typehints when specified and place them in
# descriptions of the relevant function/method.
autodoc_typehints = "description"

# Don't show class signature with the class' name.
autodoc_class_signature = "separated"
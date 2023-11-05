# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html


import os
import sys

import autoapi
from importlib.metadata import version

# Define path to the code to be documented **relative to where conf.py (this file) is kept**
sys.path.insert(0, os.path.abspath('../src/'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "libradtranpy"
copyright = "2023, Sylvie Dagoret-Campagne"
author = "Sylvie Dagoret-Campagne"
release = version("libradtranpy")
# for example take major/minor
version = ".".join(release.split(".")[:2])

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration



extensions = [ 'sphinx.ext.autodoc',
               'sphinx.ext.duration',
               'sphinx.ext.doctest',
               'sphinx.ext.intersphinx',
               'sphinx.ext.autosummary',
               'sphinx.ext.mathjax', 
               'sphinx.ext.napoleon', 
               'sphinx.ext.viewcode',
               'sphinx_tabs.tabs',
               'numpydoc',
               'nbsphinx',
#               'sphinx.ext.graphviz',
#               'sphinx.ext.inheritance_diagram',
#               "autoapi.extension",
               ]

#extensions = [ 'sphinx.ext.autodoc',
#               'sphinx.ext.doctest',
#               'sphinx.ext.intersphinx',
#               'sphinx_design',
#              ]

nbsphinx_execute = 'never'
nbsphinx_allow_errors = True
source_suffix = ['.rst']


templates_path = []
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store','**.ipynb_checkpoints']

master_doc = "index"  # This assumes that sphinx-build is called from the root directory
html_show_sourcelink = False  # Remove 'view source code' from top of page (for html, not python)
add_module_names = False # Remove namespaces from class/method signatures

autoapi_type = "python"
autoapi_dirs = ["../src"]
autoapi_ignore = ["*/__main__.py", "*/_version.py"]
autoapi_add_toc_tree_entry = False
autoapi_member_order = "bysource"

html_theme = "sphinx_rtd_theme"

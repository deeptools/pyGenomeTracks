# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('../'))

from unittest.mock import Mock

import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)


# import mock

MOCK_MODULES = ['pygenometracks.tracksClass',  # (required for both)
                # For plotTracks in plotTracks
                'matplotlib',
                # For plotTracks in utilities
                'numpy',
                'tqdm',
                'intervaltree']

for mod_name in MOCK_MODULES:
    sys.modules[mod_name] = Mock()

autodoc_mock_imports = MOCK_MODULES
# -- Project information -----------------------------------------------------

project = 'pyGenomeTracks'
copyright = '2020, Lucille Lopez-Delisle, Leily Rabbani, Joachim Wolff, Thomas Manke, , Bjoern Gruening, Fidel Ramirez'
author = 'Lucille Lopez-Delisle, Leily Rabbani, Joachim Wolff, Thomas Manke,  Bjoern Gruening, Fidel Ramirez'


# Copied from deeptools
def get_version():
    import re
    try:
        f = open("../pygenometracks/_version.py")
    except EnvironmentError:
        return None
    for line in f.readlines():
        mo = re.match("__version__ = '([^']+)'", line)
        if mo:
            ver = mo.group(1)
            return ver
    return None


version = get_version()


# The full version, including alpha/beta/rc tags
release = version


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinxarg.ext']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'alabaster'
# Copied from deeptools
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

if not on_rtd:  # only import and set the theme if we're building docs locally
    import sphinx_rtd_theme
    html_theme = 'sphinx_rtd_theme'
    html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

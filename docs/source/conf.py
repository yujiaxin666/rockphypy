# Configuration file for the Sphinx documentation builder.
#

# -- Path setup ----------------------------------------
import os
import sys
sys.path.insert(0, os.path.abspath('../..' ))
#sys.path.insert(0, os.path.abspath('..'))
#os.path.abspath(os.path.join('..', '..'))


# from rockphypy import __version__
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'rockphypy'
copyright = '2023, jiaxin yu'
author = 'jiaxin yu'
# The short X.Y version
#version = ""
# The full version, including alpha/beta/rc tags
#release = __version__
# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    #"sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    'sphinx.ext.inheritance_diagram',
    "autoapi.extension",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinx.ext.extlinks",
    "numpydoc",
    "nbsphinx",
    #'matplotlib.sphinxext.only_directives',
    "matplotlib.sphinxext.plot_directive",
    "IPython.sphinxext.ipython_directive",
    "IPython.sphinxext.ipython_console_highlighting",
   # "sphinx.ext.inheritance_diagram",
    ]


autoapi_type = "python"
autoapi_dirs = ["../../rockphypy"]
autoapi_python_class_content = "both"
autoapi_options = [
    "members",
    "undoc-members",  # this is temporary until we add docstrings across the codebase
    "show-inheritance",
    "show-module-summary",
    "special-members",
    "imported-members",
    "inherited-members",
]

autosummary_generate = True  # Turn on sphinx.ext.autosummary
numpydoc_show_class_members = False 
autodoc_typehints = 'none'
napoleon_use_rtype = False



templates_path = ['_templates']

# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = ".rst"

exclude_patterns = ['_build',  "build", "**.ipynb_checkpoints"]

language = 'English'


# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"
# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
## adjust the white space between functions 

html_css_files = [
    'custom.css',
]

html_theme_options = {
    # 'logo': 'logo.png',
    # 'logo': 'logo.png',
    "github_user": "yujiaxin666",
    "github_repo": "RockPhysicsPy",
    "github_button": True,
    "github_banner": True,
    "travis_button": False,
    "show_powered_by": False,
    "font_family": '-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,\
        "Helvetica Neue",Arial,sans-serif,"Apple Color Emoji","Segoe UI Emoji"\
        ,"Segoe UI Symbol"',
    "font_size": "15px",
    "code_font_size": "13px",
    "head_font_family": '-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,\
        "Helvetica Neue",Arial,sans-serif,"Apple Color Emoji","Segoe UI Emoji"\
        ,"Segoe UI Symbol"',
    "page_width": "940px",
}

html_show_sphinx = False
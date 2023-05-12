# Configuration file for the Sphinx documentation builder.
#

# -- Path setup ----------------------------------------
import os
import sys
sys.path.insert(0, os.path.abspath('../..' ))
#sys.path.insert(0, os.path.abspath('..'))
#os.path.abspath(os.path.join('..', '..'))

from sphinx_gallery.sorting import FileNameSortKey
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

#autodoc_mock_imports = ['kde_diffusion']


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
    'sphinx_gallery.gen_gallery',
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

#autoapi_python_use_implicit_namespaces = False
intersphinx_mapping = {
    'numpy': ('https://numpy.org/doc/stable/', None),
    'python': ('https://docs.python.org/{.major}'.format(sys.version_info), None),
    'matplotlib': ('https://matplotlib.org/', None),
    'sphinx': ('http://www.sphinx-doc.org/en/stable', None),
    'pandas': ('https://pandas.pydata.org/pandas-docs/stable/', None),
}

autoapi_keep_files = True
autosummary_generate = False  # Turn on sphinx.ext.autosummary
numpydoc_show_class_members = False 
autodoc_typehints = 'none'
napoleon_use_rtype = False

python_version = "3.8"

templates_path = ['_templates']

# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = ['.rst']

exclude_patterns = ['_build',  "build", "**.ipynb_checkpoints"]

language = 'en'


# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"
highlight_language = 'python3'
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



# -- Sphinx Gallery Options
sphinx_gallery_conf = {
    # path to your examples scripts
    "examples_dirs": [
        "../../examples/getting_started",
        "../../examples/advanced_examples"
    ],
    # path where to save gallery generated examples
    "gallery_dirs": [
        'getting_started',
        'advanced_examples'
    ],
     # specify that examples should be ordered according to filename
    'within_subsection_order': FileNameSortKey,
    # this case sphinx_gallery and numpy in a tuple of strings.
    "filename_pattern": r"\.py",
    # Remove the "Download all examples" button from the top level gallery
    "download_all_examples": False,
    # directory where function granular galleries are stored
    "backreferences_dir": "gallery/generated/backreferences",
    # Modules for which function level galleries are created.
    "doc_module": "rockphypy"
    
}
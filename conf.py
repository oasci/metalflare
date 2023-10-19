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

# For our custom sphinx extensions
sys.path.insert(0, os.path.abspath("./website/source/exts"))


# -- Project information -----------------------------------------------------

project = "harbor"
copyright = "2023, OASCI"
author = "OASCI"
html_logo = "website/source/_static/images/harbor-logo-whale.svg"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "autoapi.extension",
    "myst_parser",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosectionlabel",
    "sphinx_multiversion",
    "sphinx_design",
    "sphinxcontrib.mermaid",
    "sphinxemoji.sphinxemoji",
    "sphinx_autodoc_typehints",
    "sphinx_copybutton",
    "sphinx_togglebutton",
    "sphinxcontrib.bibtex",
    "sphinx_inline_tabs",
]

suppress_warnings = ["autosectionlabel.*"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = [".venv", "website/html/"]

# Updating master docs
root_doc = "index"

# Add mappings
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}

# Include __init__ docstring for classes
autoclass_content = "both"

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

html_theme = "furo"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["website/source/_static"]

# Including sphinx multiversion
templates_path = [
    "website/source/_templates",
]
smv_branch_whitelist = r"main"  # Only include the main branch
html_sidebars = {
    "**": [
        "sidebar/scroll-start.html",
        "sidebar/brand.html",
        "sidebar/search.html",
        "sidebar/navigation.html",
        "sidebar/ethical-ads.html",
        "sidebar/scroll-end.html",
        "versions.html",
    ],
}

# Manually copy over files to the root. These can then be referenced outside of the
# download directive.
# html_extra_path = []


# autoapi
autoapi_type = "python"
autoapi_generate_api_docs = True
autoapi_dirs = ["02-methods/harbor"]
autoapi_add_toctree_entry = True
autoapi_python_class_content = "both"
autoapi_keep_files = False
autodoc_typehints = "description"

# bibtex
bibtex_bibfiles = ["website/source/refs.bib"]

# Header buttons
html_theme_options = {
    "sidebar_hide_name": True,
    "source_repository": "https://github.com/oasci/harbor",
    "source_branch": "main",
    "source_directory": "/",
    "top_of_page_button": ["edit", "save", "launch"],
    "path_to_docs": "/",
    "repository_url": "https://github.com/oasci/harbor",
    "repository_branch": "main",
    "launch_buttons": {
        "binderhub_url": "https://mybinder.org",
        "colab_url": "https://colab.research.google.com/",
        "deepnote_url": "https://deepnote.com/",
        "notebook_interface": "jupyterlab",
        "thebe": False,
    },
    "footer_icons": [
        {
            "name": "GitHub",
            "url": "https://github.com/oasci/harbor",
            "html": """
                <svg stroke="currentColor" fill="currentColor" stroke-width="0" viewBox="0 0 16 16">
                    <path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0 0 16 8c0-4.42-3.58-8-8-8z"></path>
                </svg>
            """,
            "class": "",
        },
    ],
}

myst_enable_extensions = [
    "dollarmath",
    "amsmath",
    "deflist",
    "html_admonition",
    "html_image",
    "colon_fence",
    "smartquotes",
    "replacements",
    "substitution",
]

togglebutton_hint = ""

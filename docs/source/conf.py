# Configuration file for the Sphinx documentation builder.
#
import os
import sys
sys.path.insert(0, os.path.abspath('../../src/'))


# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Reaction RDF Converter'
copyright = '2025, CWRU-SDLE'
author = 'Quynh D. Tran, Holly Schreiber, Owen Schessler, Laura S. Bruckman, Roger H. French'
release = '0.1.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',          # ESSENTIAL: This enables automodule directives
    'sphinx.ext.viewcode',         # Add source code links
    'sphinx.ext.napoleon',         # Support for Google/NumPy style docstrings
    'sphinx.ext.autosummary',      # Generate summary tables
    'sphinx.ext.intersphinx',      # Cross-reference external docs
    'sphinx.ext.githubpages',      # GitHub Pages compatibility
]

# Autodoc configuration
autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'special-members': '__init__',
    'undoc-members': True,
    'exclude-members': '__weakref__'
}

# Mock imports for dependencies that might not be available during doc build
autodoc_mock_imports = [
    'ord_schema',
    'rdkit',
    'protobuf',
]

# Napoleon settings for your docstring style
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = True

# Intersphinx mapping for cross-references
intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'pandas': ('https://pandas.pydata.org/docs/', None),
    'ord_schema': ('https://docs.open-reaction-database.org/en/latest/', None),
    'rdkit': ('https://www.rdkit.org/docs/index.html', None),
    'protobufs': ('https://googleapis.dev/python/protobuf/latest/', None),
}

# Autosummary settings
autosummary_generate = True
autosummary_imported_members = True

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# Theme options
html_theme_options = {
    'canonical_url': '',
    'analytics_id': '',
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False
}

# Custom sidebar
html_sidebars = {
    '**': [
        'about.html',
        'navigation.html',
        'relations.html',
        'searchbox.html',
        'donate.html',
    ]
}

# Additional HTML options
html_show_sourcelink = True
html_show_sphinx = True
html_show_copyright = True

# Custom CSS/JS files (if you create them)
html_css_files = [
    # 'custom.css',  # Uncomment if you create custom CSS
]

html_js_files = [
    # 'custom.js',   # Uncomment if you create custom JS
]

# Logo and favicon (add these files to _static/ if desired)
# html_logo = '_static/logo.png'
# html_favicon = '_static/favicon.ico'

# Project description for metadata
html_title = f'{project} v{release}'
html_short_title = project

# Social/project links
html_context = {
    'display_github': True,
    'github_user': 'quynhdtran17',
    'github_repo': 'rxn_rdf_converter',
    'github_version': 'main',
    'conf_py_path': '/docs/source/',
}

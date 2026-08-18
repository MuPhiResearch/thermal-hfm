# Configuration file for the Sphinx documentation builder.
#
# Thermal EDFM/pEDFM Module - MRST 2026a
# ----------------------------------------

import os
import sys

# -- Project information -----------------------------------------------------
project = 'Thermal HFM'
copyright = '2025-2026, Reza Najafi Silab, David Egya, Florian Doster'
author = 'Reza Najafi Silab, David Egya, Florian Doster'
release = '1.0'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.mathjax',
]

# Code blocks default to MATLAB syntax highlighting
highlight_language = 'matlab'

# Template and static paths
templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_css_files = ['custom.css']

html_logo = '_static/logo.jpg'

html_theme_options = {
    'logo_only': True,
    'navigation_depth': 4,
    'collapse_navigation': False,
    'titles_only': False,
}

# -- MathJax -----------------------------------------------------------------
mathjax3_config = {
    'tex': {
        'macros': {
            'lR': r'\lambda_{\mathrm{R}}',
            'lF': r'\lambda_{\mathrm{F}}',
            'CI': r'\mathrm{CI}',
        }
    }
}

# Number figures/tables so they can be cross-referenced
numfig = True

# -- Options for LaTeX / PDF output ------------------------------------------
# Preamble providing the TikZ setup and colours used by chart.tex, which is
# pulled into the LaTeX/PDF build via "\input{chart.tex}" from models.rst.
latex_elements = {
    'preamble': r'''
\usepackage{tikz}
\usepackage{adjustbox}
\usepackage{float}
\usetikzlibrary{arrows.meta,shadows,positioning,calc}
% Macros used by the figure caption
\providecommand{\lR}{\lambda_{\mathrm{R}}}
\providecommand{\lF}{\lambda_{\mathrm{F}}}
% Class-hierarchy figure colours (chart.tex)
\definecolor{hier-parent}{RGB}{207,216,220}
\definecolor{hier-hfm}{RGB}{197,225,165}
\definecolor{hier-state}{RGB}{206,194,233}
\definecolor{hier-preproc}{RGB}{255,224,178}
\definecolor{hier-ext}{RGB}{236,239,241}
\definecolor{primary}{RGB}{21,101,192}
\definecolor{accent-teal}{RGB}{0,137,123}
\definecolor{accent-orange}{RGB}{239,108,0}
''',
}

# Copy chart.tex into the LaTeX build dir so "\input{chart.tex}" resolves.
latex_additional_files = ['chart.tex']

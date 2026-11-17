# Sphinx configuration for the PyLightcurve 4.0 documentation.
import shutil
import sys
from pathlib import Path

DOCS_DIR = Path(__file__).resolve().parent
REPO_DIR = DOCS_DIR.parent
sys.path.insert(0, str(REPO_DIR))

project = 'PyLightcurve'
author = 'Angelos Tsiaras and collaborators'
copyright = '2016-present, Angelos Tsiaras and collaborators'

release = (REPO_DIR / 'pylightcurve' / '__version__.txt').read_text().strip()
version = '.'.join(release.split('.')[:2])

LATEST_DOCS_URL = 'https://github.com/ucl-exoplanets/pylightcurve'

(DOCS_DIR / '_static').mkdir(exist_ok=True)
shutil.copy2(REPO_DIR / 'logo.png', DOCS_DIR / '_static' / 'logo.png')

extensions = [
    'myst_parser',                # Markdown pages
    'sphinx.ext.autodoc',         # API reference from docstrings
    'sphinx.ext.autosummary',     # summary tables in the API reference
    'sphinx.ext.napoleon',        # NumPy-style docstrings
    'sphinx.ext.mathjax',         # LaTeX maths
    'sphinx.ext.intersphinx',     # links to numpy / astropy / emcee docs
    'sphinx.ext.viewcode',        # "[source]" links
    'sphinx_copybutton',          # copy button on code blocks
    'sphinx_design',              # cards on the home page
]

source_suffix = {'.md': 'markdown', '.rst': 'restructuredtext'}
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

myst_enable_extensions = ['dollarmath', 'amsmath', 'colon_fence', 'deflist']
myst_heading_anchors = 3

autosummary_generate = True
autodoc_member_order = 'bysource'
autodoc_default_options = {
    'members': True,
    'undoc-members': False,
    'show-inheritance': False,
}
autodoc_typehints = 'none'

napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_rtype = False

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
    'astropy': ('https://docs.astropy.org/en/stable/', None),
    'emcee': ('https://emcee.readthedocs.io/en/stable/', None),
    'matplotlib': ('https://matplotlib.org/stable/', None),
}

html_theme = 'pydata_sphinx_theme'
html_title = f'PyLightcurve {release}'
html_logo = '_static/logo.png'
html_static_path = ['_static']
html_theme_options = {
    'announcement': (
        f'This is the documentation of PyLightcurve {version}, the long-term support '
        f'version. For the latest version, see the <a href="https://github.com/ucl-exoplanets/pylightcurve" target="_blank">main Pylightcurve repository</a>.'
    ),
    'github_url': 'https://github.com/atsiaras/pylightcurve-4',
    'icon_links': [
        {
            'name': 'PyPI',
            'url': 'https://pypi.org/project/pylightcurve/',
            'icon': 'fa-brands fa-python',
        },
    ],
    'navbar_align': 'left',
    'show_toc_level': 2,
    'navigation_with_keys': False,
}
html_context = {
    'github_user': 'atsiaras',
    'github_repo': 'pylightcurve-4',
    'github_version': 'main',
    'doc_path': 'docs',
}

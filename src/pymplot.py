"""
pymplot — convenience shim for the src/api.py function API.

Usage:
    import pymplot
    pymplot.showmatrix(data)

    # or with direct imports:
    from pymplot import showmatrix, showwiggle
    showmatrix(data)
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'src'))

from api import (
    showmatrix,
    showcontour,
    showslice,
    showslicon,
    showvolume,
    showvolcon,
    showwiggle,
    showcolorbar,
    showgraph,
)

__all__ = [
    'showmatrix', 'showcontour', 'showslice', 'showslicon',
    'showvolume', 'showvolcon', 'showwiggle', 'showcolorbar', 'showgraph',
]

"""
pymplot — function API for plotting regular gridded data.

Usage
-----
    import pymplot as pm
    import numpy as np

    data2d = np.random.randn(100, 200)
    pm.showmatrix(data2d)
    pm.showmatrix(data2d, colormap='gray', legend=True, outfile='out.pdf')

    data3d = np.random.randn(50, 80, 60)
    pm.showslice(data3d)
    pm.showvolume(data3d, outfile='vol.pdf')

All keyword arguments map directly to the corresponding CLI parameters.
"""

import os
import sys

# Ensure the src directory is on sys.path so the module files are importable
_src_dir = os.path.dirname(os.path.abspath(__file__))
if _src_dir not in sys.path:
    sys.path.insert(0, _src_dir)

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
    'showmatrix',
    'showcontour',
    'showslice',
    'showslicon',
    'showvolume',
    'showvolcon',
    'showwiggle',
    'showcolorbar',
    'showgraph',
]

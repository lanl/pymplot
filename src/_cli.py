"""
Console-script entry points for the pymplot CLI tools.
Each function is wired to an x_show* command by pyproject.toml.
"""

import sys
import os

# Ensure src/ is on the path (needed when invoked as an installed entry point)
_src_dir = os.path.dirname(os.path.abspath(__file__))
if _src_dir not in sys.path:
    sys.path.insert(0, _src_dir)

import argparse
from argparse import RawTextHelpFormatter


def _run(program, module_name):
    import importlib
    mod = importlib.import_module(module_name)
    from module_getarg import getarg
    print()
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser = getarg(parser, program)
    args = parser.parse_args()
    mod.run(args)


def showmatrix():   _run('matrix',   'showmatrix')
def showcontour():  _run('contour',  'showcontour')
def showwiggle():   _run('wiggle',   'showwiggle')
def showgraph():    _run('graph',    'showgraph')
def showslice():    _run('slice',    'showslice')
def showslicon():   _run('slicon',   'showslicon')
def showvolume():   _run('volume',   'showvolume')
def showvolcon():   _run('volcon',   'showvolcon')
def showcolorbar(): _run('colorbar', 'showcolorbar')

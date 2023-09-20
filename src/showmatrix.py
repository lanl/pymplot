#!/usr/bin/python3

from pylab import *
import matplotlib as mplt
import numpy as np
import matplotlib.pyplot as plt
import os
import math
import argparse
from module_getarg import getarg
from argparse import RawTextHelpFormatter
from module_io import *

# this is to ignore warnings
import warnings

warnings.filterwarnings("ignore", module="matplotlib")

# tag
program = 'matrix'
print()

## read arguments
parser = argparse.ArgumentParser(description='''
                                purpose:
                                    Plot a 2D array as an image.
                                ''',
                                 formatter_class=RawTextHelpFormatter)
parser = getarg(parser, program)
args = parser.parse_args()

# input
data, n1, n2, dmin, dmax = read_array(args, which='fore', dim=2)
backdata, _, _, backdmin, backdmax = read_array(args, which='back', dim=2)
mask, _, _, _, _ = read_array(args, which='mask', dim=2)
if mask is not None:
    data = data * mask
if mask is not None and backdata is not None:
    backdata = backdata * mask

# dimension information
d1 = float(args.d1)
d2 = float(args.d2)

## limit of axis
from module_range import *
sp1beg, sp1end, x1beg, x1end, n1beg, n1end = set_range(args.o1, n1, d1, args.x1beg, args.x1end)
sp2beg, sp2end, x2beg, x2end, n2beg, n2end = set_range(args.o2, n2, d2, args.x2beg, args.x2end)
data = data[n1beg:n1end, n2beg:n2end]
if args.background is not None:
    backdata = backdata[n1beg:n1end, n2beg:n2end]

## figure size
from module_size import *
figheight, figwidth = set_size(args, n1beg, n1end, n2beg, n2end)

## begin plot
if not args.imageonly:
    frameon = True
else:
    frameon = False

fig = plt.figure(figsize=(figwidth, figheight), frameon=frameon)
ax = fig.add_axes([0, 0, 1, 1])

if not args.imageonly:
    ax.set_axis_on()
else:
    ax.set_axis_off()

# show the image
im = ax.imshow(data, zorder=2)
im.set_extent([0, figwidth, figheight, 0])
if args.background is not None:
    imb = ax.imshow(backdata, zorder=1)
    imb.set_extent([0, figwidth, figheight, 0])

## set frame
from module_frame import *
set_frame(args)

## set clip
from module_clip import *
cmin, cmax = set_clip(args, data, 'fore', dmin, dmax)
if args.norm == 'log':
    if cmin > np.floor(cmax) or cmax < np.ceil(cmin):
        print(' Error: Values in dataset have same order of magnitude. ')
        exit()
im.set_clim(cmin, cmax)

if args.background is not None:
    backcmin, backcmax = set_clip(args, backdata, 'back', backdmin, backdmax)
    if args.norm == 'log':
        if backcmin > np.floor(backcmax) or backcmax < np.ceil(backcmin):
            print(' Error: Values in dataset have same order of magnitude. ')
            exit()
    imb.set_clim(backcmin, backcmax)

## set colormap
from module_colormap import set_colormap, set_colormap_alpha
colormap = set_colormap(args)
colormap = set_colormap_alpha(args, colormap, cmin, cmax)
im.set_cmap(colormap)
if args.background is not None:
    colormap = set_colormap(args, 'background')
    colormap = set_colormap_alpha(args, colormap, backcmin, backcmax, 'background')
    imb.set_cmap(colormap)

## set interpolation
im.set_interpolation(args.interp)
if args.background is not None:
    imb.set_interpolation(args.backinterp)

## set font
from module_font import *
font, fontbold = set_font(args)

## set tick
from module_tick import *
set_tick(args, font, x1beg, x1end, n1beg, n1end, d1, figheight, x2beg, x2end, n2beg, n2end, d2, figwidth)

## set grid line
from module_gridline import *
set_gridline(args)

## set title
from module_title import *
set_title(args, fontbold)

## set annotation
from module_annotation import *
set_annotation(args, font, x1beg, n1end - n1beg, d1, figheight, x2beg, n2end - n2beg, d2, figwidth)

## set colorbar
from module_colorbar import *
if args.legend and not args.backlegend:
    set_colorbar(args, im, font, cmin, cmax, figheight, figwidth, fig)
if args.backlegend:
    set_colorbar(args, imb, font, backcmin, backcmax, figheight, figwidth, fig)

## set shading plot
if args.shading is not None:

    from matplotlib.colors import LightSource
    angle = args.shading_angle.split(',')
    az = np.float32(angle[0])
    alt = np.float32(angle[1])
    ls = LightSource(azdeg=az, altdeg=alt)
    rgb = ls.shade(data, cmap=colormap, vert_exag=args.shading_scale, blend_mode=args.shading)
    im = ax.imshow(rgb, cmap=colormap, extent=[0, figwidth, figheight, 0], zorder=2)

## axis invert
if args.reverse1:
    ax.invert_yaxis()
if args.reverse2:
    ax.invert_xaxis()

## output
output(args)

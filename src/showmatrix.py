#!/usr/bin/python3

## dependencies
from pylab import *
import matplotlib as mplt
import numpy as np
import matplotlib.pyplot as plt
import os
import math
import argparse
from module_getarg import *
from argparse import RawTextHelpFormatter

# this is to ignore warnings
import warnings
warnings.filterwarnings("ignore", module="matplotlib")

## read arguments
# assign description to the help doc
parser = argparse.ArgumentParser(description='''
    Functionality: Read a 2D array from raw binary file and plot as an image. 
    Author:        K.G. @ 2016.05, 2016.06, 2016.08, 2016.10, 2017.03, 
                          2019.10, 2020.05''',
                                 formatter_class=RawTextHelpFormatter)

# arguments -- general
parser = getarg_general(parser)

# arguments -- color
parser = getarg_color(parser)

# arguments -- axis
parser = getarg_axis(parser, 2)

# arguments -- tick
parser = getarg_tick(parser, 2)

# arguments -- colorbar
parser = getarg_colorbar(parser)

# arguments -- title
parser = getarg_title(parser)

# arguments -- annotation
parser = getarg_annotation(parser)

# array for all arguments passed to script
args = parser.parse_args()

## input data
infile = args.infile[0]

if not os.path.exists(infile):
    print()
    print(' Error: File', infile, 'does not exists; Exit. ')
    print()
    exit()

if len(args.background) != 0 and (not os.path.exists(args.background)):
    print()
    print(' Error: File', args.background, 'does not exists. Exit. ')
    print()
    exit()

fsize = os.path.getsize(infile)
datatype = args.datatype
if datatype == 'double':
    fsize = fsize / 8
if datatype == 'float':
    fsize = fsize / 4
if datatype == 'int':
    fsize = fsize / 2

n1 = args.n1
if args.n2 == 0:
    n2 = int(fsize * 1.0 / n1)
else:
    n2 = args.n2

# data type
from module_datatype import *
dt = set_datatype(args)

data = np.empty([n1, n2])
data = fromfile(infile, dtype=dt, count=n1 * n2)
if len(args.background) != 0:
    backdata = np.empty([n1, n2])
    backdata = fromfile(args.background, dtype=dt, count=n1 * n2)

# reshape array
if args.transpose == 0:
    data = data.reshape((n2, n1))
    data = data.transpose()
else:
    data = data.reshape(n1, n2)

if (args.flip1 == 1):
    data = np.flipud(data)
if (args.flip2 == 1):
    data = np.fliplr(data)

if len(args.background) != 0:
    if args.transpose == 0:
        backdata = backdata.reshape((n2, n1))
        backdata = backdata.transpose()
    else:
        backdata = backdata.reshape(n1, n2)

    if (args.flip1 == 1):
        backdata = np.flipud(backdata)
    if (args.flip2 == 1):
        backdata = np.fliplr(backdata)

if len(args.mask) != 0:
    mk = np.empty([n1, n2])
    mk = fromfile(args.mask, dtype=np.float32, count=n1 * n2)

    if args.transpose == 0:
        mk = mk.reshape((n2, n1))
        mk = mk.transpose()
    else:
        mk = mk.reshape(n1, n2)

    if (args.flip1 == 1):
        mk = np.flipud(mk)
    if (args.flip2 == 1):
        mk = np.fliplr(mk)

    # foreground image masking
    data = np.ma.masked_array(data, mask=(abs(mk - 1) > 1.0e-5))

    # background image masking
    if len(args.background) != 0:
        backdata = np.ma.masked_array(backdata, mask=(abs(mk - 1) > 1.0e-5))

# scale
data = data * eval(args.colorscale)

# data min and max values
if isnan(sum(data)) == True:
    udata = data[~isnan(data)]
    if udata.shape == (0, ):
        print(' Error: Input dataset is all NaN; Exit. ')
        exit()
    else:
        dmin = udata.min()
        dmax = udata.max()
else:
    dmin = data.min()
    dmax = data.max()

# print value range of input matrix
print()
print('input <<    ', infile)
print('shape       ', data.shape)
print('value range ', dmin, ' -- ', dmax)

# dimension information
d1 = float(args.d1)
d2 = float(args.d2)

## limit of axis
from module_range import *
sp1beg, sp1end, x1beg, x1end, n1beg, n1end = set_range(args.f1, n1, d1, args.x1beg, args.x1end)
sp2beg, sp2end, x2beg, x2end, n2beg, n2end = set_range(args.f2, n2, d2, args.x2beg, args.x2end)
data = data[n1beg:n1end, n2beg:n2end]
if args.norm == 'log': data = np.log10(data)
if len(args.background) != 0:
    backdata = backdata[n1beg:n1end, n2beg:n2end]

## figure size
from module_size import *
figheight, figwidth = set_size(args, n1beg, n1end, n2beg, n2end)

## begin plot
if args.imageonly == 0:
    frameon = True
else:
    frameon = False

fig = plt.figure(figsize=(figwidth, figheight), frameon=frameon)
ax = fig.add_axes([0, 0, 1, 1])

if args.imageonly == 0:
    ax.set_axis_on()
else:
    ax.set_axis_off()

im = ax.imshow(data, zorder=2)
im.set_extent([0, figwidth, figheight, 0])
if len(args.background) != 0:
    imb = ax.imshow(backdata, zorder=1)
    imb.set_extent([0, figwidth, figheight, 0])

## set frame
from module_frame import *
set_frame(args)

## set clip
from module_clip import *
cmin, cmax = set_clip(args, data, 'fore')
if args.norm == 'log':
    if cmin > np.floor(cmax) or cmax < np.ceil(cmin):
        print('error: values in dataset have same order of magnitude')
        exit()
im.set_clim(cmin, cmax)

if len(args.background) != 0:
    backcmin, backcmax = set_clip(args, backdata, 'back')
    if args.norm == 'log':
        if backcmin > np.floor(backcmax) or backcmax < np.ceil(backcmin):
            print('error: values in dataset have same order of magnitude')
            exit()
    imb.set_clim(backcmin, backcmax)

## set colormap
from module_colormap import set_colormap, set_colormap_alpha
colormap = set_colormap(args)
colormap = set_colormap_alpha(args, colormap, cmin, cmax)
im.set_cmap(colormap)
if len(args.background) != 0:
    colormap = set_colormap(args, 'background')
    colormap = set_colormap_alpha(args, colormap, backcmin, backcmax, 'background')
    imb.set_cmap(colormap)

## set interpolation
im.set_interpolation(args.interp)
if len(args.background) != 0:
    imb.set_interpolation(args.backinterp)

## set font
from module_font import *
font, fontbold = set_font(args)

## set tick
from module_tick import *
set_tick(args,font, \
    x1beg,x1end,n1beg,n1end,d1,figheight, \
    x2beg,x2end,n2beg,n2end,d2,figwidth)

## set grid line
from module_gridline import *
set_gridline(args)

## set title
from module_title import *
set_title(args, fontbold)

## set annotation
from module_annotation import *
set_annotation(args,font, \
    x1beg,n1end-n1beg,d1,figheight, \
    x2beg,n2end-n2beg,d2,figwidth)

## set colorbar
from module_colorbar import *
if args.legend == 1 and args.backlegend == 0:
    set_colorbar(args, im, font, cmin, cmax, figheight, figwidth, fig)
if args.backlegend == 1:
    set_colorbar(args, imb, font, backcmin, backcmax, figheight, figwidth, fig)

## axis invert
if args.reverse1 == 1:
    ax.invert_yaxis()
if args.reverse2 == 1:
    ax.invert_xaxis()

## output
from module_output import *
output(args)

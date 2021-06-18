#!/usr/bin/python3

## dependencies
from pylab import *
import matplotlib as mplt
import numpy as np
import matplotlib.pyplot as plt
import os
import math
import argparse
from module_getarg import getarg
from argparse import RawTextHelpFormatter
from module_contour import *
from module_io import *

# ignore warnings
import warnings
warnings.filterwarnings("ignore", module="matplotlib")

# tag
program = 'slicon'
print()

## read arguments
parser = argparse.ArgumentParser(description='''
                                purpose:
                                    Plot a 3D array as three orthogonal slices of contours.
                                ''',
                                 formatter_class=RawTextHelpFormatter)
parser = getarg(parser, program)
args = parser.parse_args()

# input
data, n1, n2, n3, dmin, dmax = read_array(args, which='fore', dim=3)
backdata, _, _, _, backdmin, backdmax = read_array(args, which='back', dim=3)
mask, _, _, _, _, _ = read_array(args, which='mask', dim=3)
if mask is not None:
    data = data * mask
if mask is not None and backdata is not None:
    backdata = backdata * mask

d1 = float(args.d1)
d2 = float(args.d2)
d3 = float(args.d3)

## set clip
from module_clip import *
cmin, cmax = set_clip(args, data, 'fore', dmin, dmax)
if args.norm == 'log':
    if cmin > np.floor(cmax) or cmax < np.ceil(cmin):
        print('error: values in dataset have same order of magnitude')
        exit()

if args.background is not None:
    backcmin, backcmax = set_clip(args, backdata, 'back', backdmin, backdmax)
else:
    backcmin = cmin
    backcmax = cmax

## limit of axis
from module_range import *
sp1beg, sp1end, x1beg, x1end, n1beg, n1end = set_range(args.o1, n1, d1, args.x1beg, args.x1end)
sp2beg, sp2end, x2beg, x2end, n2beg, n2end = set_range(args.o2, n2, d2, args.x2beg, args.x2end)
sp3beg, sp3end, x3beg, x3end, n3beg, n3end = set_range(args.o3, n3, d3, args.x3beg, args.x3end)

## set slice
from module_utility import *
# axis 1
if args.slice1 is None:
    sl1 = rounddecimalbase((x1end + x1beg) / 2.0, d1)
else:
    sl1 = eval(args.slice1)

# axis 2
if args.slice2 is None:
    sl2 = rounddecimalbase((x2end + x2beg) / 2.0, d2)
else:
    sl2 = eval(args.slice2)

# axis 3
if args.slice3 is None:
    sl3 = rounddecimalbase((x3end + x3beg) / 2.0, d3)
else:
    sl3 = eval(args.slice3)

# slice index
slice1 = int(round((sl1 - sp1beg) / d1))
slice2 = int(round((sl2 - sp2beg) / d2))
slice3 = int(round((sl3 - sp3beg) / d3))

if slice1 <= 0 or slice1 >= n1end:
    print('error: slice 1 selection error')
    exit()
if slice2 <= 0 or slice2 >= n2end:
    print('error: slice 2 selection error')
    exit()
if slice3 <= 0 or slice3 >= n3end:
    print('error: slice 3 selection error')
    exit()

sl1 = sl1 - x1beg + 0.5 * d1
sl2 = sl2 - x2beg + 0.5 * d2
sl3 = sl3 - x3beg + 0.5 * d3

# slice yz
data12 = data[n1beg:n1end, n2beg:n2end, slice3]
if args.norm == 'log':
    data12 = np.log10(data12)

# slice yx
data13 = data[n1beg:n1end, slice2, n3beg:n3end]
if args.norm == 'log':
    data13 = np.log10(data13)

# slice xz
data23 = data[slice1, n2beg:n2end, n3beg:n3end]
data23 = flipud(data23)
if args.norm == 'log':
    data23 = np.log10(data23)

# background data
if args.background is not None:

    # slice yz
    backdata12 = backdata[n1beg:n1end, n2beg:n2end, slice3]

    # slice yx
    backdata13 = backdata[n1beg:n1end, slice2, n3beg:n3end]

    # slice xz
    backdata23 = backdata[slice1, n2beg:n2end, n3beg:n3end]
    backdata23 = flipud(backdata23)

else:

    backdata12 = []
    backdata13 = []
    backdata23 = []

## set figure size
# inch per point
ipp = 0.0138889

# default longest axis of three figures is 5 inch
figbase = 5.0
golden_ratio = 1.0 / 1.61803398875
nmax = max(n1end - n1beg, n2end - n2beg, n3end - n3beg)

# if figure width/height or vice versa larger than 6 then use golden ratio
limit = 6.0

if args.size1 is None:
    ratio = float(n1end - n1beg) / nmax
    if ratio < 1.0 / limit:
        ratio = golden_ratio
    size1 = figbase * ratio
else:
    size1 = float(args.size1)

if args.size2 is None:
    ratio = float(n2end - n2beg) / nmax
    if ratio < 1.0 / limit:
        ratio = golden_ratio
    size2 = figbase * ratio
else:
    size2 = float(args.size2)

if args.size3 is None:
    ratio = float(n3end - n3beg) / nmax
    if ratio < 1.0 / limit:
        ratio = golden_ratio
    size3 = figbase * ratio
else:
    size3 = float(args.size3)

## set font
from module_font import *
font, fontbold = set_font(args)

# if automatically determine contour levels
mcontour = int(args.mcontour) + 1
if mcontour <= 0:
    print('sublevel contour specification error')
    exit()

# linear data norm
if args.norm == 'linear':

    if args.contours is None:

        if args.contourlevel is None:
            ctrd = nice((cmax - cmin) / 10.0)
        else:
            ctrd = float(args.contourlevel)

        if args.contourbeg is None:
            ctrbeg = nice(cmin)
            base = 0.5
            nb = 0
            while nb <= 10 and ctrbeg > cmin + ctrd:
                base = base / 10.0
                ctrbeg = nice(cmin, base)
                nb = nb + 1
        else:
            ctrbeg = float(args.contourbeg)
        if args.contourend is None:
            ctrend = cmax
        else:
            ctrend = float(args.contourend)

        # contour levels
        levels = np.arange(ctrbeg, ctrend + 0.5 * abs(ctrd) / float(mcontour), ctrd / float(mcontour))
        levels = np.append(levels, levels[-1])

    # if explicitly specify contour levels
    else:

        # major levels
        levels = args.contours[0].split(',')
        for i in range(0, size(levels)):
            levels[i] = float(levels[i])

        # calculate minor levels
        tls = levels[0]
        for i in range(0, size(levels) - 1):
            tls = np.append(tls, np.linspace(levels[int(i)], levels[int(i) + 1], mcontour + 1))
        levels = unique(tls)
        levels = np.append(levels, levels[-1])

# log data norm
if args.norm == 'log':

    if args.contours is None:

        if args.contourbeg is None:
            ctrbeg = np.floor(cmin)
        else:
            ctrbeg = float(args.contourbeg)
        if args.contourend is None:
            ctrend = np.ceil(cmax) + 1
        else:
            ctrend = float(args.contourend)

        if args.contourlevel is None:
            ctrd = max(1, int((ctrbeg - ctrend) / 5.0))
        else:
            ctrd = int(args.contourlevel)

        # contour levels
        levels = np.arange(ctrbeg, ctrend + ctrd, ctrd)
        levels = np.append(levels, levels[-1] + 1)
        nl = len(levels)
        mlevels = []
        for i in range(0, nl - 1):
            mlevels = np.append(mlevels,
                                np.log10(np.linspace(10**levels[i], 10**levels[i + 1], args.mcontour + 2)))
        levels = np.unique(mlevels)
        levels = np.append(levels, levels[-1])

    # if explicitly specify contour levels
    else:

        # major levels
        levels = args.contours[0].split(',')
        for i in range(0, size(levels)):
            levels[i] = float(levels[i])

        # calculate minor levels
        tls = levels[0]
        for i in range(0, size(levels) - 1):
            tls = np.append(tls,
                            np.log10(np.linspace(10**levels[int(i)], 10**levels[int(i) + 1], mcontour + 1)))
        levels = unique(tls)
        levels = np.append(levels, levels[-1])

# contour font size
if args.clabelsize is None:
    clabelsize = min(float(args.label1size), float(args.label2size)) - 1
else:
    clabelsize = float(args.clabelsize)

# contour widths
if args.mcontourwidth is None:
    mw = 0.25 * float(args.contourwidth)
else:
    mw = float(args.mcontourwidth)

lw = np.array([mw for i in range(0, size(levels))])
lw[0:-1:mcontour] = float(args.contourwidth)

ls = np.array([args.mcontourstyle for i in range(0, size(levels))])
ls[0:-1:mcontour] = args.contourstyle

lc0 = args.contourcolor[0].split(',')
lc = ['k' for i in range(0, size(levels))]
lc[0:size(lc0)] = lc0

## plot
figheight = size1 + size2 + float(args.slicegap)
figwidth = size2 + size3 + float(args.slicegap)
fig = plt.figure(figsize=(figwidth, figheight))

# image 21
ax21locx = 0.0
ax21locy = 0.0
im21height = size1 / figheight
im21width = size3 / figwidth

ax21 = fig.add_axes([ax21locx, ax21locy, im21width, im21height])
im21 = add_contour(args, size3, size1, n1beg, n1end, n3beg, n3end, ax21, data13, backdata13, levels, lc, lw,
                   ls, clabelsize, mcontour, font, cmin, cmax, backcmin, backcmax)

#image 11
ax11locx = ax21locx
ax11locy = ax21locy + im21height + float(args.slicegap) / figheight
im11height = size2 / figheight
im11width = size3 / figwidth

ax11 = fig.add_axes([ax11locx, ax11locy, im11width, im11height])
im11 = add_contour(args, size3, size2, n2beg, n2end, n3beg, n3end, ax11, np.flipud(data23),
                   np.flipud(backdata23), levels, lc, lw, ls, clabelsize, mcontour, font, cmin, cmax,
                   backcmin, backcmax)
ax11.invert_yaxis()

# image 22
ax22locx = ax21locx + im21width + float(args.slicegap) / figwidth
ax22locy = ax21locy
im22height = size1 / figheight
im22width = size2 / figwidth

ax22 = fig.add_axes([ax22locx, ax22locy, im22width, im22height])
im22 = add_contour(args, size2, size1, n1beg, n1end, n2beg, n2end, ax22, data12, backdata12, levels, lc, lw,
                   ls, clabelsize, mcontour, font, cmin, cmax, backcmin, backcmax)

# 3d image at image 12
if args.topright is not None:
    volslice = plt.imread(args.topright)
    ax12 = fig.add_axes([ax22locx, ax11locy, im22width, im11height])
    ax12.imshow(volslice, aspect=1, interpolation='kaiser')
    ax12.set_axis_off()

## set tick
from module_tick import *

# axes labels
label1size = float(args.label1size)
label2size = float(args.label2size)
label3size = float(args.label3size)

ax11.set_ylabel(args.label2, fontsize=label2size, labelpad=float(args.label2pad)*100, fontproperties=font)
ax21.set_ylabel(args.label1, fontsize=label1size, labelpad=float(args.label1pad)*100, fontproperties=font)
ax21.set_xlabel(args.label3, fontsize=label3size, labelpad=float(args.label3pad)*100, fontproperties=font)
ax22.set_xlabel(args.label2, fontsize=label2size, labelpad=float(args.label2pad)*100, fontproperties=font)
ax21.yaxis.set_label_position('left')

l = ax11.yaxis.get_label()
l.set_fontproperties(font)
l.set_fontsize(label2size)
l = ax21.yaxis.get_label()
l.set_fontproperties(font)
l.set_fontsize(label1size)
l = ax21.xaxis.get_label()
l.set_fontproperties(font)
l.set_fontsize(label3size)
l = ax22.xaxis.get_label()
l.set_fontproperties(font)
l.set_fontsize(label2size)

# ticks on/off
ax11.get_yaxis().set_tick_params(which='both', direction='out')
ax21.get_xaxis().set_tick_params(which='both', direction='out')
ax21.get_yaxis().set_tick_params(which='both', direction='out')
ax22.get_xaxis().set_tick_params(which='both', direction='out')

# tick location, label and minor tick location
tick_1_location, tick_1_label, tick_1_minor = define_tick(args.ticks1, args.tick1beg, args.tick1end,
                                                          args.tick1d, args.mtick1, x1beg, x1end,
                                                          n1end - n1beg + 1, d1, size1, args.tick1format)

tick_2_location, tick_2_label, tick_2_minor = define_tick(args.ticks2, args.tick2beg, args.tick2end,
                                                          args.tick2d, args.mtick2, x2beg, x2end,
                                                          n2end - n2beg + 1, d2, size2, args.tick2format)

tick_3_location, tick_3_label, tick_3_minor = define_tick(args.ticks3, args.tick3beg, args.tick3end,
                                                          args.tick3d, args.mtick3, x3beg, x3end,
                                                          n3end - n3beg + 1, d3, size3, args.tick3format)

# if tick font size and family not speciefied, then inherit from axis labels
if args.tick1size is None:
    tick1size = label1size - 2
else:
    tick1size = float(args.tick1size)

if args.tick2size is None:
    tick2size = label2size - 2
else:
    tick2size = float(args.tick2size)

if args.tick3size is None:
    tick3size = label3size - 2
else:
    tick3size = float(args.tick3size)

ax11.yaxis.set_ticks(tick_2_location)
ax11.yaxis.set_ticklabels(tick_2_label,
                          fontsize=tick2size,
                          fontproperties=font,
                          rotation=float(args.tick2rot))

ax21.yaxis.set_ticks(tick_1_location)
ax21.yaxis.set_ticklabels(tick_1_label,
                          fontsize=tick1size,
                          fontproperties=font,
                          rotation=float(args.tick1rot))
ax21.xaxis.set_ticks(tick_3_location)
ax21.xaxis.set_ticklabels(tick_3_label,
                          fontsize=tick3size,
                          fontproperties=font,
                          rotation=float(args.tick3rot))

ax22.xaxis.set_ticks(tick_2_location)
ax22.xaxis.set_ticklabels(tick_2_label,
                          fontsize=tick2size,
                          fontproperties=font,
                          rotation=float(args.tick2rot))

# major and minor ticks sytle
tick_major_length = float(args.tickmajorlen)
tick_major_width = float(args.tickmajorwid)

ax11.tick_params('both', length=tick_major_length, width=tick_major_width, which='major')
ax21.tick_params('both', length=tick_major_length, width=tick_major_width, which='major')
ax22.tick_params('both', length=tick_major_length, width=tick_major_width, which='major')

# minor tick positions
ax11.yaxis.set_ticks(tick_2_minor, minor=True)
ax21.yaxis.set_ticks(tick_1_minor, minor=True)
ax21.xaxis.set_ticks(tick_3_minor, minor=True)
ax22.xaxis.set_ticks(tick_2_minor, minor=True)

# minor ticks style
if args.tickminorlen is None:
    tick_minor_length = 0.5 * tick_major_length
else:
    tick_minor_length = float(args.tickminorlen)

if args.tickminorwid is None:
    tick_minor_width = 0.75 * tick_major_width
else:
    tick_minor_width = float(args.tickminorwid)

ax11.tick_params('both', length=tick_minor_length, width=tick_minor_width, which='minor')
ax21.tick_params('both', length=tick_minor_length, width=tick_minor_width, which='minor')
ax22.tick_params('both', length=tick_minor_length, width=tick_minor_width, which='minor')

ax11.tick_params('x', labelbottom=0, bottom=0, labeltop=0, top=0, which='both')
ax11.tick_params('y', left=True, right=0, which='both')

ax21.tick_params('x', bottom=True, top=0, which='both')
ax21.tick_params('y', left=True, right=0, which='both')

ax22.tick_params('y', labelleft=0, left=0, right=0, labelright=0, which='both')
ax22.tick_params('x', bottom=True, top=0, which='both')

if args.tick1label == 0:
    ax21.yaxis.set_ticklabels([])
if args.tick2label == 0:
    ax11.yaxis.set_ticklabels([])
    ax22.xaxis.set_ticklabels([])
if args.tick3label == 0:
    ax21.xaxis.set_ticklabels([])

for l in ax11.yaxis.get_ticklabels():
    l.set_fontproperties(font)
    l.set_fontsize(tick2size)
for l in ax21.yaxis.get_ticklabels():
    l.set_fontproperties(font)
    l.set_fontsize(tick1size)
for l in ax21.xaxis.get_ticklabels():
    l.set_fontproperties(font)
    l.set_fontsize(tick3size)
for l in ax22.xaxis.get_ticklabels():
    l.set_fontproperties(font)
    l.set_fontsize(tick2size)

## curve annotation
from module_annotation import set_default
from matplotlib.patches import *
if args.curve is not None:

    curvefile = args.curve[0].split(",")
    nf = len(curvefile)

    curvestyle = set_default(args.curvestyle, ',', nf, 'scatter.')
    curvecolor = set_default(args.curvecolor, ',', nf, 'k')
    curvefacecolor = set_default(args.curvefacecolor, ',', nf, 'k')
    curveedgecolor = set_default(args.curveedgecolor, ',', nf, 'none')
    curvesize = set_default(args.curvesize, ',', nf, 1.0, 'float')
    curvewidth = set_default(args.curvewidth, ',', nf, 1.0, 'float')
    curveorder = set_default(args.curveorder, ',', nf, 9, 'int')

    for i in range(0, nf):

        curve = np.loadtxt(curvefile[i], ndmin=2)  # using ndmin=2 to ensure read as 2d array
        nsp = len(curve)
        x1 = curve[0:nsp, 0]
        x2 = curve[0:nsp, 1]
        x3 = curve[0:nsp, 2]
        px1 = (x1 - x1beg + 0.5 * d1) / (n1 * d1) * size1
        px2 = (x2 - x2beg + 0.5 * d2) / (n2 * d2) * size2
        px3 = (x3 - x3beg + 0.5 * d3) / (n3 * d3) * size3
        curve[0:nsp, 0] = px1
        curve[0:nsp, 1] = px2
        curve[0:nsp, 2] = px3

        # plot scatter points on the figure
        if 'scatter' in curvestyle[i]:

            ax11.scatter(px3,
                         px2,
                         marker=curvestyle[i][7:],
                         facecolor=curvefacecolor[i],
                         zorder=curveorder[i],
                         edgecolor=curveedgecolor[i],
                         s=curvesize[i])
            ax11.set_xlim(0, size3)
            ax11.set_ylim(0, size2)

            ax21.scatter(px3,
                         px1,
                         marker=curvestyle[i][7:],
                         facecolor=curvefacecolor[i],
                         zorder=curveorder[i],
                         edgecolor=curveedgecolor[i],
                         s=curvesize[i])
            ax21.set_xlim(0, size3)
            ax21.set_ylim(size1, 0)

            ax22.scatter(px2,
                         px1,
                         marker=curvestyle[i][7:],
                         facecolor=curvefacecolor[i],
                         zorder=curveorder[i],
                         edgecolor=curveedgecolor[i],
                         s=curvesize[i])
            ax22.set_xlim(0, size2)
            ax22.set_ylim(size1, 0)

        # plot line on the figure
        if 'line' in curvestyle[i]:

            extra = Line2D(px3,
                           px2,
                           linestyle=curvestyle[i][4:],
                           zorder=curveorder[i],
                           color=curvecolor[i],
                           linewidth=curvewidth[i])
            ax11.add_artist(extra)

            extra = Line2D(px3,
                           px1,
                           linestyle=curvestyle[i][4:],
                           zorder=curveorder[i],
                           color=curvecolor[i],
                           linewidth=curvewidth[i])
            ax21.add_artist(extra)

            extra = Line2D(px2,
                           px1,
                           linestyle=curvestyle[i][4:],
                           zorder=curveorder[i],
                           color=curvecolor[i],
                           linewidth=curvewidth[i])
            ax22.add_artist(extra)

        # plot polygon on the figure
        if 'polygon' in curvestyle[i]:

            curve = np.delete(curve, 2, 1)

            curve[0:nsp, 0] = px3
            curve[0:nsp, 1] = px2
            extra = Polygon(curve,
                            fill=False,
                            zorder=curveorder[i],
                            edgecolor=curvecolor[i],
                            linewidth=curvewidth[i])
            ax11.add_artist(extra)

            curve[0:nsp, 0] = px3
            curve[0:nsp, 1] = px1
            extra = Polygon(curve,
                            fill=False,
                            zorder=curveorder[i],
                            edgecolor=curvecolor[i],
                            linewidth=curvewidth[i])
            ax21.add_artist(extra)

            curve[0:nsp, 0] = px2
            curve[0:nsp, 1] = px1
            extra = Polygon(curve,
                            fill=False,
                            zorder=curveorder[i],
                            edgecolor=curvecolor[i],
                            linewidth=curvewidth[i])
            ax22.add_artist(extra)

## set grid line
if args.grid1:
    # grid line width
    if args.grid1width is None:
        grid1width = float(args.tickmajorwid)
    else:
        grid1width = float(args.grid1width)
    # add grid
    ax21.grid(which='major', axis='x', linestyle=args.grid1style, color=args.grid1color, linewidth=grid1width)
    ax22.grid(which='major', axis='x', linestyle=args.grid1style, color=args.grid1color, linewidth=grid1width)

if args.grid2:
    # grid line width
    if args.grid2width is None:
        grid2width = float(args.tickmajorwid)
    else:
        grid2width = float(args.grid2width)
    # add grid
    ax11.grid(which='major', axis='x', linestyle=args.grid2style, color=args.grid2color, linewidth=grid2width)
    ax22.grid(which='major', axis='y', linestyle=args.grid2style, color=args.grid2color, linewidth=grid2width)

if args.grid3:
    # grid line width
    if args.grid3width is None:
        grid3width = float(args.tickmajorwid)
    else:
        grid3width = float(args.grid3width)
    # add grid
    ax11.grid(which='major', axis='y', linestyle=args.grid3style, color=args.grid3color, linewidth=grid3width)
    ax21.grid(which='major', axis='y', linestyle=args.grid3style, color=args.grid3color, linewidth=grid3width)

## set title
if args.title is not None:

    if args.titlesize is None:
        title_font_size = max(float(args.label1size), float(args.label2size), float(args.label3size)) + 2
    else:
        title_font_size = float(args.titlesize)

    if args.titlex is None:
        title_x = 0.5 * (size2 + size3)
    else:
        title_x = float(args.titlex)

    if args.titley is None:
        title_y = size2 + 0.2
    else:
        title_y = float(args.titley)

    t = ax11.text(title_x,
                  title_y,
                  args.title,
                  horizontalalignment='center',
                  fontproperties=fontbold,
                  fontweight='bold')
    t.set_fontsize(title_font_size)

## set colorbar
if args.background is not None:
    cmin, cmax = set_clip(args, backdata, 'back', backdmin, backdmax)
if (args.legend and cmin != cmax) and (args.contourfill or args.overlay or args.background is not None):

    lloc = args.lloc

    # legend location
    if lloc == 'left':
        lorient = 'vertical'
        lleft = 'on'
        lright = 0
        ltop = 0
        lbottom = 0
        lrotate = 90
    if lloc == 'right':
        lorient = 'vertical'
        lleft = 0
        lright = 'on'
        ltop = 0
        lbottom = 0
        lrotate = 270
    if lloc == 'top':
        lorient = 'horizontal'
        lleft = 0
        lright = 0
        ltop = 'on'
        lbottom = 0
        lrotate = 0
    if lloc == 'bottom':
        lorient = 'horizontal'
        lleft = 0
        lright = 0
        ltop = 0
        lbottom = 'on'
        lrotate = 0

    if args.unitpad is None:
        if lloc == 'right':
            upad = 30.0
        if lloc == 'bottom':
            upad = 5.0
        if lloc == 'left':
            upad = 10.0
        if lloc == 'top':
            upad = 10.0
    else:
        upad = float(args.unitpad)

    # legend pad
    if args.lpad is None:
        lpad = 0.1
    else:
        lpad = float(args.lpad)

    # set colorbar axis location
    if lloc in ['left', 'right']:

        if args.lheight is None:
            lheight = size1 + size2 + float(args.slicegap)
        else:
            lheight = float(args.lheight)

        if args.lwidth is None:
            lwidth = 0.2
        else:
            lwidth = float(args.lwidth)

        if lloc == 'right':
            cbx = ax22locx + im22width + lpad / figwidth
            cby = ax22locy + (figheight - lheight) / 2.0 / figheight
        else:
            tlen = 0
            for i in tick_1_label:
                if len(i) > tlen: tlen = len(i)
            for i in tick_2_label:
                if len(i) > tlen: tlen = len(i)
            tlen = tlen + 1
            cbx = ax11locx - (max(label1size, label2size) * ipp + tlen * max(tick1size, tick2size) * ipp +
                              tick_major_length * ipp + lpad) / figwidth
            cby = ax21locy + (figheight - lheight) / 2.0 / figheight

    if args.lloc in ['top', 'bottom']:

        if args.lheight is None:
            lheight = 0.2
        else:
            lheight = float(args.lheight)

        if args.lwidth is None:
            lwidth = size2 + size3 + float(args.slicegap)
        else:
            lwidth = float(args.lwidth)

        if lloc == 'top':
            cbx = ax11locx + (figwidth - lwidth) / 2.0 / figwidth
            cby = ax11locy + im11height + lpad / figheight
        else:
            cbx = ax21locx + (figwidth - lwidth) / 2.0 / figwidth
            cby = ax21locy - (max(label2size, label3size) * ipp + 3.0 * max(tick2size, tick3size) * ipp +
                              tick_major_length * ipp + lpad) / figheight

    # add colorbar by add_axes
    cax = fig.add_axes([cbx, cby, lwidth / figwidth, lheight / figheight])
    cb = fig.colorbar(im11, cax=cax, orientation=lorient)

    # set colorbar label and styles
    if args.unitsize is None:
        lufs = min(float(args.label1size), float(args.label2size), float(args.label3size)) - 1
    else:
        lufs = float(args.unitsize)

    if args.unit is not None:
        if lloc == 'left' or lloc == 'right':
            cb.ax.set_ylabel(legend_units, rotation=lrotate, labelpad=upad)
        if lloc == 'top' or lloc == 'bottom':
            cb.ax.set_xlabel(legend_units, rotation=lrotate, labelpad=upad)

    cb.ax.yaxis.label.set_fontproperties(font)
    cb.ax.yaxis.label.set_fontsize(lufs)
    cb.ax.xaxis.label.set_fontproperties(font)
    cb.ax.xaxis.label.set_fontsize(lufs)

    # tick font size
    if args.lticksize is None:
        ltfs = lufs - 1
    else:
        ltfs = float(args.lticksize)

    if args.norm == 'linear':

        # set colorbar major ticks
        if args.ld is None:
            ld = nice((cmax - cmin) / 5.0)
        else:
            ld = float(args.ld)

        if args.ltickbeg is None:
            ltickbeg = nice(cmin, 0.5)
            base = 0.5
            nb = 0
            while nb <= 10 and ltickbeg > cmin + ld:
                base = base / 10.0
                ltickbeg = nice(cmin, base)
                nb = nb + 1
            if abs(ltickbeg) < abs(cmax) and orderm(ltickbeg) + 2 < orderm(cmax):
                ltickbeg = 0.0
        else:
            ltickbeg = float(args.ltickbeg)
        if args.ltickend is None:
            ltickend = cmax
        else:
            ltickend = float(args.ltickend)

        # scalar
        maxtick = max(abs(ltickbeg), abs(ltickend))
        if maxtick >= 1.0e4 or maxtick <= 1.0e-3:
            scalar = int(floor(log10(maxtick)))
            cscale = pow(10, scalar)
        else:
            cscale = 1.0

        ticks = np.arange(ltickbeg, ltickend + ld, ld)
        pminv = cmin
        pmaxv = cmax

        tbeg = max(pminv, ltickbeg)
        tend = min(pmaxv, ltickend)

        # set tick positions on colorbar
        ticks = np.asarray(
            [i for i in ticks if i >= tbeg - 1.0e-10 * abs(tbeg) and i <= tend + 1.0e-10 * abs(tend)])
        if lloc == 'left' or lloc == 'right':
            cb.ax.yaxis.set_ticks(ticks)
            cb.ax.yaxis.set_ticks_position(lloc)
        else:
            cb.ax.xaxis.set_ticks(ticks)
            cb.ax.xaxis.set_ticks_position(lloc)

        # set tick labels on colorbar
        tick_labels = ['' for i in range(0, len(ticks))]
        for i in range(0, len(ticks)):
            tick_labels[i] = ('%f' % (ticks[i] / cscale)).rstrip('0').rstrip('.')
        if lloc == 'left' or lloc == 'right':
            if cscale != 1:
                if lloc == 'left':
                    tick_labels[-1] = r'$\mathregular{10^{%i}}\times$' % scalar + '\n' + tick_labels[-1]
                else:
                    tick_labels[-1] = r'$\mathregular{\times 10^{%i}}$' % scalar + '\n' + tick_labels[-1]
            cb.ax.set_yticklabels(tick_labels)
        else:
            if cscale != 1:
                if lloc == 'top':
                    tick_labels[-1] = r'$\mathregular{\times 10^{%i}}$' % scalar + '\n' + tick_labels[-1]
                else:
                    tick_labels[-1] = tick_labels[-1] + '\n' + r'$\mathregular{\times 10^{%i}}$' % scalar
            cb.ax.set_xticklabels(tick_labels)

        # colorbar minor ticks
        if args.lmtick != 0:
            # extend tail and head
            pticks = np.append(ticks, ticks[0] - ld)
            pticks = np.append(pticks, ticks[-1] + ld)
            # sort all major ticks
            pticks = np.sort(pticks)
            # get pseudo-location of minor ticks
            nt = len(pticks)
            mticks = []
            for i in range(0, nt - 1):
                mticks = np.append(mticks, np.linspace(pticks[i], pticks[i + 1], args.lmtick + 2))
            mticks = [i for i in mticks if (i not in pticks)]
            mticks = np.asarray(
                [i for i in mticks if i >= tbeg - 1.0e-10 * abs(tbeg) and i <= tend + 1.0e-10 * abs(tend)])
            # set minor ticks
            if lloc == 'left' or lloc == 'right':
                cb.ax.yaxis.set_ticks(mticks, minor=True)
            else:
                cb.ax.xaxis.set_ticks(mticks, minor=True)

    if args.norm == 'log':

        # set colorbar major ticks
        if args.ltickbeg is None:
            ltickbeg = np.floor(cmin)
        else:
            ltickbeg = float(args.ltickbeg)
        if args.ltickend is None:
            ltickend = np.ceil(cmax)
        else:
            ltickend = float(args.ltickend)
        if args.ld is None:
            ld = max(1, round((ltickend - ltickbeg) / 5.0))
        else:
            ld = int(args.ld)

        ticks = np.arange(ltickbeg, ltickend + 1, ld)
        pminv = cmin
        pmaxv = cmax

        tbeg = max(pminv, ltickbeg)
        tend = min(pmaxv, ltickend)

        # set tick positions on colorbar
        ticks = np.asarray(
            [i for i in ticks if i >= tbeg - 1.0e-10 * abs(tbeg) and i <= tend + 1.0e-10 * abs(tend)])

        # set tick positions on colorbar
        if lloc == 'left' or lloc == 'right':
            cb.ax.yaxis.set_ticks(ticks)
            cb.ax.yaxis.set_ticks_position(lloc)
        else:
            cb.ax.xaxis.set_ticks(ticks)
            cb.ax.xaxis.set_ticks_position(lloc)

        # set tick labels on colorbar
        tick_labels = ['' for i in range(0, len(ticks))]
        for i in range(0, len(ticks)):
            tick_labels[i] = '$\mathregular{10^{%i}}$' % (ticks[i])
        if lloc == 'left' or lloc == 'right':
            cb.ax.set_yticklabels(tick_labels)
        else:
            cb.ax.set_xticklabels(tick_labels)

        # colorbar minor ticks
        if args.lmtick != 0:
            # extend tail and head
            pticks = np.append(ticks, ticks[0] - ld)
            pticks = np.append(pticks, ticks[-1] + ld)
            # sort all major ticks
            pticks = np.sort(pticks)
            # get pseudo-location of minor ticks
            nt = len(pticks)
            mticks = []
            for i in range(0, nt - 1):
                mticks = np.append(mticks,
                                   np.log10(np.linspace(10**pticks[i], 10**pticks[i + 1], args.lmtick + 2)))
            mticks = np.asarray(
                [i for i in mticks if i >= tbeg - 1.0e-10 * abs(tbeg) and i <= tend + 1.0e-10 * abs(tend)])
            # set minor ticks
            if lloc == 'left' or lloc == 'right':
                cb.ax.yaxis.set_ticks(mticks, minor=True)
            else:
                cb.ax.xaxis.set_ticks(mticks, minor=True)

    cb.ax.tick_params(
        direction='out',
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        top=ltop,  # ticks along the top axis
        bottom=lbottom,  # ticks along the bottom axis
        labeltop=ltop,  # labels along the top axis
        labelbottom=lbottom)  # labels along the bottom axis
    cb.ax.tick_params(
        direction='out',
        axis='y',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        left=lleft,  # ticks along the left axis
        right=lright,  # ticks along the right axis
        labelleft=lleft,  # labels along the left axis
        labelright=lright)  # labels along the right axis

    # set colorbar tick styles: font size, family, and ticks direction
    if lloc == 'left' or lloc == 'right':
        for l in cb.ax.yaxis.get_ticklabels():
            l.set_fontproperties(font)
            l.set_fontsize(ltfs)
    else:
        for l in cb.ax.xaxis.get_ticklabels():
            l.set_fontproperties(font)
            l.set_fontsize(ltfs)

    if lloc == 'left' or lloc == 'right':
        cb.ax.yaxis.set_label_position(lloc)
    if lloc == 'top' or lloc == 'bottom':
        cb.ax.xaxis.set_label_position(lloc)

    # make colorbar solid color continuous
    cb.solids.set_edgecolor("face")

    # colorbar reverse
    if lloc in ['left', 'right']:
        if args.lreverse == 1:
            cb.ax.invert_yaxis()
    else:
        if args.lreverse == 1:
            cb.ax.invert_xaxis()

## set slice line
scale1 = size1 / ((n1end - n1beg) * d1)
scale2 = size2 / ((n2end - n2beg) * d2)
scale3 = size3 / ((n3end - n3beg) * d3)

# plot slice lines
if args.sliceline1:
    sl1pos = sl1 * scale1
    linex = [0, size3]
    linez = [sl1pos, sl1pos]
    extra = Line2D(linex,
                   linez,
                   linestyle=args.sliceline1style,
                   zorder=3,
                   color=args.sliceline1color,
                   lw=float(args.sliceline1width))
    ax21.add_artist(extra)
    linex = [0, size2]
    linez = [sl1pos, sl1pos]
    extra = Line2D(linex,
                   linez,
                   linestyle=args.sliceline1style,
                   zorder=3,
                   color=args.sliceline1color,
                   lw=float(args.sliceline1width))
    ax22.add_artist(extra)

if args.sliceline2:
    sl2pos = sl2 * scale2
    linex = [0, size3]
    linez = [sl2pos, sl2pos]
    extra = Line2D(linex,
                   linez,
                   linestyle=args.sliceline2style,
                   zorder=3,
                   color=args.sliceline2color,
                   lw=float(args.sliceline2width))
    ax11.add_artist(extra)
    linex = [sl2pos, sl2pos]
    linez = [0, size1]
    extra = Line2D(linex,
                   linez,
                   linestyle=args.sliceline2style,
                   zorder=3,
                   color=args.sliceline2color,
                   lw=float(args.sliceline2width))
    ax22.add_artist(extra)

if args.sliceline3:
    sl3pos = sl3 * scale3
    linex = [sl3pos, sl3pos]
    linez = [0, size1]
    extra = Line2D(linex,
                   linez,
                   linestyle=args.sliceline3style,
                   zorder=3,
                   color=args.sliceline3color,
                   lw=float(args.sliceline3width))
    ax21.add_artist(extra)
    linex = [sl3pos, sl3pos]
    linez = [0, size2]
    extra = Line2D(linex,
                   linez,
                   linestyle=args.sliceline3style,
                   zorder=3,
                   color=args.sliceline3color,
                   lw=float(args.sliceline3width))
    ax11.add_artist(extra)

## output
output(args)

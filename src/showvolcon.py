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
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib import colors
from matplotlib.colors import colorConverter
from argparse import RawTextHelpFormatter
from module_io import *

# ignore warnings
import warnings
warnings.filterwarnings("ignore", module="matplotlib")

# tag
program = 'volcon'
print()

# read arguments
parser = argparse.ArgumentParser(description='''
                                purpose:
                                    Plot a 3D array as a full or cropped volume of contours.
                                ''',
                                 formatter_class=RawTextHelpFormatter)
parser = getarg(parser, program)
args = parser.parse_args()

# input
data, n1, n2, n3, dmin, dmax = read_array(args, which='fore', dim=3)
mask, _, _, _, _, _ = read_array(args, which='mask', dim=3)
if mask is not None:
    data = data * mask

d1 = float(args.d1)
d2 = float(args.d2)
d3 = float(args.d3)

## limit of axis
from module_range import *
sp1beg, sp1end, x1beg, x1end, n1beg, n1end = set_range(args.o1, n1, d1, args.x1beg, args.x1end)
sp2beg, sp2end, x2beg, x2end, n2beg, n2end = set_range(args.o2, n2, d2, args.x2beg, args.x2end)
sp3beg, sp3end, x3beg, x3end, n3beg, n3end = set_range(args.o3, n3, d3, args.x3beg, args.x3end)

n1 = n1end - n1beg + 1
n2 = n2end - n2beg + 1
n3 = n3end - n3beg + 1

## set slice
# axis 1
if args.slice1 is None:
    sl1 = (x1end + x1beg) / 2.0
else:
    sl1 = eval(args.slice1)

# axis 2
if args.slice2 is None:
    sl2 = (x2end + x2beg) / 2.0
else:
    sl2 = eval(args.slice2)

# axis 3
if args.slice3 is None:
    sl3 = (x3end + x3beg) / 2.0
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

# recalculate slice position
sl1 = x1beg + (slice1 - n1beg) * d1
sl2 = x2beg + (slice2 - n2beg) * d2
sl3 = x3beg + (slice3 - n3beg) * d3

takeslice = (args.slice1) != 0 or len(args.slice2) != 0 or len(args.slice3 is not None)

# select data based on angle
if args.angle is None:
    angle1 = 40.0
    angle2 = 10.0
else:
    angles = args.angle[0].split(',')
    if len(angles) == 1:
        angle1 = float(args.angle[0])
        angle2 = 0.0
    if len(angles) == 2:
        angle1 = float(angles[0])
        angle2 = float(angles[1])

# check angle ranges
octant = args.octant
if not (octant in ['--+', '---', '-+-', '-++']):
    print('error: octant should be one of --+, ---, -+-, -++')
    exit()
if angle1 < 0 or angle1 > 90 or angle2 < 0 or angle2 > 90:
    print('error: angles should be in [0,90]')
    exit()

# reassign angle
angles = [angle1, angle2]
octant_1 = (octant == '--+')
octant_2 = (octant == '---')
octant_3 = (octant == '-+-')
octant_4 = (octant == '-++')

if octant_1:  # --+

    # note the global index of slices
    # for different view angles, the indices are different

    # slice data
    # slice xy
    data12 = data[n1beg:slice1 + 1, n2beg:slice2 + 1, slice3]
    if args.norm == 'log': data12 = np.log10(data12)

    # slice xz
    data13 = data[n1beg:slice1 + 1, slice2, slice3:n3end]
    if args.norm == 'log': data13 = np.log10(data13)

    # slice yz
    data23 = data[slice1, n2beg:slice2 + 1, slice3:n3end]
    if args.norm == 'log': data23 = np.log10(data23)

    # face data
    # slice xy
    fdata12 = data[n1beg:n1end, n2beg:n2end, n3end - 1]
    if args.norm == 'log': fdata12 = np.log10(fdata12)

    # slice xz
    fdata13 = data[n1beg:n1end, n2beg, n3beg:n3end]
    if args.norm == 'log': fdata13 = np.log10(fdata13)

    # slice yz
    fdata23 = data[n1beg, n2beg:n2end, n3beg:n3end]
    if args.norm == 'log': fdata23 = np.log10(fdata23)

    # mask array with nan
    if takeslice:
        mask = np.zeros([n1, n2, n3])
        mask[n1beg:slice1, n2beg:slice2, slice3 + 1:n3end] = 1
        mask12 = mask[n1beg:n1end, n2beg:n2end, n3end - 1]
        mask13 = mask[n1beg:n1end, n2beg, n3beg:n3end]
        mask23 = mask[n1beg, n2beg:n2end, n3beg:n3end]
    else:
        mask = []
        mask12 = []
        mask13 = []
        mask23 = []

if octant_2:  # ---

    # slice data
    # slice xy
    data12 = data[n1beg:slice1 + 1, n2beg:slice2 + 1, slice3]
    if args.norm == 'log': data12 = np.log10(data12)

    # slice xz
    data13 = data[n1beg:slice1 + 1, slice2, n3beg:slice3 + 1]
    if args.norm == 'log': data13 = np.log10(data13)

    # slice yz
    data23 = data[slice1, n2beg:slice2 + 1, n3beg:slice3 + 1]
    if args.norm == 'log': data23 = np.log10(data23)

    # face data
    # slice xy
    fdata12 = data[n1beg:n1end, n2beg:n2end, n3beg]
    if args.norm == 'log': fdata12 = np.log10(fdata12)

    # slice xz
    fdata13 = data[n1beg:n1end, n2beg, n3beg:n3end]
    if args.norm == 'log': fdata13 = np.log10(fdata13)

    # slice yz
    fdata23 = data[n1beg, n2beg:n2end, n3beg:n3end]
    if args.norm == 'log': fdata23 = np.log10(fdata23)

    if takeslice:
        mask = np.zeros([n1, n2, n3])
        mask[n1beg:slice1, n2beg:slice2, n3beg:slice3] = 1
        mask12 = mask[n1beg:n1end, n2beg:n2end, n3beg]
        mask13 = mask[n1beg:n1end, n2beg, n3beg:n3end]
        mask23 = mask[n1beg, n2beg:n2end, n3beg:n3end]
    else:
        mask = []
        mask12 = []
        mask13 = []
        mask23 = []

if octant_3:  # -+-

    # slice data
    # slice xy
    data12 = data[n1beg:slice1 + 1, slice2:n2end, slice3]
    if args.norm == 'log': data12 = np.log10(data12)

    # slice xz
    data13 = data[n1beg:slice1 + 1, slice2, n3beg:slice3 + 1]
    if args.norm == 'log': data13 = np.log10(data13)

    # slice yz
    data23 = data[slice1, slice2:n2end, n3beg:slice3 + 1]
    if args.norm == 'log': data23 = np.log10(data23)

    # face data
    # slice xy
    fdata12 = data[n1beg:n1end, n2beg:n2end, n3beg]
    if args.norm == 'log': fdata12 = np.log10(fdata12)

    # slice xz
    fdata13 = data[n1beg:n1end, n2end - 1, n3beg:n3end]
    if args.norm == 'log': fdata13 = np.log10(fdata13)

    # slice yz
    fdata23 = data[n1beg, n2beg:n2end, n3beg:n3end]
    if args.norm == 'log': fdata23 = np.log10(fdata23)

    if takeslice:
        mask = np.zeros([n1, n2, n3])
        mask[n1beg:slice1, slice2 + 1:n2end, n3beg:slice3] = 1
        mask12 = mask[n1beg:n1end, n2beg:n2end, n3beg]
        mask13 = mask[n1beg:n1end, n2end - 1, n3beg:n3end]
        mask23 = mask[n1beg, n2beg:n2end, n3beg:n3end]
    else:
        mask = []
        mask12 = []
        mask13 = []
        mask23 = []

if octant_4:  # -++

    # slice data
    # slice xy
    data12 = data[n1beg:slice1 + 1, slice2:n2end, slice3]
    if args.norm == 'log': data12 = np.log10(data12)

    # slice xz
    data13 = data[n1beg:slice1 + 1, slice2, slice3:n3end]
    if args.norm == 'log': data13 = np.log10(data13)

    # slice yz
    data23 = data[slice1, slice2:n2end, slice3:n3end]
    if args.norm == 'log': data23 = np.log10(data23)

    # face data
    # slice xy
    fdata12 = data[n1beg:n1end, n2beg:n2end, n3end - 1]
    if args.norm == 'log': fdata12 = np.log10(fdata12)

    # slice xz
    fdata13 = data[n1beg:n1end, n2end - 1, n3beg:n3end]
    if args.norm == 'log': fdata13 = np.log10(fdata13)

    # slice yz
    fdata23 = data[n1beg, n2beg:n2end, n3beg:n3end]
    if args.norm == 'log': fdata23 = np.log10(fdata23)

    if takeslice:
        mask = zeros([n1, n2, n3])
        mask[n1beg:slice1, slice2 + 1:n2end, slice3 + 1:n3end] = 1
        mask12 = mask[n1beg:n1end, n2beg:n2end, n3end - 1]
        mask13 = mask[n1beg:n1end, n2end - 1, n3beg:n3end]
        mask23 = mask[n1beg, n2beg:n2end, n3beg:n3end]
    else:
        mask = []
        mask12 = []
        mask13 = []
        mask23 = []

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


## projection operation
# axonometric projection
def isometric_projection(x, y, z, angles):
    alpha = angles[0]
    beta = angles[1]
    bx = x * np.cos(alpha * np.pi / 180.0) - y * np.cos(beta * np.pi / 180.0)
    by = z + x * np.sin(alpha * np.pi / 180.0) + y * np.sin(beta * np.pi / 180.0)
    return bx, by


## plot contour and axis on deformed rectangular mesh
def project_contour(args, ax, data, px, py, colormap, cmin, cmax, font, mask=[]):

    # determine input dimension
    n1, n2 = data.shape

    # generate shape functions
    s1 = np.linspace(1, 0, n1)
    s2 = np.linspace(0, 1, n1)
    s3 = np.linspace(1, 0, n2)
    s4 = np.linspace(0, 1, n2)

    shapes = [[[s1[i] * s3[j], s2[i] * s3[j], s2[i] * s4[j], s1[i] * s4[j]] for j in range(0, n2)]
              for i in range(0, n1)]

    # triangulartion
    nnx = np.zeros([n1, n2])
    nny = np.zeros([n1, n2])
    for i in range(0, n1):
        for j in range(0, n2):
            nnx[i, j] = sum(shapes[i][j] * np.asarray(px))
            nny[i, j] = sum(shapes[i][j] * np.asarray(py))
    nnx = nnx.flatten()
    nny = nny.flatten()
    data = data.flatten()
    mesh = mplt.tri.Triangulation(nnx, nny)

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
                mlevels = np.append(
                    mlevels, np.log10(np.linspace(10**levels[i], 10**levels[i + 1], args.mcontour + 2)))
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
                tls = np.append(
                    tls, np.log10(np.linspace(10**levels[int(i)], 10**levels[int(i) + 1], mcontour + 1)))
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

    lw = np.asarray([mw for i in range(0, size(levels))])
    lw[0:-1:mcontour] = float(args.contourwidth)

    ls = np.array([args.mcontourstyle for i in range(0, size(levels))])
    ls[0:-1:mcontour] = args.contourstyle

    lc0 = args.contourcolor[0].split(',')
    lc = ['k' for i in range(0, size(levels))]
    lc[0:size(lc0)] = lc0

    # show filled contours if necessary
    # fill contours use full mesh since no need to worry about overlap
    if args.contourfill:
        if args.norm == 'linear':
            cf = ax.tricontourf(mesh,
                                data,
                                levels[0:size(levels) - 1],
                                cmap=colormap,
                                antialiased=True,
                                extend='both')
        if args.norm == 'log':
            from matplotlib.colors import LogNorm
            cf = ax.tricontourf(mesh,
                                data,
                                levels[0:size(levels) - 1],
                                cmap=colormap,
                                antialiased=True,
                                norm=LogNorm(cmin, cmax))
        for l in cf.collections:
            l.set_edgecolor('face')
            l.set_linewidth(0.025)

    # show ordinary contours by default
    # triangular mesh mask
    if len(mask) > 0:
        mask = mask.flatten()
        mask = mask[mesh.triangles].mean(axis=1)
        mask = np.where(mask > 0.33333333, 1, 0)
        mesh.set_mask(mask)

    # contours use masked mesh since otherwise overlap
    cs = ax.tricontour(mesh,
                       data,
                       levels[0:size(levels) - 1],
                       colors=lc,
                       linewidths=lw,
                       linestyles=ls,
                       antialiased=True)
    #        inline=True,

    if clabelsize != 0:

        # choose label levels
        lvl = cs.levels[::mcontour]

        # set format
        clabels = ['' for i in range(0, size(lvl))]
        if args.norm == 'linear':
            for i in range(0, size(lvl)):
                if lvl[i] != 0 and (abs(lvl[i]) < 1.0e-3 or abs(lvl[i]) > 1.0e3):
                    scalar = int(floor(log10(abs(lvl[i]))))
                    cscale = pow(10, scalar)
                    clabels[i] = ('%f' % (lvl[i] / cscale)
                                  ).rstrip('0').rstrip('.') + r'$\mathregular{\times 10^{%i}}$' % scalar
                else:
                    clabels[i] = ('%f' % (lvl[i])).rstrip('0').rstrip('.')

        if args.norm == 'log':
            for i in range(0, size(lvl)):
                clabels[i] = r'$\mathregular{10^{%i}}$' % (lvl[i])

        fmt = {}
        for l, s in zip(cs.levels[::mcontour], clabels):
            fmt[l] = s

        # place contour labels
        clabels = ax.clabel(cs, cs.levels[::mcontour], inline=True, fmt=fmt,
                            fontsize=clabelsize)  #, fontproperties=font)
        #        ax.get_label.set_fontproperties(font)
        #        ax.get_label.set_size(font)
        for txt in clabels:
            txt.set_fontproperties(font)
            txt.set_fontsize(clabelsize)
            txt.set_color(args.clabelcolor)
            if args.clabelbackcolor is not None:
                txt.set_backgroundcolor(args.clabelbackcolor)


## plot image and axis on deformed rectangular mesh
def project_image(ax, data, px, py, colormap, cmin, cmax):

    # determine input dimension
    n1, n2 = data.shape
    n1 = n1 + 1
    n2 = n2 + 1

    # generate shape functions
    s1 = np.linspace(1, 0, n1)
    s2 = np.linspace(0, 1, n1)
    s3 = np.linspace(1, 0, n2)
    s4 = np.linspace(0, 1, n2)

    shapes = [[[s1[i] * s3[j], s2[i] * s3[j], s2[i] * s4[j], s1[i] * s4[j]] for j in range(0, n2)]
              for i in range(0, n1)]

    # do projection
    nnx = [[sum(shapes[i][j] * np.asarray(px)) for j in range(0, n2)] for i in range(0, n1)]
    nny = [[sum(shapes[i][j] * np.asarray(py)) for j in range(0, n2)] for i in range(0, n1)]

    # pcolor fast is the fastest way, using polygon collections is slow
    im = ax.pcolorfast(nnx, nny, data, linewidth=0, edgecolor='none', antialiased=True)
    im.set_cmap(colormap)
    im.set_clim(cmin, cmax)
    im.set_rasterized(True)


from module_utility import *


def project_axis(ax, p1x, p1y, p2x, p2y, ticks, tickbeg, tickend, tickd, mtick, xbeg, xend, ns, d, font,
                 tick_major_len, tick_major_width, tick_minor_len, tick_minor_width, axislen, tick_orient,
                 ticklabel_orient, tick_size, label, label_orient, label_size, label_pad, tick_format):

    # regular ticks
    if ticks is None:

        # major tick interval
        if tickd is None:
            tick_interval = nice((xend - xbeg) / 5.0)
            if tick_interval == 0:
                tick_interval = 1.0e10
        else:
            tick_interval = float(tickd)

        # tick begin location
        if tickbeg is None:
            tick_beg = nice(xbeg)
            base = 0.5
            nb = 0
            if tick_interval > 0:
                while nb <= 10 and tick_beg > xbeg + tick_interval:
                    base = base / 10.0
                    tick_beg = nice(xbeg, base)
                    nb = nb + 1
            else:
                while nb <= 10 and tick_beg < xbeg + tick_interval:
                    base = base / 10.0
                    tick_beg = nice(xbeg, base)
                    nb = nb + 1
        else:
            tick_beg = float(tickbeg)

        # tick end location
        if tickend is None:
            tick_end = tick_beg + (round((xend - xbeg) / tick_interval) + 2) * tick_interval
            if tick_interval > 0:
                while tick_end < xend:
                    tick_end = tick_end + abs(tick_interval)
            else:
                while tick_end > xend:
                    tick_end = tick_end - abs(tick_interval)
        else:
            tick_end = float(tickend)

        # regular major and minor tick locations
        tick = np.arange(tick_beg, tick_end + 0.1 * abs(tick_interval), tick_interval)
        minor_tick_interval = tick_interval / (mtick + 1.0)
        minor_tick = np.arange(tick_beg, tick_end + 0.1 * abs(minor_tick_interval), minor_tick_interval)

        # some ticks might out of axis range
        if d > 0:
            tick = np.asarray([i for i in tick if i >= xbeg and i <= xend])
            minor_tick = np.asarray([i for i in minor_tick if i >= xbeg and i <= xend and (not i in tick)])
        if d < 0:
            tick = np.asarray([i for i in tick if i <= xbeg and i >= xend])
            minor_tick = np.asarray([i for i in minor_tick if i <= xbeg and i >= xend and (not i in tick)])

        # linearly scale the ticks to figure canvas
        if ns == 1:
            print('error: contour plot requires at least two points')
        if ns == 2:
            # if only one sample point, then tick location is 0.5
            tick_location = np.asarray([0, axislen])
            ntick = 2
        if ns > 2:
            # if multiple sample points, then scale to apparent axis length
            # !! for contour plot remove the 0.5*d extra length
            # tick_location=(tick-xbeg+0.5*d)/((ns-1)*d)*axislen
            tick_location = (tick - xbeg) / ((ns - 2) * d) * axislen
            minor_tick_location = (minor_tick - xbeg) / ((ns - 2) * d) * axislen
            t = tick_location

        # set major tick location and labels, note some major ticks might be out of axis range
        tl = []
        tick_label = []
        for i in range(0, len(tick)):
            if tick_location[i] >= 0 and tick_location[i] <= axislen + 1.0e-10:
                tl.append(tick_location[i])
                if tick_format == 'sci' or tick_format == 'plain':
                    tick_label.append(('%f' % tick[i]).rstrip('0').rstrip('.'))
                else:
                    tick_label.append((tick_format % tick[i]))
        tick_location = tl

    # irregular ticks
    else:

        # get contents from user-specified ticks
        ticks = ticks[0].split(':')
        location = [0 for i in range(0, len(ticks))]
        label = ['' for i in range(0, len(ticks))]

        # set tick locations
        for i in range(0, len(ticks)):
            t = ticks[i].split(',')
            # !!
            # location[i]=(float(t[0])+0.5*d)/((ns-1)*d)*axislen
            location[i] = (float(t[0])) / ((ns - 2) * d) * axislen
            label[i] = t[1]

        # sort according to tick location
        yx = list(zip(location, label))
        yx.sort()
        tick_location = [location for location, label in yx]
        tick_label = [label for location, label in yx]

        # minor ticks
        mtick = mtick + 1
        minor_tick_location = np.linspace(tick_location[0], tick_location[1], mtick + 1)
        minor_tick_location = minor_tick_location[1:mtick]
        for i in range(1, len(tick_location) - 1):
            t = np.linspace(tick_location[i], tick_location[i + 1], mtick + 1)
            minor_tick_location = np.append(minor_tick_location, t[1:mtick])

    # major ticks
    # projection of ticks to target position
    # !! ensure no zero encountered
    if tick_location[-1] == axislen:
        tick_location[-1] = axislen - 1.0e-10
    tx = [p1x + i * (p2x - p1x) / axislen for i in tick_location]
    ty = [p1y + i * (p2y - p1y) / axislen for i in tick_location]

    # in the following, the ticks are assumed to be always perpendicular with
    # the axis, whatever the pointing of the axis is
    if tick_orient == 'counterclock':
        # outward pointing -- 90 degree counterclockwise
        # vector is (x1,y1) ----> (x2,y2)
        # outward means the tick should counterclockwise rotate this vector to
        # get ticks
        tx2 = -(p2y - ty)
        ty2 = p2x - tx

    if tick_orient == 'clock':
        # inward pointing -- 90 degree clockwise
        # vector is (x1,y1) ----> (x2,y2)
        # inward means the tick should clockwise rotate this vector to
        # get ticks
        tx2 = p2y - ty
        ty2 = -(p2x - tx)

    # direct = will cause txx/tyy always equal to tx2/ty2
    txx = [i for i in tx2]
    tyy = [i for i in ty2]

    ticklen = tick_major_len * ipp
    for i in range(0, len(tick_location)):
        epx = tx[i] + tx2[i] / sqrt(tx2[i]**2 + ty2[i]**2) * ticklen
        epy = ty[i] + ty2[i] / sqrt(tx2[i]**2 + ty2[i]**2) * ticklen
        tx2[i] = epx
        ty2[i] = epy

    # form lines
    majortick = []
    e1 = list(zip(tx, ty))
    e2 = list(zip(tx2, ty2))
    for i in range(0, len(tick_location)):
        majortick.append([e1[i], e2[i]])

    majortick = LineCollection(majortick, linewidths=tick_major_width, colors='k')
    ax.add_collection(majortick)

    # get center location of the axis
    hx = 0.5 * (p1x + p2x)
    hy = 0.5 * (p1y + p2y)

    # get normal vector starting from center location of the axis
    if tick_orient == 'counterclock':
        # outward pointing -- 90 degree counterclockwise
        hx2 = -(p2y - hy)
        hy2 = p2x - hx
    if tick_orient == 'clock':
        # inward pointing -- 90 degree clockwise
        hx2 = p2y - hy
        hy2 = -(p2x - hx)

    # angle of axis [-180,180]
    axis_angle = np.arctan2(p2y - p1y, p2x - p1x + 1.0e-10) * 180.0 / np.pi

    # angle of tick label
    if ticklabel_orient == 'positive':
        ticklabel_angle = axis_angle
    if ticklabel_orient == 'negative':
        ticklabel_angle = axis_angle - 180.0
    if ticklabel_orient == 'counterclock':
        ticklabel_angle = axis_angle + 90.0
    if ticklabel_orient == 'clock':
        ticklabel_angle = axis_angle - 90.0

    # maximum tick label length
    ticklabellen = [0 for i in range(0, len(tick_location))]
    for i in range(0, len(tick_location)):
        ticklines = tick_label[i].split('\n')
        if ticklabel_orient == 'positive' or ticklabel_orient == 'negative':
            ticklabellen[i] = len(ticklines) * tick_size * ipp
        if ticklabel_orient == 'counterclock' or ticklabel_orient == 'clock':
            ticklabellen[i] = 0.6 * max([len(j) for j in ticklines]) * tick_size * ipp
    maxlen = max(ticklabellen)

    # tick label alignment
    pad = ticklen + 0.1
    if ticklabel_orient == 'positive' or ticklabel_orient == 'negative':
        ha = 'center'
        va = 'center'
        pad = pad + 0.5 * tick_size * ipp
    if ticklabel_orient == 'clock':
        if tick_orient == 'clock':
            ha = 'left'
            va = 'center'
        if tick_orient == 'counterclock':
            ha = 'right'
            va = 'center'
    if ticklabel_orient == 'counterclock':
        if tick_orient == 'clock':
            ha = 'right'
            va = 'center'
        if tick_orient == 'counterclock':
            ha = 'left'
            va = 'center'

    # tick labels
    for i in range(0, len(tick_location)):
        epx = tx[i] + txx[i] / sqrt(txx[i]**2 + tyy[i]**2) * pad
        epy = ty[i] + tyy[i] / sqrt(txx[i]**2 + tyy[i]**2) * pad
        ax.text(epx,
                epy,
                tick_label[i],
                ha=ha,
                va=va,
                rotation=ticklabel_angle,
                fontproperties=font,
                size=tick_size)

    # angle of axis label
    if label_orient == 'positive':
        label_angle = axis_angle
    if label_orient == 'negative':
        label_angle = axis_angle - 180.0
    if label_orient == 'counterclock':
        label_angle = axis_angle + 90.0
    if label_orient == 'clock':
        label_angle = axis_angle - 90.0

    # axis label alignment
    pad = ticklen + 0.05 + maxlen + 0.1
    if label_orient == 'positive' or label_orient == 'negative':
        ha = 'center'
        va = 'center'
        pad = pad + 0.6 * label_size * ipp + label_pad
    if label_orient == 'clock':
        ha = 'left'
        va = 'center'
        pad = pad + label_pad
    if label_orient == 'counterclock':
        ha = 'right'
        va = 'center'
        pad = pad + label_pad

    # maximum axis label length
    labellines = label.split('\n')
    if ticklabel_orient == 'positive' or ticklabel_orient == 'negative':
        labellen = len(labellines) * label_size * ipp
    if ticklabel_orient == 'counterclock' or ticklabel_orient == 'clock':
        labellen = 0.6 * max([len(i) for i in labellines]) * label_size * ipp

    # axis label position
    labelx = hx + hx2 / sqrt(hx2**2 + hy2**2) * pad
    labely = hy + hy2 / sqrt(hx2**2 + hy2**2) * pad
    if label is not None:
        ax.text(labelx,
                labely,
                label,
                ha=ha,
                va=va,
                rotation=label_angle,
                fontproperties=font,
                size=label_size)

    # minor tick locations
    # projection of ticks to target position
    tx = [p1x + i * (p2x - p1x) / axislen for i in minor_tick_location]
    ty = [p1y + i * (p2y - p1y) / axislen for i in minor_tick_location]

    if tick_orient == 'counterclock':
        # outward pointing -- 90 degree counterclockwise
        tx2 = -(p2y - ty)
        ty2 = p2x - tx

    if tick_orient == 'clock':
        # inward pointing -- 90 degree clockwise
        tx2 = p2y - ty
        ty2 = -(p2x - tx)

    ticklen = tick_minor_len * ipp
    for i in range(0, len(minor_tick_location)):
        epx = tx[i] + tx2[i] / sqrt(tx2[i]**2 + ty2[i]**2) * ticklen
        epy = ty[i] + ty2[i] / sqrt(tx2[i]**2 + ty2[i]**2) * ticklen
        tx2[i] = epx
        ty2[i] = epy

    # form lines
    minortick = []
    e1 = list(zip(tx, ty))
    e2 = list(zip(tx2, ty2))
    for i in range(0, len(minor_tick_location)):
        minortick.append([e1[i], e2[i]])

    minortick = LineCollection(minortick, linewidths=tick_minor_width, colors='k')
    ax.add_collection(minortick)


## color converter
def cc(arg):
    return colorConverter.to_rgba(arg)


## compose volume
if octant_1:

    angle1 = angle1
    angles = [angle1, angle2]

    #
    #    A sketch of cavalier projection. 
    #    The relative point positions also apply to isometric projection
    #
    #
    #         8 *------------------------* 6
    #          /.                       /|
    #         / .                      / |
    #        /  .        #............#  |
    #       /   .       .:            :  |
    #      /    .      . :            :  |
    #   2 /     *(7) .  :      (4)    :  * 5 (size2)
    #    *-----------#   #............# /
    #    | size1     : .            .  /
    #    |           :.            .  /
    #    |           #............#  /
    #    |                        | /
    #    |                        |/ ) angle
    #    *------------------------*  ----
    #    1                        3 (size3)
    #

    # vertices positions for faces
    p1x, p1y = isometric_projection(0, size3, 0, angles)
    p2x, p2y = isometric_projection(0, size3, size1, angles)
    p3x, p3y = isometric_projection(0, 0, 0, angles)
    p4x, p4y = isometric_projection(0, 0, size1, angles)
    p5x, p5y = isometric_projection(size2, 0, 0, angles)
    p6x, p6y = isometric_projection(size2, 0, size1, angles)
    p7x, p7y = isometric_projection(size2, size3, 0, angles)
    p8x, p8y = isometric_projection(size2, size3, size1, angles)

    # set slice cut plane
    if x1end != x1beg:
        scale1 = size1 / (x1end - x1beg)
    else:
        scale1 = 1.0e2
    if x2end != x2beg:
        scale2 = size2 / (x2end - x2beg)
    else:
        scale2 = 1.0e2
    if x3end != x3beg:
        scale3 = size3 / (x3end - x3beg)
    else:
        scale3 = 1.0e2
    sl1 = size1 - (sl1 - x1beg) / ((n1end - n1beg - 1) * d1) * size1
    sl2 = (sl2 - x2beg) / ((n2end - n2beg - 1) * d2) * size2
    sl3 = size3 - (sl3 - x3beg) / ((n3end - n3beg - 1) * d3) * size3

    # vertices positions for slices
    q1x, q1y = isometric_projection(0, sl3, sl1, angles)
    q2x, q2y = isometric_projection(0, sl3, size1, angles)
    q3x, q3y = isometric_projection(0, 0, sl1, angles)
    q4x, q4y = isometric_projection(0, 0, size1, angles)
    q5x, q5y = isometric_projection(sl2, 0, sl1, angles)
    q6x, q6y = isometric_projection(sl2, 0, size1, angles)
    q7x, q7y = isometric_projection(sl2, sl3, sl1, angles)
    q8x, q8y = isometric_projection(sl2, sl3, size1, angles)

    # face 13 coordinates
    face13x = [p2x, p1x, p3x, p4x]
    face13y = [p2y, p1y, p3y, p4y]

    # face 23 coordinates
    face23x = [p2x, p8x, p6x, p4x]
    face23y = [p2y, p8y, p6y, p4y]

    # face 12 coordinates
    face12x = [p4x, p3x, p5x, p6x]
    face12y = [p4y, p3y, p5y, p6y]

    # slice 13 coordinates
    slice13x = [q8x, q7x, q5x, q6x]
    slice13y = [q8y, q7y, q5y, q6y]

    # slice 23 coordinates
    slice23x = [q1x, q7x, q5x, q3x]
    slice23y = [q1y, q7y, q5y, q3y]

    # slice 12 coordinates
    slice12x = [q2x, q1x, q7x, q8x]
    slice12y = [q2y, q1y, q7y, q8y]

    # axis limits
    xlim = [p1x - 0.01, p5x + 0.01]
    ylim = [p3y - 0.01, p8y + 0.01]

    # frames
    if args.frametop:
        topline = [[(p2x, p2y), (p8x, p8y), (p6x, p6y)]]
        topline = LineCollection(topline, color='k')
        topline2 = [[(q2x, q2y), (q8x, q8y), (q6x, q6y)]]
        topline2 = LineCollection(topline2, color='k')
    if args.framebottom:
        bottomline = [[(p1x, p1y), (p3x, p3y), (p5x, p5y)]]
        bottomline = LineCollection(bottomline, color='k')
        bottomline2 = [[(q1x, q1y), (q3x, q3y), (q5x, q5y)]]
        bottomline2 = LineCollection(bottomline2, color='k')
    if args.frameleft:
        leftline = [[(p1x, p1y), (p2x, p2y)]]
        leftline = LineCollection(leftline, color='k')
        leftline2 = [[(q1x, q1y), (q2x, q2y)]]
        leftline2 = LineCollection(leftline2, color='k')
    if args.frameright:
        rightline = [[(p5x, p5y), (p6x, p6y)]]
        rightline = LineCollection(rightline, color='k')
        rightline2 = [[(q5x, q5y), (q6x, q6y)]]
        rightline2 = LineCollection(rightline2, color='k')
    if args.centerframe:
        if takeslice:
            centerline = [[(p2x, p2y), (q2x, q2y)], [(p3x, p3y), (q3x, q3y)], [(p6x, p6y), (q6x, q6y)]]
            centerline = LineCollection(centerline, color='k', zorder=1)
        else:
            centerline = [[(p4x, p4y), (p2x, p2y)], [(p4x, p4y), (p3x, p3y)], [(p4x, p4y), (p6x, p6y)]]
            centerline = LineCollection(centerline, color='k', zorder=1)
        centerline2 = [[(q7x, q7y), (q1x, q1y)], [(q7x, q7y), (q5x, q5y)], [(q7x, q7y), (q8x, q8y)]]
        centerline2 = LineCollection(centerline2, color='k', zorder=2, linestyle='dashed')

if octant_2:

    #
    #    A sketch of cavalier projection. 
    #    The relative point positions also apply to isometric projection
    #
    #
    #         8 *------------------------* 6
    #          /.      face 23          /|
    #         / .                      / |
    #        /  .        #............#  |
    #       /   .       .:            :  |
    #      /    .      . :            :  |
    #   2 /     *(7) .  :      (4)    :  * 5 (size3)
    #    *-----------#   #............# /
    #    | size1     : .            .  /
    #    |           :.            .  /   face 13
    #    |           #............#  /
    #    |                        | /
    #    |   face 12              |/ ) angle
    #    *------------------------*  ----
    #    1                        3 (size2)
    #

    # vertices positions for faces
    p1x, p1y = isometric_projection(0, size2, 0, angles)
    p2x, p2y = isometric_projection(0, size2, size1, angles)
    p3x, p3y = isometric_projection(0, 0, 0, angles)
    p4x, p4y = isometric_projection(0, 0, size1, angles)
    p5x, p5y = isometric_projection(size3, 0, 0, angles)
    p6x, p6y = isometric_projection(size3, 0, size1, angles)
    p7x, p7y = isometric_projection(size3, size2, 0, angles)
    p8x, p8y = isometric_projection(size3, size2, size1, angles)

    # set slice cut plane
    if x1end != x1beg:
        scale1 = size1 / (x1end - x1beg)
    else:
        scale1 = 1.0e2
    if x2end != x2beg:
        scale2 = size2 / (x2end - x2beg)
    else:
        scale2 = 1.0e2
    if x3end != x3beg:
        scale3 = size3 / (x3end - x3beg)
    else:
        scale3 = 1.0e2
    sl1 = size1 - (sl1 - x1beg) / ((n1end - n1beg - 1) * d1) * size1
    sl2 = (sl2 - x2beg) / ((n2end - n2beg - 1) * d2) * size2
    sl3 = (sl3 - x3beg) / ((n3end - n3beg - 1) * d3) * size3

    # vertices positions for slices
    q1x, q1y = isometric_projection(0, sl2, sl1, angles)
    q2x, q2y = isometric_projection(0, sl2, size1, angles)
    q3x, q3y = isometric_projection(0, 0, sl1, angles)
    q4x, q4y = isometric_projection(0, 0, size1, angles)
    q5x, q5y = isometric_projection(sl3, 0, sl1, angles)
    q6x, q6y = isometric_projection(sl3, 0, size1, angles)
    q7x, q7y = isometric_projection(sl3, sl2, sl1, angles)
    q8x, q8y = isometric_projection(sl3, sl2, size1, angles)

    # face 13 coordinates
    face13x = [p4x, p3x, p5x, p6x]
    face13y = [p4y, p3y, p5y, p6y]

    # face 23 coordinates
    face23x = [p4x, p2x, p8x, p6x]
    face23y = [p4y, p2y, p8y, p6y]

    # face 12 coordinates
    face12x = [p4x, p3x, p1x, p2x]
    face12y = [p4y, p3y, p1y, p2y]

    # slice 13 coordinates
    slice13x = [q2x, q1x, q7x, q8x]
    slice13y = [q2y, q1y, q7y, q8y]

    # slice 23 coordinates
    slice23x = [q3x, q1x, q7x, q5x]
    slice23y = [q3y, q1y, q7y, q5y]

    # slice 12 coordinates
    slice12x = [q6x, q5x, q7x, q8x]
    slice12y = [q6y, q5y, q7y, q8y]

    # axis limits
    xlim = [p1x - 0.01, p5x + 0.01]
    ylim = [p3y - 0.01, p8y + 0.01]

    # frames
    if args.frametop:
        topline = [[(p2x, p2y), (p8x, p8y), (p6x, p6y)]]
        topline = LineCollection(topline, color='k')
        topline2 = [[(q2x, q2y), (q8x, q8y), (q6x, q6y)]]
        topline2 = LineCollection(topline2, color='k')
    if args.framebottom:
        bottomline = [[(p1x, p1y), (p3x, p3y), (p5x, p5y)]]
        bottomline = LineCollection(bottomline, color='k')
        bottomline2 = [[(q1x, q1y), (q3x, q3y), (q5x, q5y)]]
        bottomline2 = LineCollection(bottomline2, color='k')
    if args.frameleft:
        leftline = [[(p2x, p2y), (p1x, p1y)]]
        leftline = LineCollection(leftline, color='k')
        leftline2 = [[(q2x, q2y), (q1x, q1y)]]
        leftline2 = LineCollection(leftline2, color='k')
    if args.frameright:
        rightline = [[(p6x, p6y), (p5x, p5y)]]
        rightline = LineCollection(rightline, color='k')
        rightline2 = [[(q6x, q6y), (q5x, q5y)]]
        rightline2 = LineCollection(rightline2, color='k')
    if args.centerframe:
        if takeslice:
            centerline = [[(p2x, p2y), (q2x, q2y)], [(p3x, p3y), (q3x, q3y)], [(p6x, p6y), (q6x, q6y)]]
            centerline = LineCollection(centerline, color='k', zorder=1)
        else:
            centerline = [[(p4x, p4y), (p2x, p2y)], [(p4x, p4y), (p6x, p6y)], [(p4x, p4y), (p3x, p3y)]]
            centerline = LineCollection(centerline, color='k', zorder=1)
        centerline2 = [[(q7x, q7y), (q1x, q1y)], [(q7x, q7y), (q5x, q5y)], [(q7x, q7y), (q8x, q8y)]]
        centerline2 = LineCollection(centerline2, color='k', zorder=2, linestyle='dashed')

if octant_3:

    #
    #    A sketch of cavalier projection. 
    #    The relative point positions also apply to isometric projection. 
    #
    #
    #         8 *------------------------* 6
    #          /.                       /|
    #         / .                      / |
    #        /  .        #............#  |
    #       /   .       .:            :  |
    #      /    .      . :            :  |
    #   2 /     *(7) .  :      (4)    :  * 5 (size2)
    #    *-----------#   #............# /
    #    | size1     : .            .  /
    #    |           :.            .  /
    #    |           #............#  /
    #    |                        | /
    #    |                        |/ ) angle
    #    *------------------------*  ----
    #    1                        3 (size3)
    #

    # vertices positions for faces
    p1x, p1y = isometric_projection(0, size3, 0, angles)
    p2x, p2y = isometric_projection(0, size3, size1, angles)
    p3x, p3y = isometric_projection(0, 0, 0, angles)
    p4x, p4y = isometric_projection(0, 0, size1, angles)
    p5x, p5y = isometric_projection(size2, 0, 0, angles)
    p6x, p6y = isometric_projection(size2, 0, size1, angles)
    p7x, p7y = isometric_projection(size2, size3, 0, angles)
    p8x, p8y = isometric_projection(size2, size3, size1, angles)

    # set slice cut plane
    if x1end != x1beg:
        scale1 = size1 / (x1end - x1beg)
    else:
        scale1 = 1.0e2
    if x2end != x2beg:
        scale2 = size2 / (x2end - x2beg)
    else:
        scale2 = 1.0e2
    if x3end != x3beg:
        scale3 = size3 / (x3end - x3beg)
    else:
        scale3 = 1.0e2
    sl1 = size1 - (sl1 - x1beg) / ((n1end - n1beg - 1) * d1) * size1
    sl2 = size2 - (sl2 - x2beg) / ((n2end - n2beg - 1) * d2) * size2
    sl3 = (sl3 - x3beg) / ((n3end - n3beg - 1) * d3) * size3

    # vertices positions for slices
    q1x, q1y = isometric_projection(0, sl3, sl1, angles)
    q2x, q2y = isometric_projection(0, sl3, size1, angles)
    q3x, q3y = isometric_projection(0, 0, sl1, angles)
    q4x, q4y = isometric_projection(0, 0, size1, angles)
    q5x, q5y = isometric_projection(sl2, 0, sl1, angles)
    q6x, q6y = isometric_projection(sl2, 0, size1, angles)
    q7x, q7y = isometric_projection(sl2, sl3, sl1, angles)
    q8x, q8y = isometric_projection(sl2, sl3, size1, angles)

    # face 13 coordinates
    face13x = [p4x, p3x, p1x, p2x]
    face13y = [p4y, p3y, p1y, p2y]

    # face 23 coordinates
    face23x = [p6x, p4x, p2x, p8x]
    face23y = [p6y, p4y, p2y, p8y]

    # face 12 coordinates
    face12x = [p6x, p5x, p3x, p4x]
    face12y = [p6y, p5y, p3y, p4y]

    # slice 13 coordinates
    slice13x = [q6x, q5x, q7x, q8x]
    slice13y = [q6y, q5y, q7y, q8y]

    # slice 23 coordinates
    slice23x = [q5x, q3x, q1x, q7x]
    slice23y = [q5y, q3y, q1y, q7y]

    # slice 12 coordinates
    slice12x = [q8x, q7x, q1x, q2x]
    slice12y = [q8y, q7y, q1y, q2y]

    # axis limits
    xlim = [p1x - 0.01, p5x + 0.01]
    ylim = [p3y - 0.01, p8y + 0.01]

    # frames
    if args.frametop:
        topline = [[(p2x, p2y), (p8x, p8y), (p6x, p6y)]]
        topline = LineCollection(topline, color='k')
        topline2 = [[(q2x, q2y), (q8x, q8y), (q6x, q6y)]]
        topline2 = LineCollection(topline2, color='k')
    if args.framebottom:
        bottomline = [[(p1x, p1y), (p3x, p3y), (p5x, p5y)]]
        bottomline = LineCollection(bottomline, color='k')
        bottomline2 = [[(q1x, q1y), (q3x, q3y), (q5x, q5y)]]
        bottomline2 = LineCollection(bottomline2, color='k')
    if args.frameleft:
        leftline = [[(p1x, p1y), (p2x, p2y)]]
        leftline = LineCollection(leftline, color='k')
        leftline2 = [[(q1x, q1y), (q2x, q2y)]]
        leftline2 = LineCollection(leftline2, color='k')
    if args.frameright:
        rightline = [[(p5x, p5y), (p6x, p6y)]]
        rightline = LineCollection(rightline, color='k')
        rightline2 = [[(q5x, q5y), (q6x, q6y)]]
        rightline2 = LineCollection(rightline2, color='k')
    if args.centerframe:
        if takeslice:
            centerline = [[(p2x, p2y), (q2x, q2y)], [(p3x, p3y), (q3x, q3y)], [(p6x, p6y), (q6x, q6y)]]
            centerline = LineCollection(centerline, color='k', zorder=1)
        else:
            centerline = [[(p4x, p4y), (p2x, p2y)], [(p4x, p4y), (p6x, p6y)], [(p4x, p4y), (p3x, p3y)]]
            centerline = LineCollection(centerline, color='k', zorder=1)
        centerline2 = [[(q7x, q7y), (q1x, q1y)], [(q7x, q7y), (q5x, q5y)], [(q7x, q7y), (q8x, q8y)]]
        centerline2 = LineCollection(centerline2, color='k', zorder=2, linestyle='dashed')

if octant_4:

    #
    #    A sketch of cavalier projection. 
    #    The relative point positions also apply to isometric projection.
    #
    #
    #         8 *------------------------* 6
    #          /.      face 23          /|
    #         / .                      / |
    #        /  .        #............#  |
    #       /   .       .:            :  |
    #      /    .      . :            :  |
    #   2 /     *(7) .  :      (4)    :  * 5 (size3)
    #    *-----------#   #............# /
    #    | size1     : .            .  /
    #    |           :.            .  /   face 13
    #    |           #............#  /
    #    |                        | /
    #    |   face 12              |/ ) angle
    #    *------------------------*  ----
    #    1                        3 (size2)
    #

    # vertices positions for faces
    p1x, p1y = isometric_projection(0, size2, 0, angles)
    p2x, p2y = isometric_projection(0, size2, size1, angles)
    p3x, p3y = isometric_projection(0, 0, 0, angles)
    p4x, p4y = isometric_projection(0, 0, size1, angles)
    p5x, p5y = isometric_projection(size3, 0, 0, angles)
    p6x, p6y = isometric_projection(size3, 0, size1, angles)
    p7x, p7y = isometric_projection(size3, size2, 0, angles)
    p8x, p8y = isometric_projection(size3, size2, size1, angles)

    # set slice cut plane
    if x1end != x1beg:
        scale1 = size1 / (x1end - x1beg)
    else:
        scale1 = 1.0e2
    if x2end != x2beg:
        scale2 = size2 / (x2end - x2beg)
    else:
        scale2 = 1.0e2
    if x3end != x3beg:
        scale3 = size3 / (x3end - x3beg)
    else:
        scale3 = 1.0e2
    sl1 = size1 - (sl1 - x1beg) / ((n1end - n1beg - 1) * d1) * size1
    sl2 = size2 - (sl2 - x2beg) / ((n2end - n2beg - 1) * d2) * size2
    sl3 = size3 - (sl3 - x3beg) / ((n3end - n3beg - 1) * d3) * size3

    # vertices positions for slices
    q1x, q1y = isometric_projection(0, sl2, sl1, angles)
    q2x, q2y = isometric_projection(0, sl2, size1, angles)
    q3x, q3y = isometric_projection(0, 0, sl1, angles)
    q4x, q4y = isometric_projection(0, 0, size1, angles)
    q5x, q5y = isometric_projection(sl3, 0, sl1, angles)
    q6x, q6y = isometric_projection(sl3, 0, size1, angles)
    q7x, q7y = isometric_projection(sl3, sl2, sl1, angles)
    q8x, q8y = isometric_projection(sl3, sl2, size1, angles)

    # face 13 coordinates
    face13x = [p6x, p5x, p3x, p4x]
    face13y = [p6y, p5y, p3y, p4y]

    # face 23 coordinates
    face23x = [p8x, p6x, p4x, p2x]
    face23y = [p8y, p6y, p4y, p2y]

    # face 12 coordinates
    face12x = [p2x, p1x, p3x, p4x]
    face12y = [p2y, p1y, p3y, p4y]

    # slice 13 coordinates
    slice13x = [q8x, q7x, q1x, q2x]
    slice13y = [q8y, q7y, q1y, q2y]

    # slice 23 coordinates
    slice23x = [q7x, q5x, q3x, q1x]
    slice23y = [q7y, q5y, q3y, q1y]

    # slice 12 coordinates
    slice12x = [q8x, q7x, q5x, q6x]
    slice12y = [q8y, q7y, q5y, q6y]

    # axis limits
    xlim = [p1x - 0.01, p5x + 0.01]
    ylim = [p3y - 0.01, p8y + 0.01]

    # frames
    if args.frametop:
        topline = [[(p2x, p2y), (p8x, p8y), (p6x, p6y)]]
        topline = LineCollection(topline, color='k')
        topline2 = [[(q2x, q2y), (q8x, q8y), (q6x, q6y)]]
        topline2 = LineCollection(topline2, color='k')
    if args.framebottom:
        bottomline = [[(p1x, p1y), (p3x, p3y), (p5x, p5y)]]
        bottomline = LineCollection(bottomline, color='k')
        bottomline2 = [[(q1x, q1y), (q3x, q3y), (q5x, q5y)]]
        bottomline2 = LineCollection(bottomline2, color='k')
    if args.frameleft:
        leftline = [[(p2x, p2y), (p1x, p1y)]]
        leftline = LineCollection(leftline, color='k')
        leftline2 = [[(q2x, q2y), (q1x, q1y)]]
        leftline2 = LineCollection(leftline2, color='k')
    if args.frameright:
        rightline = [[(p6x, p6y), (p5x, p5y)]]
        rightline = LineCollection(rightline, color='k')
        rightline2 = [[(q6x, q6y), (q5x, q5y)]]
        rightline2 = LineCollection(rightline2, color='k')
    if args.centerframe:
        if takeslice:
            centerline = [[(p2x, p2y), (q2x, q2y)], [(p3x, p3y), (q3x, q3y)], [(p6x, p6y), (q6x, q6y)]]
            centerline = LineCollection(centerline, color='k', zorder=1)
        else:
            centerline = [[(p4x, p4y), (p2x, p2y)], [(p4x, p4y), (p6x, p6y)], [(p4x, p4y), (p3x, p3y)]]
            centerline = LineCollection(centerline, color='k', zorder=1)
        centerline2 = [[(q7x, q7y), (q1x, q1y)], [(q7x, q7y), (q5x, q5y)], [(q7x, q7y), (q8x, q8y)]]
        centerline2 = LineCollection(centerline2, color='k', zorder=2, linestyle='dashed')

## set colormap
from module_colormap import set_colormap
colormap = set_colormap(args)

## set clip
from module_clip import *
data = np.concatenate((fdata12.flatten(), fdata13.flatten(), fdata23.flatten()))
if args.slice1 is not None:
    data = np.concatenate((data, data23.flatten()))
if args.slice2 is not None:
    data = np.concatenate((data, data13.flatten()))
if args.slice3 is not None:
    data = np.concatenate((data, data12.flatten()))
cmin, cmax = set_clip(args, data, 'fore', dmin, dmax)
if args.norm == 'log':
    if cmin > np.floor(cmax) or cmax < np.ceil(cmin):
        print('error: values in dataset have same order of magnitude')
        exit()

## plot the projection faces and slices
angle1 = angle1

figheight = p8y - p3y
figwidth = max(p5x, p6x) - min(p1x, p2x)

fig = plt.figure(figsize=(figwidth, figheight))
ax = fig.add_axes([0, 0, 1, 1])

project_contour(args, ax, fdata13, face13x, face13y, colormap, cmin, cmax, font, mask13)
project_contour(args, ax, fdata23, face23x, face23y, colormap, cmin, cmax, font, mask23)
project_contour(args, ax, fdata12, face12x, face12y, colormap, cmin, cmax, font, mask12)

if args.frametop:
    ax.add_collection(topline)
if args.framebottom:
    ax.add_collection(bottomline)
if args.frameleft:
    ax.add_collection(leftline)
if args.frameright:
    ax.add_collection(rightline)
if args.centerframe:
    ax.add_collection(centerline)

if args.slice1 is not None or args.slice2 is not None or args.slice3 is not None:

    project_contour(args, ax, data13, slice13x, slice13y, colormap, cmin, cmax, font)
    project_contour(args, ax, data23, slice23x, slice23y, colormap, cmin, cmax, font)
    project_contour(args, ax, data12, slice12x, slice12y, colormap, cmin, cmax, font)

    if args.frametop:
        ax.add_collection(topline2)
    if args.framebottom:
        ax.add_collection(bottomline2)
    if args.frameleft:
        ax.add_collection(leftline2)
    if args.frameright:
        ax.add_collection(rightline2)
    if args.centerframe:
        ax.add_collection(centerline2)

## set ticks
# major ticks style
tick_major_length = float(args.tickmajorlen)
tick_major_width = float(args.tickmajorwid)

# minor ticks style
if args.tickminorlen is None:
    tick_minor_length = 0.5 * tick_major_length
else:
    tick_minor_length = float(args.tickminorlen)

if args.tickminorwid is None:
    tick_minor_width = 0.75 * tick_major_width
else:
    tick_minor_width = float(args.tickminorwid)

# tick font size
label_1_size = float(args.label1size)
label_2_size = float(args.label2size)
label_3_size = float(args.label3size)

if args.tick1size is None:
    tick_1_size = label_1_size - 2
else:
    tick_1_size = float(args.tick1size)

if args.tick2size is None:
    tick_2_size = label_2_size - 2
else:
    tick_2_size = float(args.tick2size)

if args.tick3size is None:
    tick_3_size = label_3_size - 2
else:
    tick_3_size = float(args.tick3size)

if octant_1:

    # axis 1
    if args.axis1loc == 'left' or args.axis1loc == 'both':
        project_axis(ax, p2x, p2y, p1x, p1y, args.ticks1, args.tick1beg, args.tick1end, args.tick1d,
                     args.mtick1, x1beg, x1end, n1, d1, font, tick_major_length, tick_major_width,
                     tick_minor_length, tick_minor_width, size1, 'clock', 'counterclock', tick_1_size,
                     args.label1, 'negative', label_1_size, float(args.label1pad), args.tick1format)

    if args.axis1loc == 'right' or args.axis1loc == 'both':
        project_axis(ax, p6x, p6y, p5x, p5y, args.ticks1, args.tick1beg, args.tick1end, args.tick1d,
                     args.mtick1, x1beg, x1end, n1, d1, font, tick_major_length, tick_major_width,
                     tick_minor_length, tick_minor_width, size1, 'counterclock', 'counterclock', tick_1_size,
                     args.label1, 'positive', label_1_size, float(args.label1pad), args.tick1format)

    # axis 2
    if args.axis2loc == 'top' or args.axis2loc == 'both':
        project_axis(ax, p2x, p2y, p8x, p8y, args.ticks2, args.tick2beg, args.tick2end, args.tick2d,
                     args.mtick2, x2beg, x2end, n2, d2, font, tick_major_length, tick_major_width,
                     tick_minor_length, tick_minor_width, size2, 'counterclock', 'positive', tick_2_size,
                     args.label2, 'positive', label_2_size, float(args.label2pad), args.tick2format)

    if args.axis2loc == 'bottom' or args.axis2loc == 'both':
        project_axis(ax, p3x, p3y, p5x, p5y, args.ticks2, args.tick2beg, args.tick2end, args.tick2d,
                     args.mtick2, x2beg, x2end, n2, d2, font, tick_major_length, tick_major_width,
                     tick_minor_length, tick_minor_width, size2, 'clock', 'positive', tick_2_size,
                     args.label2, 'positive', label_2_size, float(args.label2pad), args.tick2format)

    # axis 3
    if args.axis3loc == 'top' or args.axis3loc == 'both':
        project_axis(ax, p8x, p8y, p6x, p6y, args.ticks3, args.tick3beg, args.tick3end, args.tick3d,
                     args.mtick3, x3beg, x3end, n3, d3, font, tick_major_length, tick_major_width,
                     tick_minor_length, tick_minor_width, size3, 'counterclock', 'positive', tick_3_size,
                     args.label3, 'positive', label_3_size, float(args.label3pad), args.tick3format)

    if args.axis3loc == 'bottom' or args.axis3loc == 'both':
        project_axis(ax, p1x, p1y, p3x, p3y, args.ticks3, args.tick3beg, args.tick3end, args.tick3d,
                     args.mtick3, x3beg, x3end, n3, d3, font, tick_major_length, tick_major_width,
                     tick_minor_length, tick_minor_width, size3, 'clock', 'positive', tick_3_size,
                     args.label3, 'positive', label_3_size, float(args.label3pad), args.tick3format)

if octant_2:

    # axis 1
    if args.axis1loc == 'left' or args.axis1loc == 'both':
        project_axis(ax, p2x, p2y, p1x, p1y, args.ticks1, args.tick1beg, args.tick1end, args.tick1d,
                     args.mtick1, x1beg, x1end, n1, d1, font, tick_major_length, tick_major_width,
                     tick_minor_length, tick_minor_width, size1, 'clock', 'counterclock', tick_1_size,
                     args.label1, 'negative', label_1_size, float(args.label1pad), args.tick1format)

    if args.axis1loc == 'right' or args.axis1loc == 'both':
        project_axis(ax, p6x, p6y, p5x, p5y, args.ticks1, args.tick1beg, args.tick1end, args.tick1d,
                     args.mtick1, x1beg, x1end, n1, d1, font, tick_major_length, tick_major_width,
                     tick_minor_length, tick_minor_width, size1, 'counterclock', 'counterclock', tick_1_size,
                     args.label1, 'positive', label_1_size, float(args.label1pad), args.tick1format)

    # axis 2
    if args.axis2loc == 'top' or args.axis2loc == 'both':
        project_axis(ax, p6x, p6y, p8x, p8y, args.ticks2, args.tick2beg, args.tick2end, args.tick2d,
                     args.mtick2, x2beg, x2end, n2, d2, font, tick_major_length, tick_major_width,
                     tick_minor_length, tick_minor_width, size2, 'clock', 'negative', tick_2_size,
                     args.label2, 'negative', label_2_size, float(args.label2pad), args.tick2format)

    if args.axis2loc == 'bottom' or args.axis2loc == 'both':
        project_axis(ax, p3x, p3y, p1x, p1y, args.ticks2, args.tick2beg, args.tick2end, args.tick2d,
                     args.mtick2, x2beg, x2end, n2, d2, font, tick_major_length, tick_major_width,
                     tick_minor_length, tick_minor_width, size2, 'counterclock', 'negative', tick_2_size,
                     args.label2, 'negative', label_2_size, float(args.label2pad), args.tick2format)

    # axis 3
    if args.axis3loc == 'top' or args.axis3loc == 'both':
        project_axis(ax, p2x, p2y, p8x, p8y, args.ticks3, args.tick3beg, args.tick3end, args.tick3d,
                     args.mtick3, x3beg, x3end, n3, d3, font, tick_major_length, tick_major_width,
                     tick_minor_length, tick_minor_width, size3, 'counterclock', 'positive', tick_3_size,
                     args.label3, 'positive', label_3_size, float(args.label3pad), args.tick3format)

    if args.axis3loc == 'bottom' or args.axis3loc == 'both':
        project_axis(ax, p3x, p3y, p5x, p5y, args.ticks3, args.tick3beg, args.tick3end, args.tick3d,
                     args.mtick3, x3beg, x3end, n3, d3, font, tick_major_length, tick_major_width,
                     tick_minor_length, tick_minor_width, size3, 'clock', 'positive', tick_3_size,
                     args.label3, 'positive', label_3_size, float(args.label3pad), args.tick3format)

if octant_3:

    # axis 1
    if args.axis1loc == 'left' or args.axis1loc == 'both':
        project_axis(ax, p2x, p2y, p1x, p1y, args.ticks1, args.tick1beg, args.tick1end, args.tick1d,
                     args.mtick1, x1beg, x1end, n1, d1, font, tick_major_length, tick_major_width,
                     tick_minor_length, tick_minor_width, size1, 'clock', 'counterclock', tick_1_size,
                     args.label1, 'negative', label_1_size, float(args.label1pad), args.tick1format)

    if args.axis1loc == 'right' or args.axis1loc == 'both':
        project_axis(ax, p6x, p6y, p5x, p5y, args.ticks1, args.tick1beg, args.tick1end, args.tick1d,
                     args.mtick1, x1beg, x1end, n1, d1, font, tick_major_length, tick_major_width,
                     tick_minor_length, tick_minor_width, size1, 'counterclock', 'counterclock', tick_1_size,
                     args.label1, 'positive', label_1_size, float(args.label1pad), args.tick1format)

    # axis 2
    if args.axis2loc == 'top' or args.axis2loc == 'both':
        project_axis(ax, p8x, p8y, p2x, p2y, args.ticks2, args.tick2beg, args.tick2end, args.tick2d,
                     args.mtick2, x2beg, x2end, n2, d2, font, tick_major_length, tick_major_width,
                     tick_minor_length, tick_minor_width, size2, 'clock', 'negative', tick_2_size,
                     args.label2, 'negative', label_2_size, float(args.label2pad), args.tick2format)

    if args.axis2loc == 'bottom' or args.axis2loc == 'both':
        project_axis(ax, p5x, p5y, p3x, p3y, args.ticks2, args.tick2beg, args.tick2end, args.tick2d,
                     args.mtick2, x2beg, x2end, n2, d2, font, tick_major_length, tick_major_width,
                     tick_minor_length, tick_minor_width, size2, 'counterclock', 'negative', tick_2_size,
                     args.label2, 'negative', label_2_size, float(args.label2pad), args.tick2format)

    # axis 3
    if args.axis3loc == 'top' or args.axis3loc == 'both':
        project_axis(ax, p6x, p6y, p8x, p8y, args.ticks3, args.tick3beg, args.tick3end, args.tick3d,
                     args.mtick3, x3beg, x3end, n3, d3, font, tick_major_length, tick_major_width,
                     tick_minor_length, tick_minor_width, size3, 'clock', 'negative', tick_3_size,
                     args.label3, 'negative', label_3_size, float(args.label3pad), args.tick3format)

    if args.axis3loc == 'bottom' or args.axis3loc == 'both':
        project_axis(ax, p3x, p3y, p1x, p1y, args.ticks3, args.tick3beg, args.tick3end, args.tick3d,
                     args.mtick3, x3beg, x3end, n3, d3, font, tick_major_length, tick_major_width,
                     tick_minor_length, tick_minor_width, size3, 'counterclock', 'negative', tick_3_size,
                     args.label3, 'negative', label_3_size, float(args.label3pad), args.tick3format)

if octant_4:

    # axis 1
    if args.axis1loc == 'left' or args.axis1loc == 'both':
        project_axis(ax, p2x, p2y, p1x, p1y, args.ticks1, args.tick1beg, args.tick1end, args.tick1d,
                     args.mtick1, x1beg, x1end, n1, d1, font, tick_major_length, tick_major_width,
                     tick_minor_length, tick_minor_width, size1, 'clock', 'counterclock', tick_1_size,
                     args.label1, 'negative', label_1_size, float(args.label1pad), args.tick1format)

    if args.axis1loc == 'right' or args.axis1loc == 'both':
        project_axis(ax, p6x, p6y, p5x, p5y, args.ticks1, args.tick1beg, args.tick1end, args.tick1d,
                     args.mtick1, x1beg, x1end, n1, d1, font, tick_major_length, tick_major_width,
                     tick_minor_length, tick_minor_width, size1, 'counterclock', 'counterclock', tick_1_size,
                     args.label1, 'positive', label_1_size, float(args.label1pad), args.tick1format)

    # axis 2
    if args.axis2loc == 'top' or args.axis2loc == 'both':
        project_axis(ax, p8x, p8y, p6x, p6y, args.ticks2, args.tick2beg, args.tick2end, args.tick2d,
                     args.mtick2, x2beg, x2end, n2, d2, font, tick_major_length, tick_major_width,
                     tick_minor_length, tick_minor_width, size2, 'counterclock', 'positive', tick_2_size,
                     args.label2, 'positive', label_2_size, float(args.label2pad), args.tick2format)

    if args.axis2loc == 'bottom' or args.axis2loc == 'both':
        project_axis(ax, p1x, p1y, p3x, p3y, args.ticks2, args.tick2beg, args.tick2end, args.tick2d,
                     args.mtick2, x2beg, x2end, n2, d2, font, tick_major_length, tick_major_width,
                     tick_minor_length, tick_minor_width, size2, 'clock', 'positive', tick_2_size,
                     args.label2, 'positive', label_2_size, float(args.label2pad), args.tick2format)

    # axis 3
    if args.axis3loc == 'top' or args.axis3loc == 'both':
        project_axis(ax, p8x, p8y, p2x, p2y, args.ticks3, args.tick3beg, args.tick3end, args.tick3d,
                     args.mtick3, x3beg, x3end, n3, d3, font, tick_major_length, tick_major_width,
                     tick_minor_length, tick_minor_width, size3, 'clock', 'negative', tick_3_size,
                     args.label3, 'negative', label_3_size, float(args.label3pad), args.tick3format)

    if args.axis3loc == 'bottom' or args.axis3loc == 'both':
        project_axis(ax, p5x, p5y, p3x, p3y, args.ticks3, args.tick3beg, args.tick3end, args.tick3d,
                     args.mtick3, x3beg, x3end, n3, d3, font, tick_major_length, tick_major_width,
                     tick_minor_length, tick_minor_width, size3, 'counterclock', 'negative', tick_3_size,
                     args.label3, 'negative', label_3_size, float(args.label3pad), args.tick3format)

# remove original figure frames
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)

# remove original figure ticks and labels
ax.tick_params(which='both',
               top='off',
               bottom='off',
               labeltop='off',
               labelbottom='off',
               left='off',
               right='off',
               labelleft='off',
               labelright='off')

## set colorbar
# when using filled contour, the colorbar coding could be a little tricky
# therefore dismissed -- little use anyway
# colorbar for filled contour is also dismissed in 2D case (showcontour)
# if args.legend and cmin != cmax:
if False:

    lloc = args.lloc

    # colorbar values
    cp = 512
    cinterval = (cmax - cmin) / (cp - 1)
    temp = np.linspace(cmin, cmax, cp)
    cval = np.zeros([cp, 1])
    for i in range(0, cp):
        cval[i] = temp[i]

    if lloc in ['left', 'right']:

        if args.lheight is None:
            lheight = figheight
        else:
            lheight = float(args.lheight)
            if lheight > figheight:
                lheight = figheight

        if args.lwidth is None:
            lwidth = 0.2
        else:
            lwidth = float(args.lwidth)

        if lloc == 'right':

            # colorbar location
            if args.axis1loc == 'right' or args.axis1loc == 'both':
                pad = tick_1_size * ipp + float(args.tickmajorlen) * ipp + 0.5 + label_1_size * ipp
            else:
                pad = 0.2
            if args.lpad is None:
                cbpad = 0.0
            else:
                cbpad = float(args.lpad)
            pad = pad + cbpad
            cx = [p5x + pad, p5x + pad, p5x + pad + lwidth, p5x + pad + lwidth]
            dl = (figheight - lheight) / 2.0
            cy = [p3y + dl, p8y - dl, p8y - dl, p3y + dl]

    if args.lloc in ['top', 'bottom']:

        if args.lheight is None:
            lheight = 0.2
        else:
            lheight = float(args.lheight)

        if args.lwidth is None:
            lwidth = figwidth
        else:
            lwidth = float(args.lwidth)
            if lwidth > figwidth:
                lwidth = figwidth

        if lloc == 'bottom':

            # colorbar location
            if args.axis2loc == 'bottom' or args.axis2loc == 'both' or args.axis3loc == 'bottom' or args.axis3loc == 'both':
                pad = max(tick_2_size, tick_3_size) * ipp + float(args.tickmajorlen) * ipp + 0.2 + max(
                    label_2_size, label_3_size) * ipp + 0.1 * abs(cos(angle1 * np.pi / 180.0))
            else:
                pad = 0.2
            if args.lpad is None:
                cbpad = 0.0
            else:
                cbpad = float(args.lpad)
            pad = pad + cbpad
            dl = (figwidth - lwidth) / 2.0
            cx = [p1x + dl, p5x - dl, p5x - dl, p1x + dl]
            cy = [p3y - pad, p3y - pad, p3y - pad - lheight, p3y - pad - lheight]

    # create colorbar
    cb = project_image(ax, cval, cx, cy, colormap, cmin, cmax)

    # crate colorbar frame
    line = [[(cx[0], cy[0]), (cx[1], cy[1]), (cx[2], cy[2]), (cx[3], cy[3]), (cx[0], cy[0])]]
    line = LineCollection(line, linewidth=1.0, color='k')
    ax.add_collection(line)

    if args.unitsize is None:
        lufs = min(float(args.label1size), float(args.label2size), float(args.label3size)) - 1
    else:
        lufs = float(args.unitsize)

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

        # set ticks
        ticks = np.arange(ltickbeg, ltickend + ld, ld)
        tbeg = max(cmin, ltickbeg)
        tend = min(cmax, ltickend)

        # set tick positions on colorbar
        ticks = np.asarray(
            [i for i in ticks if i >= tbeg - 1.0e-10 * abs(tbeg) and i <= tend + 1.0e-10 * abs(tend)])
        tick_labels = ['' for i in range(0, len(ticks))]
        for i in range(0, len(ticks)):
            tick_labels[i] = ('%f' % (ticks[i] / cscale)).rstrip('0').rstrip('.')

        # set minor ticks
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
        tbeg = max(cmin, ltickbeg)
        tend = min(cmax, ltickend)

        # set tick positions on colorbar
        ticks = np.asarray(
            [i for i in ticks if i >= tbeg - 1.0e-10 * abs(tbeg) and i <= tend + 1.0e-10 * abs(tend)])
        tick_labels = ['' for i in range(0, len(ticks))]
        for i in range(0, len(ticks)):
            tick_labels[i] = '$\mathregular{10^{%i}}$' % (ticks[i])

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

    # add ticks by drawing lines
    if lloc == 'right':

        # add major ticks
        cbeg = cy[3]
        cend = cy[2]
        ticks = [(i - cmin) / (cmax - cmin) * (cend - cbeg) + cbeg for i in ticks]

        tx = [cx[2] for i in range(0, len(ticks))]
        ty = [i for i in ticks]
        last_tick = ty[-1]
        ticklen = float(args.tickmajorlen) * ipp
        tx2 = [i + ticklen for i in tx]
        ty2 = [i for i in ticks]
        ttx = [i + ticklen + 0.05 for i in tx]
        tty = [i for i in ticks]

        majortick = []
        e1 = list(zip(tx, ty))
        e2 = list(zip(tx2, ty2))
        for i in range(0, len(ticks)):
            majortick.append([e1[i], e2[i]])

        tick_major_width = 0.5
        majortick = LineCollection(majortick, linewidths=tick_major_width, colors='k')
        ax.add_collection(majortick)

        # add tick labels
        if args.norm == 'linear' and cscale != 1.0:
            tick_labels[-1] = r'$\mathregular{\times 10^{%i}}$' % scalar + '\n' + tick_labels[-1]
        for i in range(0, len(ticks)):
            ax.text(ttx[i], tty[i], tick_labels[i], fontproperties=font, size=ltfs, ha='left', va='center')

        # add minor ticks
        if args.lmtick != 0:
            mticks = [(i - cmin) / (cmax - cmin) * (cend - cbeg) + cbeg for i in mticks]

            tx = [cx[2] for i in range(0, len(mticks))]
            ty = [i for i in mticks]
            ticklen = float(args.tickmajorlen) * ipp * 0.5
            tx2 = [i + ticklen for i in tx]
            ty2 = [i for i in mticks]
            ttx = [i + ticklen + 0.05 for i in tx]
            tty = [i for i in mticks]

            minortick = []
            e1 = list(zip(tx, ty))
            e2 = list(zip(tx2, ty2))
            for i in range(0, len(mticks)):
                minortick.append([e1[i], e2[i]])

            tick_minor_width = 0.5 * 0.75
            minortick = LineCollection(minortick, linewidths=tick_minor_width, colors='k')
            ax.add_collection(minortick)

        # set unit
        if args.unit is not None:
            if args.unitpad is None:
                upad = 0.05
            else:
                upad = float(args.unitpad)
            if args.norm == 'linear':
                maxlen = max([len(i) for i in tick_labels]) * 0.75
            if args.norm == 'log':
                maxlen = 3.5
            ux = cx[2] + ticklen * 2.0 + 0.025 + maxlen * ltfs * ipp + upad
            uy = 0.5 * (cbeg + cend)
            ct = ax.text(ux,
                         uy,
                         args.unit,
                         size=lufs,
                         fontproperties=font,
                         rotation=270,
                         ha='left',
                         va='center')
            ct.set_size(lufs)

    # add ticks by drawing lines
    if lloc == 'bottom':

        # add major ticks
        cbeg = cx[0]
        cend = cx[1]
        ticks = [(i - cmin) / (cmax - cmin) * (cend - cbeg) + cbeg for i in ticks]

        ty = [cy[3] for i in range(0, len(ticks))]
        tx = [i for i in ticks]
        last_tick = tx[-1]
        ticklen = float(args.tickmajorlen) * ipp
        ty2 = [i - ticklen for i in ty]
        tx2 = [i for i in ticks]
        tty = [i - ticklen - 0.05 for i in ty]
        ttx = [i for i in ticks]

        majortick = []
        e1 = list(zip(tx, ty))
        e2 = list(zip(tx2, ty2))
        for i in range(0, len(ticks)):
            majortick.append([e1[i], e2[i]])

        tick_major_width = 0.5
        majortick = LineCollection(majortick, linewidths=tick_major_width, colors='k')
        ax.add_collection(majortick)

        # add tick labels
        if args.norm == 'linear' and cscale != 1.0:
            tick_labels[-1] = tick_labels[-1] + '\n' + r'$\mathregular{\times 10^{%i}}$' % scalar
        for i in range(0, len(ticks)):
            ct = ax.text(ttx[i],
                         tty[i],
                         tick_labels[i],
                         fontproperties=font,
                         size=ltfs,
                         ha='center',
                         va='top')
            ct.set_size(ltfs)

        # add minor ticks
        if args.lmtick != 0:
            mticks = [(i - cmin) / (cmax - cmin) * (cend - cbeg) + cbeg for i in mticks]

            ty = [cy[3] for i in range(0, len(mticks))]
            tx = [i for i in mticks]
            ticklen = float(args.tickmajorlen) * ipp * 0.5
            ty2 = [i - ticklen for i in ty]
            tx2 = [i for i in mticks]
            tty = [i - ticklen - 0.05 for i in ty]
            ttx = [i for i in mticks]

            minortick = []
            e1 = list(zip(tx, ty))
            e2 = list(zip(tx2, ty2))
            for i in range(0, len(mticks)):
                minortick.append([e1[i], e2[i]])

            tick_minor_width = 0.5 * 0.75
            minortick = LineCollection(minortick, linewidths=tick_minor_width, colors='k')
            ax.add_collection(minortick)

        # set unit
        if args.unit is not None:
            if args.unitpad is None:
                upad = 0.05
            else:
                upad = float(args.unitpad)
            if args.norm == 'linear':
                maxlen = 1.50
            if args.norm == 'log':
                maxlen = 1.75
            ux = 0.5 * (cbeg + cend)
            uy = cy[3] - ticklen - 0.025 - maxlen * ltfs * ipp - upad
            ct = ax.text(ux, uy, args.unit, size=lufs, fontproperties=font, ha='center', va='top')
            ct.set_size(lufs)

## set title
if args.title is not None:

    if args.titlesize is None:
        title_font_size = max(float(args.label1size), float(args.label2size), float(args.label3size)) + 2
    else:
        title_font_size = float(args.titlesize)

    if args.titlex is None:
        if octant_1:
            title_x = 0.5 * (p2x + p6x)
        if octant_2:
            title_x = 0.5 * (p8x + p4x)
    else:
        title_x = float(args.titlex)

    if args.titley is None:
        if args.axis2loc == 'top' or args.axis2loc == 'both' or args.axis3loc == 'top' or args.axis3loc == 'both':
            title_y = p8y + 2 * max(label_2_size, label_3_size) * ipp + 0.1
        else:
            title_y = p8y + 0.25
    else:
        title_y = float(args.titley)

    ax.text(title_x,
            title_y,
            args.title,
            ha='center',
            fontproperties=fontbold,
            fontweight='bold',
            size=title_font_size)

## set axis limits and aspect
extra0 = 0.0
extra1 = 0.0

# x axis low limit
if args.axis1loc == 'left' or args.axis1loc == 'both':
    extra0 = extra0 + float(args.tickmajorlen) * ipp * 1.25
    xlim0 = xlim[0] - extra0
else:
    xlim0 = xlim[0] - 0.05

# x axis top limit
if args.axis1loc == 'right' or args.axis1loc == 'both':
    extra1 = extra1 + float(args.tickmajorlen) * ipp * 1.25
    xlim1 = xlim[1] + extra1
else:
    xlim1 = xlim[1] + 0.1

ax.set_xlim([xlim0, xlim1])

extra0 = 0.0
extra1 = 0.0

# y axis low limit
if args.axis3loc == 'top' or args.axis3loc == 'both':
    extra0 = extra0 + float(args.tickmajorlen) * ipp * 1.25
    ylim0 = ylim[0] - extra0
else:
    ylim0 = ylim[0] - 0.1

# y axis top limit
if args.axis3loc == 'bottom' or args.axis3loc == 'both':
    extra1 = extra1 + float(args.tickmajorlen) * ipp * 1.25
    ylim1 = ylim[1] + extra1
else:
    ylim1 = ylim[1] + 0.05

ax.set_ylim([ylim0, ylim1])

# set axis to appropritate ratio
ax.set_aspect('equal')

plt.tick_params(
    axis='x',  # changes apply to the x1-axis
    which='both',  # both major and minor ticks are affected
    bottom=0,  # ticks along the bottom axis
    top=0,  # ticks along the top axis
    labelbottom=0,  # labels along the bottom axis
    labeltop=0)  # labels along the top axis
plt.tick_params(
    axis='y',  # changes apply to the x2-axis
    which='both',  # both major and minor ticks are affected
    left=0,  # ticks along the left axis
    right=0,  # ticks along the right axis
    labelleft=0,  # labels along the left axis
    labelright=0)  # labels along the right axis

## output
output(args)

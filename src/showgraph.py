#!/usr/bin/python3

from module_io import *
from math import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from module_annotation import set_default
from matplotlib.patches import *
from module_title import *
from module_tick import rigid_tick_label
from module_utility import *
from module_frame import *
from module_font import *
from module_range import *
from itertools import cycle
from pylab import *
import matplotlib as mplt
import numpy as np
import matplotlib.pyplot as plt
import os
import math
import argparse
from module_getarg import *
from argparse import RawTextHelpFormatter
# from sklearn.metrics import r2_score

# ignore warnings
import warnings
warnings.filterwarnings("ignore", module="matplotlib")

# tag
program = 'graph'
print()

# read arguments
parser = argparse.ArgumentParser(description='''
                                purpose:
                                    Plot 1D array(s) as curves or scatters.
                                ''',
                                 formatter_class=RawTextHelpFormatter)
parser = getarg(parser, program)
args = parser.parse_args()

# input data
# get names of input files
infile = args.infile[0].split(',')
ninput = len(infile)

print()
# check if input files exist
for i in range(0, ninput):
    if not os.path.exists(infile[i]):
        print('input file', infile[i], 'does not exists')
        print()
        exit()
    else:
        print('input <<    ', infile[i])

# concatanate files if necessary
if ninput > 1:

    if args.filetype == 'ascii':
        fftype = 'w'
        rftype = 'r'
    if args.filetype == 'binary':
        fftype = 'wb'
        rftype = 'rb'

    import string
    import random
    tempfile = "".join([random.choice(string.ascii_letters + string.digits) for n in range(10)])
    tempfile = '/tmp/tempall' + tempfile
    with open(tempfile, fftype) as tempall:
        for file in infile:
            with open(file, rftype) as tempsingle:
                tempall.write(tempsingle.read())
    infile = tempfile

else:

    infile = infile[0]

# plot type
if args.plottype == 1:
    ns = 1
if args.plottype == 2:
    ns = 2
if args.plottype == 3:
    ns = 3

# fpr ascii file, data can be directly read
if args.filetype == 'ascii':

    # from textual file
    data = np.loadtxt(infile, ndmin=2)
    if args.transpose:
        data = data.transpose()

    # shape of input data
    shape = data.shape
    print('shape:      ', shape)

    # get first dimension
    if args.n1 is None:
        n = np.array([shape[0]])
        nps = n
        nf = 1
    else:
        n = args.n1[0].split(',')
        n = [int(i) for i in n]
        nps = sum(n)
        if nps > shape[0]:
            print('error: first dimension specification exceeds bounds of data')
            exit()
        nf = len(n)

    # second dimension size is definite when ascii file given
    n2 = shape[1]

# for binary file, set dimensions first
if args.filetype == 'binary':

    # get number of data points in file
    n0 = os.path.getsize(infile)
    if n0 == 0:
        print('error: input file empty')
        exit()
    datatype = args.datatype
    if datatype == 'double':
        n0 = int(np.floor(n0 * 1.0 / 8))
    if datatype == 'float':
        n0 = int(np.floor(n0 * 1.0 / 4))
    if datatype == 'int':
        n0 = int(np.floor(n0 * 1.0 / 2))

    if args.n1 is None:
        nf = 1
        # if first dimension not given, all data goes to 1st dimension
        n = np.array([n0])
        # second dimension is calculated or given
        if args.n2 is None:
            n2 = ns
        else:
            n2 = int(args.n2)
        # resize first dimension
        n = np.array([int(np.floor(n * 1.0 / n2))])
        nps = n[0]
    else:
        n = args.n1[0].split(',')
        n = [int(i) for i in n]
        nps = sum(n)
        nf = len(n)
        # second dimension is calculated or given
        if args.n2 is None:
            n2 = int(np.floor(n0 * 1.0 / nps))
        else:
            n2 = int(args.n2)
            if n2 > int(np.floor(n0 * 1.0 / nps)):
                print('error: second dimension size exceeds bounds of data')
                exit()
        nmax = max(n)

# check second dimension size
if n2 < ns:
    print(n2, ns)
    print('error: second dimension size smaller than that required by plot type')
    exit()

# set index range for each set of data
nrange = np.zeros([nf, 2], dtype=int)
for i in range(0, nf - 1):
    nrange[i, 1] = int(n[i] + nrange[i, 0])
    nrange[i + 1, 0] = int(nrange[i, 1])
nrange[nf - 1, 1] = int(nps)

# read binary file here
if args.filetype == 'binary':

    # data type
    from module_datatype import *
    dt = set_datatype(args)

    # read binary file
    data = fromfile(infile, dtype=dt, count=nps * n2)
    if not args.transpose:
        data = data.reshape((n2, nps))
        data = data.transpose()
    else:
        data = data.reshape((nps, n2))

    print('shape:      ', data.shape)

# select part of the 2nd dimension
if args.select is not None:
    if (args.select)[0] == 'all':
        select = np.array([i for i in range(0, n2)])
    else:
        select = np.array([int(i) - 1 for i in args.select[0].split(',')])
    if select.min() < 0 or select.max() > n2 - 1 or len(select) % ns != 0:
        print('selection of 2nd dimension error')
        exit()
    data = data[:, select]
    nset = int(len(select) / ns)
    data0 = data[:, 0:ns]
    nf = nf * nset
    n0 = n[:]
    nrange0 = nrange[:, :]
    for i in range(1, nset):
        data0 = np.append(data0, data[:, ns * i:ns * (i + 1)], axis=0)
        n = np.append(n, n0, axis=0)
        nrange = np.append(nrange, nrange0 + nps * i, axis=0)
    nps = nps * nset
    data = data0[:, :]
    
if args.normalize2 == 'max':
    for i in range(size(data, axis=1)):
        data[:, i] = data[:, i]/max(data[:, i])

# plot
if args.projection == 'cartesian':
    fig = plt.figure()  # figsize=(figwidth,figheight))
    ax = fig.add_axes([0, 0, 1, 1])
if args.projection == 'polar':
    fig = plt.figure()
    plt.axes(projection='polar')
    ax = plt.gca()

# get plot styles
# line style
if args.linestyle is not None:
    linestyle = args.linestyle[0].split(',')
    if len(linestyle) < nf:
        l = len(linestyle)
        alinestyle = ['none' for i in range(l, nf)]
        linestyle.extend(alinestyle)
else:
    linestyle = ['none' for i in range(0, nf)]

# line width
if args.linewidth is not None:
    linewidth = args.linewidth[0].split(',')
    if len(linewidth) < nf:
        l = len(linewidth)
        alinewidth = [1.0 for i in range(l, nf)]
        linewidth.extend(alinewidth)
else:
    linewidth = [1.0 for i in range(0, nf)]
linewidth = [float(i) for i in linewidth]

# line color
defaultcolor = cycle(['blue', 'red', 'green', 'yellow', 'black', 'cyan', 'magenta'])
if args.linecolor is not None:
    linecolor = args.linecolor[0].split(',')
    nc = len(linecolor)
else:
    linecolor = []
    nc = 0
if nc < nf:
    ic = 0
    for i in cycle(defaultcolor):
        linecolor.append(i)
        ic = ic + 1
        if ic >= nf:
            break

# line transparency
if args.linealpha is not None:
    linealpha = args.linealpha[0].split(',')
    if len(linealpha) < nf:
        l = len(linealpha)
        alinealpha = [1.0 for i in range(l, nf)]
        linealpha.extend(alinealpha)
else:
    linealpha = [1.0 for i in range(0, nf)]
linealpha = [float(i) for i in linealpha]

# line close or not
if args.close is None:
    close = []
else:
    close = args.close
if len(close) > 0:
    for i in range(nf - len(close)):
        close.extend(False)
else:
    close =  [False for i in range(nf)]

# marker style
if args.plottype == 1 or args.plottype == 2:
    defaultmarker = 'none'
if args.plottype == 3:
    defaultmarker = 'o'
if args.marker is not None:
    marker = args.marker[0].split(',')
    if len(marker) < nf:
        l = len(marker)
        amarker = [defaultmarker for i in range(l, nf)]
        marker.extend(amarker)
else:
    marker = [defaultmarker for i in range(0, nf)]

# marker size
if args.markersize is not None:
    markersize = args.markersize[0].split(',')
    if len(markersize) < nf:
        l = len(markersize)
        amarkersize = [2.0 for i in range(l, nf)]
        markersize.extend(amarkersize)
else:
    markersize = [2.0 for i in range(0, nf)]
markersize = [float(i) for i in markersize]

# marker every
# mark every can be set to be more complicated according to
# matplotlib docs, see http://matplotlib.org/api/lines_api.html
#
# every=None, every point will be plotted.
# every=N, every N-th marker will be plotted starting with marker 0.
# every=(start, N), every N-th marker, starting at point start, will be plotted.
# every=slice(start, end, N), every N-th marker, starting at point start, upto
# but not including point end, will be plotted.
# every=[i, j, m, n], only markers at points i, j, m, and n will be plotted.
# every=0.1, (i.e. a float) then markers will be spaced at approximately equal
# distances along the line; the distance along the line between markers is
# determined by multiplying the display-coordinate distance of the axes
# bounding-box diagonal by the value of every.
# every=(0.5, 0.1) (i.e. a length-2 tuple of float), the same functionality
# as every=0.1 is exhibited but the first marker will be 0.5 multiplied by
# the display-cordinate-diagonal-distance along the line.
#
# here simple integer every is used
if args.markerevery is not None:
    markerevery = args.markerevery[0].split(',')
    if len(markerevery) < nf:
        l = len(markerevery)
        amarkerevery = [1 for i in range(l, nf)]
        markerevery.extend(amarkerevery)
else:
    markerevery = [1 for i in range(0, nf)]
markerevery = [int(i) for i in markerevery]

# marker face color
if args.markerfacecolor is not None:
    markerfacecolor = args.markerfacecolor[0].split(',')
    if len(markerfacecolor) < nf:
        l = len(markerfacecolor)
        amarkerfacecolor = [linecolor[i] for i in range(l, nf)]
        markerfacecolor.extend(amarkerfacecolor)
else:
    markerfacecolor = [linecolor[i] for i in range(0, nf)]

# marker edge color
if args.plottype == 1 or args.plottype == 2:
    if args.markeredgecolor is not None:
        markeredgecolor = args.markeredgecolor[0].split(',')
        if len(markeredgecolor) < nf:
            l = len(markeredgecolor)
            amarkeredgecolor = [markerfacecolor[i] for i in range(l, nf)]
            markeredgecolor.extend(amarkeredgecolor)
    else:
        markeredgecolor = [markerfacecolor[i] for i in range(0, nf)]

for i in range(0, nf):
    if markerfacecolor[i] == 'none' and markeredgecolor[i] == 'none':
        markeredgecolor[i] = linecolor[i]

if args.plottype == 3:
    if args.markeredgecolor is not None:
        markeredgecolor = args.markeredgecolor[0].split(',')
        if len(markeredgecolor) < nf:
            l = len(markeredgecolor)
            amarkeredgecolor = ['k' for i in range(l, nf)]
            markeredgecolor.extend(amarkeredgecolor)
    else:
        markeredgecolor = ['k' for i in range(0, nf)]

# marker transparency
if args.markeralpha is not None:
    markeralpha = args.markeralpha[0].split(',')
    if len(markeralpha) < nf:
        l = len(markeralpha)
        amarkeralpha = [1.0 for i in range(l, nf)]
        markeralpha.extend(amarkeralpha)
else:
    markeralpha = [1.0 for i in range(0, nf)]
markeralpha = [float(i) for i in markeralpha]

# marker size
if args.markersizemax is not None:
    markersizemax = args.markersizemax[0].split(',')
    if len(markersizemax) < nf:
        l = len(markersizemax)
        amarkersizemax = [15.0 for i in range(l, nf)]
        markersizemax.extend(amarkersizemax)
else:
    markersizemax = [15.0 for i in range(0, nf)]
markersizemax = [float(i) for i in markersizemax]

if args.markersizemin is not None:
    markersizemin = args.markersizemin[0].split(',')
    if len(markersizemin) < nf:
        l = len(markersizemin)
        amarkersizemin = [5.0 for i in range(l, nf)]
        markersizemin.extend(amarkersizemin)
else:
    markersizemin = [5.0 for i in range(0, nf)]
markersizemin = [float(i) for i in markersizemin]

# plot labels
if args.plotlabel is not None:
    plotlabel = args.plotlabel[0].split(':')
    if len(plotlabel) < nf:
        l = len(plotlabel)
        aplotlabel = ['Set ' + str(i) for i in range(l, nf)]
        plotlabel.extend(aplotlabel)
else:
    plotlabel = ['Set ' + str(i) for i in range(0, nf)]

locdict = {
    'upper_right': 1,
    'upper_left': 2,
    'lower_left': 3,
    'lower_right': 4,
    'right': 5,
    'center_left': 6,
    'center_right': 7,
    'lower_center': 8,
    'upper_center': 9,
    'center': 10
}
if args.plotlabelloc in list(locdict.keys()):
    labelloc = locdict[args.plotlabelloc]
else:
    labelloc = 2

# if both line and marker are not specified, then plot line
for i in range(0, nf):
    if linestyle[i] == 'none' and marker[i] == 'none':
        linestyle[i] = '-'

# if 1 or 2 set fill or not
if args.plottype == 1 or args.plottype == 2:

    # fill w.r.t
    if args.fillwrt is not None:
        fillwrt = args.fillwrt[0].split(',')
        if len(fillwrt) < nf:
            l = len(fillwrt)
            afillwrt = ['none' for i in range(l, nf)]
            fillwrt.extend(afillwrt)
    else:
        fillwrt = ['none' for i in range(0, nf)]

    # fill color positive
    if args.fillcolor is not None:
        fillcolor = args.fillcolor[0].split(',')
        if len(fillcolor) < nf:
            l = len(fillcolor)
            afillcolor = [linecolor[i] for i in range(l, nf)]
            fillcolor.extend(afillcolor)
    else:
        fillcolor = [linecolor[i] for i in range(0, nf)]

    # fill color transparency
    if args.fillalpha is not None:
        fillalpha = args.fillalpha[0].split(',')
        if len(fillalpha) < nf:
            l = len(fillalpha)
            afillalpha = [0.5 for i in range(l, nf)]
            fillalpha.extend(afillalpha)
    else:
        fillalpha = [0.5 for i in range(0, nf)]
    fillalpha = [float(i) for i in fillalpha]

# set colormap
if ns == 3 or ns == 4:
    from module_colormap import set_colormap
    colormap = set_colormap(args)

# set plot range
if args.plottype == 1:

    # set axis 1
    d1 = float(args.d1)

    # sampling point begin and end
    sp1beg = float(args.o1)
    sp1end = sp1beg + (max(n) - 1) * d1

    # limit of axis
    if args.x1beg is None:
        x1beg = sp1beg
    else:
        x1beg = float(args.x1beg)
    if args.x1end is None:
        x1end = sp1end
    else:
        x1end = float(args.x1end)

    if args.norm1 == 'log':
        if x1beg <= 0 and x1end <= 0:
            print('error: inappropriate axis 1 limits for log plot')
            exit()
        if x1beg <= 0 and x1end > 0:
            print('warning: axis 1 starts with inappropriate value for log plot, setting to %e' % abs(d1))
            x1beg = abs(d1)
        if x1beg > 0 and x1end <= 0:
            print('warning: axis 1 ends with inappropriate value for log plot, setting to %e' % abs(d1))
            x1end = abs(d1)

    # set axis 2
    ymin = data[~isnan(data)].min()
    ymax = data[~isnan(data)].max()
    extray = 0.05 * abs(ymax - ymin)
    if args.x2beg is None:
        x2beg = data.min() - extray
    else:
        x2beg = float(args.x2beg)
    if args.x2end is None:
        x2end = data.max() + extray
    else:
        x2end = float(args.x2end)

    if args.norm2 == 'log':
        if x2beg <= 0 and x2end <= 0:
            print('error: inappropriate axis 2 limits for log plot')
            exit()
        if x2beg <= 0 and x2end > 0:
            print('warning: axis 2 starts with inappropriate value for log plot, setting to %e' %
                  (1.0e-5 * x2end))
            x2beg = 1.0e-5 * x2end
        if x2beg > 0 and x2end <= 0:
            print('warning: axis 2 ends with inappropriate value for log plot, setting to %e' %
                  (1.0e5 * x2beg))
            x2end = 1.0e5 * x2beg

    print('value range ', ymin, ' -- ', ymax)
    print('plot range  ', x2beg, ' -- ', x2end)

if args.plottype == 2 or args.plottype == 3:

    # set axis 1
    x = data[:, 0]
    xmin = x[~isnan(x)].min()
    xmax = x[~isnan(x)].max()
    extrax = 0.05 * abs(xmax - xmin)
    if args.x1beg is None:
        x1beg = xmin - extrax
    else:
        x1beg = float(args.x1beg)
    if args.x1end is None:
        x1end = xmax + extrax
    else:
        x1end = float(args.x1end)

    if args.norm1 == 'log':
        if x1beg <= 0 and x1end <= 0:
            print('error: inappropriate axis 1 limits for log plot')
            exit()
        if x1beg <= 0 and x1end > 0:
            print('warning: axis 1 starts with inappropriate value for log plot, setting to %e' %
                  (1.0e-5 * x1end))
            x1beg = 1.0e-5 * x1end
        if x1beg > 0 and x1end <= 0:
            print('warning: axis 1 ends with inappropriate value for log plot, setting to %e' %
                  (1.0e5 * x1beg))
            x1end = 1.0e5 * x1beg

    # set axis 2
    y = data[:, 1]
    ymin = y[~isnan(y)].min()
    ymax = y[~isnan(y)].max()
    extray = 0.05 * abs(ymax - ymin)
    if args.x2beg is None:
        x2beg = ymin - extray
    else:
        x2beg = float(args.x2beg)
    if args.x2end is None:
        x2end = ymax + extray
    else:
        x2end = float(args.x2end)

    if args.norm2 == 'log':
        if x2beg <= 0 and x2end <= 0:
            print('error: inappropriate axis 2 limits for log plot')
            exit()
        if x2beg <= 0 and x2end > 0:
            print('warning: axis 2 starts with inappropriate value for log plot, setting to %e' %
                  (1.0e-5 * x2end))
            x2beg = 1.0e-5 * x2end
        if x2beg > 0 and x2end <= 0:
            print('warning: axis 2 ends with inappropriate value for log plot, setting to %e' %
                  (1.0e5 * x2beg))
            x2end = 1.0e5 * x2beg

    if args.plottype == 2:
        print('value range ', xmin, ' -- ', xmax)
        print('value range ', ymin, ' -- ', ymax)
        print('plot range  ', x1beg, ' -- ', x1end)
        print('plot range  ', x2beg, ' -- ', x2end)

    # color clip
    if args.plottype == 3:
        # scale data to log scale if necessary
        if args.norm == 'log':
            data[:, 2] = np.log10(data[:, 2])

        # marker size
        temp = data[:, 2]
        temp = temp[~isnan(temp)]
        zmin = temp.min()
        zmax = temp.max()
        msize = []
        for i in range(0, nf):
            if zmin == zmax:
                msize = np.append(msize, [markersizemin[i] for j in data[nrange[i, 0]:nrange[i, 1], 2]])
            else:
                msize = np.append(msize, [(j - zmin) / (zmax - zmin) *
                                          (markersizemax[i] - markersizemin[i]) + markersizemin[i]
                                          for j in data[nrange[i, 0]:nrange[i, 1], 2]])

        print('value range ', zmin, ' -- ', zmax)

        # set clip
        from module_clip import *
        cmin, cmax = set_clip(args, data[:, 2], 'fore', zmin, zmax)
        if args.norm == 'log':
            if cmin > np.floor(cmax) or cmax < np.ceil(cmin):
                print('error: values in dataset have same order of magnitude')
                exit()

# set axis limits
if args.direction == 'horizontal':
    if not args.reverse1:
        ax.set_xlim([x1beg, x1end])
    else:
        ax.set_xlim([x1end, x1beg])
    if not args.reverse2:
        ax.set_ylim([x2beg, x2end])
    else:
        ax.set_ylim([x2end, x2beg])
if args.direction == 'vertical':
    if not args.reverse1:
        ax.set_ylim([x1beg, x1end])
    else:
        ax.set_ylim([x1end, x1beg])
    if not args.reverse2:
        ax.set_xlim([x2beg, x2end])
    else:
        ax.set_xlim([x2end, x2beg])

# set figure size
figbase = 5.0
golden_ratio = 1.61803398875

if args.size1 is None:
    size1 = figbase
else:
    size1 = float(args.size1)

range1 = abs(x1end - x1beg)
range2 = abs(x2end - x2beg)
ratio = range1 / range2

if args.size2 is None:
    if ratio >= 1.0 / 3.0 and ratio <= 3.0:
        size2 = size1 / ratio
    else:
        size2 = figbase / golden_ratio
else:
    size2 = float(args.size2)

figheight = size2
figwidth = size1

if args.direction == 'vertical':
    figheight = size1
    figwidth = size2

# reset figure size
if args.projection == 'cartesian':
    fig.set_size_inches(figwidth, figheight)

# set font
font, fontbold = set_font(args)

# plot
# iterate through all sets
for i in range(0, nf):

    if args.plottype == 1 or args.plottype == 2:
        # line plot given y or (x,y)

        if ns == 1:
            if args.direction == 'horizontal':
                # define x
                x = np.linspace(sp1beg, sp1end, n[i])
                # define y
                y = data[nrange[i, 0]:nrange[i, 1]]
            if args.direction == 'vertical':
                # define y
                y = np.linspace(sp1beg, sp1end, n[i])
                # define x
                x = data[nrange[i, 0]:nrange[i, 1]]

        if ns == 2:
            if args.direction == 'horizontal':
                # define x
                x = data[nrange[i, 0]:nrange[i, 1], 0]
                # define y
                y = data[nrange[i, 0]:nrange[i, 1], 1]
            if args.direction == 'vertical':
                # define y
                x = data[nrange[i, 0]:nrange[i, 1], 1]
                # define x
                y = data[nrange[i, 0]:nrange[i, 1], 0]

            if close[i]:
                x = np.append(x, x[0])
                y = np.append(y, y[0])

        # line styles
        linestyles = [
            'solid',
            '-',  # solid line
            'dashed',
            '--',  # dashed line
            'dashdot',
            '-.',  # dashdot line
            'dotted',
            ':'  # dotted line
        ]

        # marker styles
        markers = [
            '.',  # point marker
            ',',  # pixel marker
            'o',  # circle marker
            'v',  # triangle_down marker
            '^',  # triangle_up marker
            '<',  # triangle_left marker
            '>',  # triangle_right marker
            '1',  # tri_down marker
            '2',  # tri_up marker
            '3',  # tri_left marker
            '4',  # tri_right marker
            's',  # square marker
            'p',  # pentagon marker
            '*',  # star marker
            'h',  # hexagon1 marker
            'H',  # hexagon2 marker
            '+',  # plus marker
            'x',  # x marker
            'D',  # diamond marker
            'd',  # thin_diamond marker
            '|',  # vline marker
            '_'  # hline marker
        ]

        # plot or scatter depending on options
        if linestyle[i] in linestyles and marker[i] == 'none':
            ax.plot(x,
                    y,
                    linestyle=linestyle[i],
                    linewidth=linewidth[i],
                    color=linecolor[i],
                    label=plotlabel[i],
                    markevery=markerevery[i],
                    antialiased=True,
                    zorder=args.plotorder)
        if linestyle[i] == 'none' and marker[i] in markers:
            ax.plot(x,
                    y,
                    linestyle='none',
                    label=plotlabel[i],
                    marker=marker[i],
                    markersize=markersize[i],
                    markerfacecolor=markerfacecolor[i],
                    markeredgecolor=markeredgecolor[i],
                    markevery=markerevery[i],
                    zorder=args.plotorder)
        if linestyle[i] in linestyles and marker[i] in markers:
            ax.plot(x,
                    y,
                    linestyle=linestyle[i],
                    linewidth=linewidth[i],
                    color=linecolor[i],
                    antialiased=True,
                    label=plotlabel[i],
                    marker=marker[i],
                    markersize=markersize[i],
                    markerfacecolor=markerfacecolor[i],
                    markeredgecolor=markeredgecolor[i],
                    markevery=markerevery[i],
                    zorder=args.plotorder)

        # add legend
        if args.plotlabel is not None:
            if args.plotlabelloc in list(locdict.keys()):
                lg = plt.legend(loc=labelloc, ncol=args.plotlabelcol)
            else:
                lg = plt.legend(bbox_to_anchor=(1.01, 1.0), loc=2, borderaxespad=0, ncol=args.plotlabelcol)
            lg.set_zorder(10)
            leg = ax.get_legend()
            ltext = leg.get_texts()
            plt.setp(ltext, fontproperties=font)
            plt.setp(ltext, fontsize=float(args.plotlabelsize))

        # add fill color if necessary
        if fillwrt[i] != 'none':
            if args.direction == 'horizontal':
                ax.fill_between(x[:],
                                float(fillwrt[i]),
                                y[:, 0],
                                facecolor=fillcolor[i],
                                edgecolor='none',
                                alpha=fillalpha[i])
            if args.direction == 'vertical':
                ax.fill_betweenx(y[:],
                                 float(fillwrt[i]),
                                 x[:, 0],
                                 facecolor=fillcolor[i],
                                 edgecolor='none',
                                 alpha=fillalpha[i])

    if args.plottype == 3:
        # scatter plot given (x,y,v)

        if args.direction == 'horizontal':
            # define x
            x = data[nrange[i, 0]:nrange[i, 1], 0]
            # define y
            y = data[nrange[i, 0]:nrange[i, 1], 1]
            # define v
            v = data[nrange[i, 0]:nrange[i, 1], 2]
        if args.direction == 'vertical':
            # define y
            x = data[nrange[i, 0]:nrange[i, 1], 1]
            # define x
            y = data[nrange[i, 0]:nrange[i, 1], 0]
            # define v
            v = data[nrange[i, 0]:nrange[i, 1], 2]

        # scatter plot with different sizes and colors
        sc = ax.scatter(
            x,
            y,
            c=v,  # variable color
            s=msize[nrange[i, 0]:nrange[i, 1]], # variable size
            vmin=cmin,  # clip min
            vmax=cmax,  # clip max
            marker=marker[i],
            zorder=args.plotorder)
        sc.set_cmap(colormap)
        sc.set_antialiaseds(True)
        if args.markeredgecolor is None:
            sc.set_edgecolor('face')
            sc.set_linewidths(0)
        else:
            sc.set_edgecolor(markeredgecolor[i])
            sc.set_linewidths(0.5)
        sc.set_alpha(markeralpha[i])

# set frame
set_frame(args)


# set tick
def define_tick(ticks,
                tickbeg,
                tickend,
                tickd,
                mtick,
                xbeg,
                xend,
                axislen,
                norm='linear',
                base=10,
                format='sci'):

    scalar = []

    if norm == 'linear':

        # regular ticks
        if ticks is None:

            # scale tick values if too large or too small
            maxtick = max(abs(xbeg), abs(xend))
            if (maxtick >= 1.0e4 or maxtick <= 1.0e-3) and format == 'sci':
                scalar = int(floor(log10(maxtick)))
                cscale = pow(10, scalar)
            else:
                cscale = 1.0

            # scale values
            xbeg = xbeg / cscale
            xend = xend / cscale

            # major tick interval
            if tickd is None:
                tick_interval = nice((xend - xbeg) / 5.0)
                if tick_interval == 0:
                    tick_interval = 1.0e10
            else:
                tick_interval = float(tickd) / cscale

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
                tick_beg = float(tickbeg) / cscale

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
                tick_end = float(tickend) / cscale

            # regular major and minor tick locations
            tick = np.arange(tick_beg, tick_end + 0.1 * abs(tick_interval), tick_interval)
            minor_tick_interval = tick_interval / (mtick + 1.0)
            minor_tick = np.arange(tick_beg, tick_end + 0.1 * abs(minor_tick_interval), minor_tick_interval)

            # some ticks might out of axis range
            tick = np.asarray([i for i in tick if i >= xbeg and i <= xend])
            minor_tick = np.asarray([i for i in minor_tick if i >= xbeg and i <= xend and (i not in tick)])

            # set major tick location and labels, note some major ticks might
            # be out of axis range
            if format == 'sci' or format == 'plain':
                tick_label = [('%f' % i).rstrip('0').rstrip('.') for i in tick]
            else:
                tick_label = [(format % i) for i in tick]

            # rescale ticks
            tick = [i * cscale for i in tick]
            minor_tick = [i * cscale for i in minor_tick]

        # irregular ticks
        else:

            # get contents from user-specified ticks
            ticks = ticks[0].split(':')
            location = [0 for i in range(0, len(ticks))]
            label = ['' for i in range(0, len(ticks))]

            # set tick locations
            for i in range(0, len(ticks)):
                t = ticks[i].split(',')
                location[i] = float(t[0])
                label[i] = t[1]

            # sort according to tick location
            yx = sorted(zip(location, label))
            tick = [location for location, label in yx]
            tick_label = [label for location, label in yx]

            # minor ticks
            if mtick != 0:
                mtick = mtick + 1
                minor_tick = np.linspace(tick[0], tick[1], mtick + 1)
                minor_tick = minor_tick[1:mtick]
                for i in range(1, len(tick) - 1):
                    t = np.linspace(tick[i], tick[i + 1], mtick + 1)
                    minor_tick = np.append(minor_tick, t[1:mtick])
            else:
                minor_tick = []

    if norm == 'log':

        # regular ticks
        if ticks is None:

            # major ticks
            if tickbeg is None:
                if xbeg <= 0:
                    print('error: invalid axis start value')
                    exit()
                else:
                    tick_beg = np.floor(math.log(xbeg, base))
            else:
                if float(tickbeg) <= 0:
                    print('error: invalid tick start value')
                    exit()
                else:
                    tick_beg = math.log(float(tickbeg), base)
            if tickend is None:
                if xend <= 0:
                    print('error: invalid axis end value')
                    exit()
                else:
                    tick_end = np.ceil(math.log(xend, base)) + 1
            else:
                if float(tickend) <= 0:
                    print('error: invalid tick end value')
                else:
                    tick_end = math.log(float(tickend), base)
            if tickd is None:
                tick_interval = max(1, round((tick_beg - tick_end) / 5.0))
            else:
                if float(tickd) <= 0:
                    print('error: invalid tick interval')
                    exit()
                else:
                    tick_interval = math.log(float(tickd), base)

            ticks = np.arange(tick_beg, tick_end + 1, tick_interval)
            ticks = [int(round(i)) for i in ticks]
            tbeg = max(math.log(xbeg, base), tick_beg)
            tend = min(math.log(xend, base), tick_end)
            tick = np.asarray(
                [i for i in ticks if i >= tbeg - 1.0e-10 * abs(tbeg) and i <= tend + 1.0e-10 * abs(tend)])

            # tick labels
            tick_label = ['' for i in range(0, len(tick))]
            for i in range(0, len(tick)):
                tick_label[i] = '$' + str(base) + '^{%i}$' % (tick[i])

            # minor ticks
            if mtick != 0:
                # extend tail and head
                pticks = np.append(ticks, ticks[0] - tick_interval)
                pticks = np.append(pticks, ticks[-1] + tick_interval)
                # sort all major ticks
                pticks = np.sort(pticks)
                # get pseudo-location of minor ticks
                nt = len(pticks)
                # print pticks
                mticks = []
                for i in range(0, nt - 1):
                    tlog = np.linspace(np.float_power(base, pticks[i]), np.float_power(base, pticks[i + 1]),
                                       mtick + 2)
                    for j in range(0, len(tlog)):
                        if np.abs(math.log(tlog[j], base) - pticks[i]) > 1.0e-2 and np.abs(
                                math.log(tlog[j], base) - pticks[min(i + 1, nt - 1)]) > 1.0e-2:
                            mticks = np.append(mticks, math.log(tlog[j], base))
                minor_tick = np.asarray([
                    i for i in mticks if i >= tbeg - 1.0e-10 * abs(tbeg) and i <= tend + 1.0e-10 * abs(tend)
                ])
                # repower minor tick
                minor_tick = [np.float_power(base, i) for i in minor_tick]
            else:
                minor_tick = []

            # repower tick: it is necessary to have the conversion np.float32(tick)
            # otherwise overflow
            tick = [np.float_power(base, i) for i in np.float32(tick)]

        # irregular ticks
        else:

            # get contents from user-specified ticks
            ticks = ticks[0].split(':')
            location = [0 for i in range(0, len(ticks))]
            label = ['' for i in range(0, len(ticks))]

            # set tick locations
            for i in range(0, len(ticks)):
                t = ticks[i].split(',')
                location[i] = math.log(float(t[0]), base)
                label[i] = t[1]

            # sort according to tick location
            yx = sorted(zip(location, label))
            tick = [location for location, label in yx]
            tick_label = [label for location, label in yx]

            # minor ticks
            if mtick != 0:
                mtick = mtick + 1
                minor_tick = np.linspace(np.float_power(base, tick[0]), np.float_power(base, tick[1]),
                                         mtick + 1)
                minor_tick = minor_tick[1:mtick]
                for i in range(1, len(tick) - 1):
                    t = np.linspace(np.float_power(base, tick[i]), np.float_power(base, tick[i + 1]),
                                    mtick + 1)
                    for j in range(1, mtick):
                        minor_tick = np.append(minor_tick, math.log(t[j], base))
            else:
                minor_tick = []

            # repower tick
            tick = [np.float_power(base, i) for i in tick]

    # return major tick location, major tick label and minor tick location
    return tick, tick_label, minor_tick, scalar

if args.projection == 'cartesian':

    if args.direction == 'horizontal':
        if args.norm1 == 'log':
            ax.set_xscale('log', nonpositive='clip', base=int(args.base1))
        if args.norm2 == 'log':
            ax.set_yscale('log', nonpositive='clip', base=int(args.base2))
    if args.direction == 'vertical':
        if args.norm1 == 'log':
            ax.set_yscale('log', nonpositive='clip', base=int(args.base1))
        if args.norm2 == 'log':
            ax.set_xscale('log', nonpositive='clip', base=int(args.base2))

    if args.ticktop:
        ticktop = 1
    else:
        ticktop = 0

    if args.tickbottom:
        tickbottom = 1
    else:
        tickbottom = 0

    if args.tickleft:
        tickleft = 1
    else:
        tickleft = 0

    if args.tickright:
        tickright = 1
    else:
        tickright = 0

    if args.direction == 'horizontal':

        # label font size
        label_1_size = float(args.label1size)
        label_2_size = float(args.label2size)

        # set axis label
        xlabel = ax.set_xlabel(args.label1, fontsize=label_1_size, labelpad=float(args.label1pad))
        ylabel = ax.set_ylabel(args.label2, fontsize=label_2_size, labelpad=float(args.label2pad))

        l = ax.yaxis.get_label()
        l.set_fontproperties(font)
        l.set_fontsize(label_1_size)
        l = ax.xaxis.get_label()
        l.set_fontproperties(font)
        l.set_fontsize(label_2_size)

        if args.label1loc is None:
            label1loc = 'bottom'
        else:
            label1loc = args.label1loc
        if args.label2loc is None:
            label2loc = 'left'
        else:
            label2loc = args.label2loc
        ax.xaxis.set_label_position(label1loc)
        ax.yaxis.set_label_position(label2loc)
        if label2loc == 'right':
            ylabel.set_rotation(270)

        # ticks on/off
        ax.get_yaxis().set_tick_params(which='both', direction='out')
        ax.get_xaxis().set_tick_params(which='both', direction='out')
        plt.tick_params(
            axis='x',  # changes apply to the x1-axis
            which='both',  # both major and minor ticks are affected
            bottom=tickbottom,  # ticks along the bottom axis
            top=ticktop,  # ticks along the top axis
            labelbottom=tickbottom,  # labels along the bottom axis
            labeltop=ticktop)  # labels along the top axis
        plt.tick_params(
            axis='y',  # changes apply to the x2-axis
            which='both',  # both major and minor ticks are affected
            left=tickleft,  # ticks along the left axis
            right=tickright,  # ticks along the right axis
            labelleft=tickleft,  # labels along the left axis
            labelright=tickright)  # labels along the right axis

        # if tick font size and family not speciefied, then inherit from axis
        # labels
        if args.tick1size is None:
            tick_1_font_size = label_1_size - 2
        else:
            tick_1_font_size = float(args.tick1size)

        if args.tick2size is None:
            tick_2_font_size = label_2_size - 2
        else:
            tick_2_font_size = float(args.tick2size)

        # axis 1
        tick1, tick1_label, minor_tick1, scalar1 = define_tick(args.ticks1, args.tick1beg, args.tick1end,
                                                               args.tick1d, args.mtick1, x1beg, x1end, size1,
                                                               args.norm1, int(args.base1), args.tick1format)
        ax.xaxis.set_ticks(tick1)
        ax.xaxis.set_ticks(minor_tick1, minor=True)
        if not args.tick1label:
            ax.xaxis.set_ticklabels([])
        else:
            ax.xaxis.set_ticklabels(tick1_label)
            for l in ax.get_xticklabels():
                l.set_rotation(float(args.tick1rot))
                l.set_fontproperties(font)
                l.set_fontsize(tick_1_font_size)

        # axis 2
        tick2, tick2_label, minor_tick2, scalar2 = define_tick(args.ticks2, args.tick2beg, args.tick2end,
                                                               args.tick2d, args.mtick2, x2beg, x2end, size2,
                                                               args.norm2, int(args.base2), args.tick2format)
        ax.yaxis.set_ticks(tick2)
        ax.yaxis.set_ticks(minor_tick2, minor=True)
        if not args.tick2label:
            ax.yaxis.set_ticklabels([])
        else:
            ax.yaxis.set_ticklabels(tick2_label)
            for l in ax.get_yticklabels():
                l.set_rotation(float(args.tick2rot))
                l.set_fontproperties(font)
                l.set_fontsize(tick_2_font_size)

        # major and minor ticks sytle
        ax.tick_params('both', length=float(args.tickmajorlen), width=float(args.tickmajorwid), which='major')

        # minor ticks style
        if args.tickminorlen is None:
            tick_minor_length = 0.5 * float(args.tickmajorlen)
        else:
            tick_minor_length = float(args.tickminorlen)

        if args.tickminorwid is None:
            tick_minor_width = 0.75 * float(args.tickmajorwid)
        else:
            tick_minor_width = float(args.tickminorwid)

        ax.tick_params('both', length=tick_minor_length, width=tick_minor_width, which='minor')

        # add axis scientific power
        ipp = 0.0138889

        if scalar1 != []:
            if not args.reverse1:
                ha = 'left'
                tt = r'$\mathregular{\times 10^{%i}}$' % scalar1
            else:
                ha = 'right'
                tt = r'$\mathregular{10^{%i} \times}$' % scalar1
            if not args.reverse2:
                py1 = x2end + 1.1 * (float(args.tickmajorlen) * ipp +
                                     tick_1_font_size * ipp) / size2 * abs(x2end - x2beg)
                py2 = x2beg - 1.15 * (float(args.tickmajorlen) * ipp +
                                      tick_1_font_size * ipp) / size2 * abs(x2end - x2beg)
            else:
                py2 = x2end + 1.1 * (float(args.tickmajorlen) * ipp +
                                     tick_1_font_size * ipp) / size2 * abs(x2end - x2beg)
                py1 = x2beg - 1.15 * (float(args.tickmajorlen) * ipp +
                                      tick_1_font_size * ipp) / size2 * abs(x2end - x2beg)
            if args.ticktop:
                px = min(tick1[-1], x1end) + 0.1 * abs(x1end - x1beg) / size1
                tx = ax.text(px, py1, tt, ha=ha, va='bottom')
                tx.set_fontproperties(font)
                tx.set_fontsize(tick_1_font_size)
            if args.tickbottom:
                px = min(tick1[-1], x1end) + 0.1 * abs(x1end - x1beg) / size1
                tx = ax.text(px, py2, tt, ha=ha, va='top')
                tx.set_fontproperties(font)
                tx.set_fontsize(tick_1_font_size)

        if scalar2 != []:
            if not args.reverse1:
                ha1 = 'right'
                ha2 = 'left'
                tt1 = r'$\mathregular{10^{%i} \times}$' % scalar2
                tt2 = r'$\mathregular{\times 10^{%i}}$' % scalar2
                px1 = x1beg - 0.1 * abs(x1end - x1beg) / size1
                px2 = x1end + 0.1 * abs(x1end - x1beg) / size1
            else:
                ha2 = 'left'
                ha1 = 'right'
                tt2 = r'$\mathregular{\times 10^{%i}}$' % scalar2
                tt1 = r'$\mathregular{10^{%i} \times}$' % scalar2
                px2 = x1beg - 0.1 * abs(x1end - x1beg) / size1
                px1 = x1end + 0.1 * abs(x1end - x1beg) / size1
            if not args.reverse2:
                # py=max(x2end+0.01*abs(x2end-x2beg)/size2,tick2[-1]+0.075*abs(x2end-x2beg)/size2)
                py = min(x2end + 0.1 * abs(x2end - x2beg) / size2,
                         tick2[-1] + 0.1 * abs(x2end - x2beg) / size2)
                va = 'bottom'
            else:
                # py=min(x2beg-0.01*abs(x2end-x2beg)/size2,tick2[0]-0.075*abs(x2end-x2beg)/size2)
                py = max(x2beg - 0.1 * abs(x2end - x2beg) / size2,
                         tick2[0] - 0.1 * abs(x2end - x2beg) / size2)
                va = 'bottom'
            if args.tickleft:
                tx = ax.text(px1, py, tt1, ha=ha1, va=va)
                tx.set_fontproperties(font)
                tx.set_fontsize(tick_2_font_size)
            if args.tickright:
                tx = ax.text(px2, py, tt2, ha=ha2, va=va)
                tx.set_fontproperties(font)
                tx.set_fontsize(tick_2_font_size)

    if args.direction == 'vertical':

        # label font size
        label_1_size = float(args.label1size)
        label_2_size = float(args.label2size)

        # set axis label
        ylabel = ax.set_ylabel(args.label1, fontsize=label_1_size, labelpad=float(args.label1pad))
        xlabel = ax.set_xlabel(args.label2, fontsize=label_2_size, labelpad=float(args.label2pad))

        l = ax.yaxis.get_label()
        l.set_fontproperties(font)
        l.set_fontsize(label_1_size)
        l = ax.xaxis.get_label()
        l.set_fontproperties(font)
        l.set_fontsize(label_2_size)

        if args.label1loc is None:
            label1loc = 'left'
        else:
            label1loc = args.label1loc
        if args.label2loc is None:
            label2loc = 'bottom'
        else:
            label2loc = args.label2loc
        ax.yaxis.set_label_position(label1loc)
        ax.xaxis.set_label_position(label2loc)
        if label1loc == 'right':
            ylabel.set_rotation(270)

        # ticks on/off
        ax.get_xaxis().set_tick_params(which='both', direction='out')
        ax.get_yaxis().set_tick_params(which='both', direction='out')
        plt.tick_params(
            axis='x',  # changes apply to the x1-axis
            which='both',  # both major and minor ticks are affected
            bottom=tickbottom,  # ticks along the bottom axis
            top=ticktop,  # ticks along the top axis
            labelbottom=tickbottom,  # labels along the bottom axis
            labeltop=ticktop)  # labels along the top axis
        plt.tick_params(
            axis='y',  # changes apply to the x2-axis
            which='both',  # both major and minor ticks are affected
            left=tickleft,  # ticks along the left axis
            right=tickright,  # ticks along the right axis
            labelleft=tickleft,  # labels along the left axis
            labelright=tickright)  # labels along the right axis

        # if tick font size and family not speciefied, then inherit from axis
        # labels
        if args.tick1size is None:
            tick_1_font_size = label_1_size - 2
        else:
            tick_1_font_size = float(args.tick1size)

        if args.tick2size is None:
            tick_2_font_size = label_2_size - 2
        else:
            tick_2_font_size = float(args.tick2size)

        # axis 1
        tick1, tick1_label, minor_tick1, scalar1 = define_tick(args.ticks1, args.tick1beg, args.tick1end,
                                                               args.tick1d, args.mtick1, x1beg, x1end, size1,
                                                               args.norm1, int(args.base1), args.tick1format)
        ax.yaxis.set_ticks(tick1)
        ax.yaxis.set_ticks(minor_tick1, minor=True)
        if not args.tick1label:
            ax.yaxis.set_ticklabels([])
        else:
            ax.yaxis.set_ticklabels(tick1_label)
            for l in ax.get_yticklabels():
                l.set_rotation(float(args.tick1rot))
                l.set_fontproperties(font)
                l.set_fontsize(tick_1_font_size)

        # axis 2
        tick2, tick2_label, minor_tick2, scalar2 = define_tick(args.ticks2, args.tick2beg, args.tick2end,
                                                               args.tick2d, args.mtick2, x2beg, x2end, size2,
                                                               args.norm2, int(args.base2), args.tick2format)
        ax.xaxis.set_ticks(tick2)
        ax.xaxis.set_ticks(minor_tick2, minor=True)
        if not args.tick2label:
            ax.xaxis.set_ticklabels([])
        else:
            ax.xaxis.set_ticklabels(tick2_label)
            for l in ax.get_xticklabels():
                l.set_rotation(float(args.tick2rot))
                l.set_fontproperties(font)
                l.set_fontsize(tick_2_font_size)

        # major and minor ticks sytle
        ax.tick_params('both', length=float(args.tickmajorlen), width=float(args.tickmajorwid), which='major')

        # minor ticks style
        if args.tickminorlen is None:
            tick_minor_length = 0.5 * float(args.tickmajorlen)
        else:
            tick_minor_length = float(args.tickminorlen)

        if args.tickminorwid is None:
            tick_minor_width = 0.75 * float(args.tickmajorwid)
        else:
            tick_minor_width = float(args.tickminorwid)

        ax.tick_params('both', length=tick_minor_length, width=tick_minor_width, which='minor')

        # add axis scientific power
        ipp = 0.0138889

        if scalar2 != []:
            if not args.reverse2:
                ha = 'left'
                tt = r'$\mathregular{\times 10^{%i}}$' % scalar2
            else:
                ha = 'right'
                tt = r'$\mathregular{10^{%i} \times}$' % scalar2
            if not args.reverse1:
                py1 = x1end + 1.1 * (float(args.tickmajorlen) * ipp +
                                     tick_2_font_size * ipp) / size1 * abs(x1end - x1beg)
                py2 = x1beg - 1.15 * (float(args.tickmajorlen) * ipp +
                                      tick_2_font_size * ipp) / size1 * abs(x1end - x1beg)
            else:
                py2 = x1end + 1.1 * (float(args.tickmajorlen) * ipp +
                                     tick_2_font_size * ipp) / size1 * abs(x1end - x1beg)
                py1 = x1beg - 1.15 * (float(args.tickmajorlen) * ipp +
                                      tick_2_font_size * ipp) / size1 * abs(x1end - x1beg)
            if args.ticktop:
                px = min(tick2[-1], x2end) + 0.1 * abs(x2end - x2beg) / size2
                tx = ax.text(px, py1, tt, ha=ha, va='bottom')
                tx.set_fontproperties(font)
                tx.set_fontsize(tick_2_font_size)
            if args.tickbottom:
                px = min(tick2[-1], x2end) + 0.1 * abs(x2end - x2beg) / size2
                tx = ax.text(px, py2, tt, ha=ha, va='top')
                tx.set_fontproperties(font)
                tx.set_fontsize(tick_2_font_size)

        if scalar1 != []:
            if not args.reverse2:
                ha1 = 'right'
                ha2 = 'left'
                tt1 = r'$\mathregular{10^{%i} \times}$' % scalar1
                tt2 = r'$\mathregular{\times 10^{%i}}$' % scalar1
                px1 = x2beg - 0.1 * abs(x2end - x2beg) / size2
                px2 = x2end + 0.1 * abs(x2end - x2beg) / size2
            else:
                ha2 = 'left'
                ha1 = 'right'
                tt2 = r'$\mathregular{\times 10^{%i}}$' % scalar1
                tt1 = r'$\mathregular{10^{%i} \times}$' % scalar1
                px2 = x2beg - 0.1 * abs(x2end - x2beg) / size2
                px1 = x2end + 0.1 * abs(x2end - x2beg) / size2
            if not args.reverse1:
                py = max(x1end, tick1[-1]) + 0.1 * abs(x1end - x1beg) / size1
                va = 'bottom'
            else:
                py = min(x1beg, tick1[0]) - 0.1 * abs(x1end - x1beg) / size1
                va = 'bottom'
            if args.tickleft:
                tx = ax.text(px1, py, tt1, ha=ha1, va=va)
                tx.set_fontproperties(font)
                tx.set_fontsize(tick_1_font_size)
            if args.tickright:
                tx = ax.text(px2, py, tt2, ha=ha2, va=va)
                tx.set_fontproperties(font)
                tx.set_fontsize(tick_1_font_size)

if args.projection == 'polar':

    ax.set_rmin(x2beg)
    ax.set_rmax(x2end)
    ax.set_theta_zero_location("E")

    # label font size
    label_1_size = float(args.label1size)
    label_2_size = float(args.label2size)

    # if tick font size and family not speciefied, then inherit from axis
    # labels
    if args.tick1size is None:
        tick_1_font_size = label_1_size - 2
    else:
        tick_1_font_size = float(args.tick1size)

    if args.tick2size is None:
        tick_2_font_size = label_2_size - 2
    else:
        tick_2_font_size = float(args.tick2size)
        
    # axis 1
    tick1, tick1_label, minor_tick1, scalar1 = define_tick(args.ticks1, args.tick1beg, args.tick1end,
                                                           args.tick1d, args.mtick1, x1beg, x1end, size1,
                                                           args.norm1, int(args.base1))
    if scalar1 != []:
        # add power to last tick if necessary
        tick1_label[-1] = tick1_label[-1] + r'$\mathregular{\times 10^{%i}}$' % scalar1
    ax.xaxis.set_ticks(tick1)
    ax.xaxis.set_ticklabels(tick1_label)

    # axis 2
    tick2, tick2_label, minor_tick2, scalar2 = define_tick(args.ticks2, args.tick2beg, args.tick2end,
                                                           args.tick2d, args.mtick2, x2beg, x2end, size2,
                                                           args.norm2, int(args.base2))
    if scalar2 != []:
        # add power to last tick if necessary
        tick2_label[-1] = tick2_label[-1] + r'$\mathregular{\times 10^{%i}}$' % scalar2
    ax.yaxis.set_ticks(tick2)
    ax.yaxis.set_ticklabels(tick2_label)

    # set tick font
    for l in ax.get_xticklabels():
        l.set_fontproperties(font)
        l.set_fontsize(tick_1_font_size)
    for l in ax.get_yticklabels():
        l.set_fontproperties(font)
        l.set_fontsize(tick_2_font_size)

    # turn off ticks if necessary
    if not args.tick1label:
        ax.xaxis.set_ticklabels([])
    if not args.tick2label:
        ax.yaxis.set_ticklabels([])

# set grid line
axx = 'x'
axy = 'y'
if args.direction == 'vertical':
    axx = 'y'
    axy = 'x'

if args.grid1:
    # grid line width
    if args.grid1width is None:
        grid1width = float(args.tickmajorwid)
    else:
        grid1width = float(args.grid1width)
    # add grid
    ax.grid(which='major', 
            axis=axx, 
            linestyle=args.grid1style, 
            color=args.grid1color, 
            linewidth=grid1width, 
            zorder=args.gridorder)

if args.grid2:
    # grid line width
    if args.grid2width is None:
        grid2width = float(args.tickmajorwid)
    else:
        grid2width = float(args.grid2width)
    # add grid
    ax.grid(which='major', 
            axis=axy, 
            linestyle=args.grid2style, 
            color=args.grid2color, 
            linewidth=grid2width, 
            zorder=args.gridorder)

if args.mgrid1 and args.mtick1 != 0:
    # minor grid line width
    if args.mgrid1width is None:
        mgrid1width = float(tick_minor_width)
    else:
        mgrid1width = float(args.mgrid1width)
    # add minor grid
    ax.grid(which='minor',
            axis=axx,
            linestyle=args.mgrid1style,
            color=args.mgrid1color,
            linewidth=mgrid1width, 
            zorder=args.gridorder)

if args.mgrid2 and args.mtick2 != 0:
    # minor grid line width
    if args.mgrid2width is None:
        mgrid2width = float(tick_minor_width)
    else:
        mgrid2width = float(args.mgrid2width)
    # add minor grid
    ax.grid(which='minor',
            axis=axy,
            linestyle=args.mgrid2style,
            color=args.mgrid2color,
            linewidth=mgrid2width, 
            zorder=args.gridorder)

# set title
set_title(args, fontbold)

# set annotation

# plot curve
if args.curve is not None:

    curvefile = args.curve[0].split(",")
    nf = len(curvefile)

    curvestyle = set_default(args.curvestyle, ',', nf, 'scatter.')
    curvecolor = set_default(args.curvecolor, ',', nf, 'k')
    curvefacecolor = set_default(args.curvefacecolor, ',', nf, 'k')
    curveedgecolor = set_default(args.curveedgecolor, ',', nf, 'none')
    curvesize = set_default(args.curvesize, ',', nf, 1.0, 'float')
    curveorder = set_default(args.curveorder, ',', nf, 9, 'int')
    curveselect = set_default(args.curveselect, ',', 2, 0, 'int')
    curvewidth = set_default(args.curvewidth, ',', nf, 1.0, 'float')

    for i in range(0, nf):
    
        curve = np.loadtxt(curvefile[i], ndmin=2)  # using ndmin=2 to ensure read as 2d array
        nsp = len(curve)
        s = curveselect
        px1 = curve[0:nsp, s[0] - 1]
        px2 = curve[0:nsp, s[1] - 1]
        if args.direction == 'vertical':
            temp = px1
            px1 = px2
            px2 = temp

        if args.plottype == 1 or args.plottype == 2:
            # plot scatter points on the figure
            if 'scatter' in curvestyle[i]:
                ax.scatter(px1,
                           px2,
                           marker=curvestyle[i][7:],
                           c=curvefacecolor[i],
                           zorder=curveorder[i],
                           edgecolor=curveedgecolor[i],
                           s=curvesize[i])
            # plot line on the figure
            if 'line' in curvestyle[i]:
                extra = Line2D(px1,
                               px2,
                               linestyle=curvestyle[i][4:],
                               zorder=curveorder[i],
                               color=curvecolor[i],
                               linewidth=curvesize[i])
                ax.add_artist(extra)
            # plot polygon on the figure
            if 'polygon' in curvestyle[i]:
                curve[0:nsp, 0] = px1
                curve[0:nsp, 1] = px2
                extra = Polygon(curve,
                                fill=False,
                                zorder=curveorder[i],
                                color=curvecolor[i],
                                linewidth=curvesize[i],
                                antialiased=True)
                ax.add_artist(extra)
        else:
            # plot scatter points on the figure
            if 'scatter' in curvestyle[i]:
                ax.scatter(px1,
                           px2,
                           marker=curvestyle[i][7:],
                           facecolor=curvefacecolor[i],
                           zorder=curveorder[i],
                           edgecolor=curveedgecolor[i],
                           s=curvesize[i],
                           antialiased=True,
                           label=plotlabel[i])
            
            # plot line on the figure
            if 'line' in curvestyle[i]:
                extra = Line2D(px1,
                               px2,
                               linestyle=curvestyle[i][4:],
                               zorder=curveorder[i],
                               color=curvecolor[i],
                               linewidth=curvewidth[i],
                               antialiased=True,
                               label=plotlabel[i])
                ax.add_artist(extra)
            # plot polygon on the figure
            if 'polygon' in curvestyle[i]:
                curve[0:nsp, 0] = px1
                curve[0:nsp, 1] = px2
                extra = Polygon(curve,
                                fill=False,
                                zorder=curveorder[i],
                                edgecolor=curvecolor[i],
                                linewidth=curvewidth[i],
                                antialiased=True,
                                label=plotlabel[i])
                ax.add_artist(extra)
            # add legend
            if args.plotlabel is not None:
                if args.plotlabelloc in list(locdict.keys()):
                    lg = plt.legend(loc=labelloc)
                else:
                    lg = plt.legend(bbox_to_anchor=(1.01, 1.0), loc=2, borderaxespad=0)
                lg.set_zorder(10)
                leg = ax.get_legend()
                ltext = leg.get_texts()
                plt.setp(ltext, fontproperties=font)
                plt.setp(ltext, fontsize=float(args.plotlabelsize))

# place text
if args.text is not None:

    # text contents
    text = args.text[0].split(":")
    nf = len(text)

    textcolor = set_default(args.textcolor, ',', nf, 'k')
    textsize = set_default(args.textsize, ',', nf, 14.0, 'float')
    textrotation = set_default(args.textrotation, ',', nf, 0, 'float')
    textbox = set_default(args.textbox, ',', nf, 0, 'int')
    textboxedgecolor = set_default(args.textboxedgecolor, ',', nf, 'k')
    textboxfacecolor = set_default(args.textboxfacecolor, ',', nf, 'w')
    textboxalpha = set_default(args.textboxalpha, ',', nf, 0.5, 'float')
    textboxstyle = set_default(args.textboxstyle, ',', nf, 'square')
    textboxpad = set_default(args.textboxpad, ',', nf, 10.0, 'float')
    textorder = set_default(args.textorder, ',', nf, 10, 'int')

    # text location
    # default coordinates according to data center
    center1 = 0.5 * (x1beg + x1end)
    center2 = 0.5 * (x2beg + x2end)
    dtextloc = [center1 for i in range(0, 2 * nf)]
    # even position set to center2
    for i in range(1, 2 * nf, 2):
        dtextloc[i] = center2
    dtextloc = reshape(dtextloc, (nf, 2))

    if args.textloc is not None:
        textloc = args.textloc[0].split(":")
        if len(textloc) < nf:
            l = len(textloc)
            for i in range(0, l):
                dtextloc[i, :] = textloc[i].split(',')
            textloc = dtextloc
        else:
            for i in range(0, nf):
                textloc[i] = textloc[i].split(',')
    else:
        textloc = dtextloc

    # iterate through all text annotation
    for i in range(0, nf):

        tloc = textloc[i]
        px1 = float(eval(str(tloc[0])))
        px2 = float(eval(str(tloc[1])))
        if args.direction == 'vertical':
            temp = px1
            px1 = px2
            px2 = temp

        if int(textbox[i]) == 1:
            # if textbox is required
            t = ax.text(px1,
                        px2,
                        text[i],
                        zorder=textorder[i],
                        color=textcolor[i],
                        fontproperties=font,
                        fontsize=textsize[i],
                        horizontalalignment='center',
                        verticalalignment='center',
                        rotation=textrotation[i],
                        bbox={
                            'boxstyle': textboxstyle[i],
                            'facecolor': textboxfacecolor[i],
                            'edgecolor': textboxedgecolor[i],
                            'alpha': textboxalpha[i]
                        })
            bb = t.get_bbox_patch()
            bb.set_boxstyle(pad=float(textboxpad[i]))

        else:
            # plain text only
            ax.text(px1,
                    px2,
                    text[i],
                    zorder=textorder[i],
                    color=textcolor[i],
                    fontproperties=font,
                    fontsize=textsize[i],
                    horizontalalignment='center',
                    verticalalignment='center',
                    rotation=textrotation[i])

# add arrows
if args.arrow is not None:

    # arrow start and ending coordinates
    arrow = args.arrow[0].split(':')
    nf = len(arrow)

    arrowfacecolor = set_default(args.arrowfacecolor, ',', nf, 'k')
    arrowedgecolor = set_default(args.arrowedgecolor, ',', nf, 'k')
    arrowstyle = set_default(args.arrowstyle, ',', nf, '->')
    arrowlinestyle = set_default(args.arrowlinestyle, ',', nf, 'solid')
    arrowconnect = set_default(args.arrowconnect, ':', nf, 'arc3')
    arrowwidth = set_default(args.arrowwidth, ',', nf, 1.0, 'float')
    arroworder = set_default(args.arroworder, ',', nf, 9, 'int')

    for i in range(0, nf):

        arrowloc = arrow[i].split(',')
        snf = len(arrowloc)

        if snf % 2 != 0:
            print('arrow tail/head coordinates specification error')
            exit()

        tail1 = float(arrowloc[0])
        tail2 = float(arrowloc[1])
        if args.direction == 'vertical':
            temp = tail1
            tail1 = tail2
            tail2 = temp
        head1 = float(arrowloc[2])
        head2 = float(arrowloc[3])
        if args.direction == 'vertical':
            temp = head1
            head1 = head2
            head2 = temp
        if arrowstyle[i] != '-':
            ax.annotate('',
                        xytext=(tail1, tail2),
                        xy=(head1, head2),
                        arrowprops=dict(arrowstyle=arrowstyle[i],
                                        connectionstyle=arrowconnect[i],
                                        facecolor=arrowfacecolor[i],
                                        edgecolor=arrowedgecolor[i],
                                        linewidth=arrowwidth[i],
                                        linestyle=arrowlinestyle[i]),
                        zorder=arroworder[i])
        else:
            x1 = [tail1, head1]
            x2 = [tail2, head2]
            extra = Line2D(x2,
                           x1,
                           linestyle=arrowlinestyle[i],
                           zorder=arroworder[i],
                           color=arrowfacecolor[i],
                           linewidth=arrowwidth[i])
            ax.add_artist(extra)

# add circle/ellipse
if args.circle is not None:

    circle = args.circle[0].split(':')
    nf = len(circle)

    circlefacecolor = set_default(args.circlefacecolor, ',', nf, 'y')
    circleedgecolor = set_default(args.circleedgecolor, ',', nf, 'y')
    circlealpha = set_default(args.circlealpha, ',', nf, 0.5, 'float')
    circlelinestyle = set_default(args.circlelinestyle, ',', nf, '-')
    circlelinewidth = set_default(args.circlelinewidth, ',', nf, 1.0, 'float')
    circleorder = set_default(args.circleorder, ',', nf, 8, 'int')

    # iterate through each circle
    for i in range(0, nf):

        scircle = circle[i].split(',')
        snf = len(scircle)

        if snf == 5:
            c1 = scircle[0]
            c2 = scircle[1]
            if args.direction == 'vertical':
                temp = c1
                c1 = c2
                c2 = temp
            r1 = scircle[2]
            r2 = scircle[3]
            angle = scircle[4]
        else:
            print('circle specification error')
            exit()

        extra = Ellipse(xy=[float(c1), float(c2)],
                        fill=True,
                        facecolor=circlefacecolor[i],
                        edgecolor=circleedgecolor[i],
                        alpha=circlealpha[i],
                        width=float(r1),
                        height=float(r2),
                        angle=float(angle),
                        linestyle=circlelinestyle[i],
                        linewidth=circlelinewidth[i],
                        zorder=circleorder[i],
                        antialiased=True)
        ax.add_artist(extra)

# add filled polygons
if args.polygon is not None:

    polygon = args.polygon[0].split(':')
    nf = len(polygon)

    polygonfacecolor = set_default(args.polygonfacecolor, ',', nf, 'y')
    polygonedgecolor = set_default(args.polygonedgecolor, ',', nf, 'y')
    polygonalpha = set_default(args.polygonalpha, ',', nf, 0.5, 'float')
    polygonlinestyle = set_default(args.polygonlinestyle, ',', nf, '-')
    polygonlinewidth = set_default(args.polygonlinewidth, ',', nf, 1.0, 'float')
    polygonorder = set_default(args.polygonorder, ',', nf, 7, 'int')

    # iterate through each polygon
    for i in range(0, nf):

        spolygon = polygon[i].split(',')
        snf = len(spolygon)

        if snf % 2 == 0:
            spolygon = reshape(spolygon, (int(snf / 2.0), 2))
        else:
            print('polygon specification error')
            exit()
        snf = int(snf / 2.0)

        spolygon = np.asanyarray(spolygon, dtype=np.float32)

        x1 = spolygon[0:snf, 0]
        x2 = spolygon[0:snf, 1]

        for j in range(0, snf):
            px1 = float(x1[j])
            px2 = float(x2[j])
            if args.direction == 'vertical':
                temp = px1
                px1 = px2
                px2 = temp
            spolygon[j, 0] = px1
            spolygon[j, 1] = px2

        extra = Polygon(spolygon,
                        fill=True,
                        facecolor=polygonfacecolor[i],
                        edgecolor=polygonedgecolor[i],
                        alpha=polygonalpha[i],
                        linestyle=polygonlinestyle[i],
                        linewidth=polygonlinewidth[i],
                        zorder=polygonorder[i],
                        antialiased=True)
        ax.add_artist(extra)


# set colorbar
def set_colorbar(args, im, font, plot_min_value, plot_max_value, figheight, figwidth, fig):

    if args.legend and plot_min_value != plot_max_value:

        # legend location
        if args.lloc is None:
            lloc = 'right'
        else:
            lloc = args.lloc

        # legend orientation and tick switc, unit rotation
        if lloc == 'left':
            lorient = 'vertical'
            lleft = 1
            lright = 0
            ltop = 0
            lbottom = 0
            lrotate = 90
        if lloc == 'right':
            lorient = 'vertical'
            lleft = 0
            lright = 1
            ltop = 0
            lbottom = 0
            lrotate = 270
        if lloc == 'top':
            lorient = 'horizontal'
            lleft = 0
            lright = 0
            ltop = 1
            lbottom = 0
            lrotate = 0
        if lloc == 'bottom':
            lorient = 'horizontal'
            lleft = 0
            lright = 0
            ltop = 0
            lbottom = 1
            lrotate = 0

        if args.unitpad is None:
            if lloc == 'right':
                lpad = 20.0
            if lloc == 'bottom':
                lpad = 0.0
            if lloc == 'left':
                lpad = -50
            if lloc == 'top':
                lpad = -50.0
        else:
            lpad = float(args.unitpad)

        if lloc == 'left' or lloc == 'right':
            if args.lheight is None:
                lheight = figheight
            else:
                lheight = float(args.lheight)
            if args.lwidth is None:
                lwidth = 0.15
            else:
                lwidth = float(args.lwidth)
        if lloc == 'top' or lloc == 'bottom':
            if args.lheight is None:
                lheight = 0.15
            else:
                lheight = float(args.lheight)
            if args.lwidth is None:
                lwidth = figwidth
            else:
                lwidth = float(args.lwidth)

        if args.lpad is None:
            cbpad = 0.05
        else:
            cbpad = float(args.lpad)
        if lloc == 'left':
            cbx = -(cbpad + lwidth) / figwidth
            cby = (figheight - lheight) / 2.0 / figheight
        if lloc == 'right':
            cbx = 1.0 + cbpad / figwidth
            cby = (figheight - lheight) / 2.0 / figheight
        if lloc == 'top':
            cbx = (figwidth - lwidth) / 2.0 / figwidth
            cby = 1.0 + cbpad / figheight
        if lloc == 'bottom':
            cbx = (figwidth - lwidth) / 2.0 / figwidth
            cby = -(cbpad + lheight) / figheight
        cax = fig.add_axes([cbx, cby, lwidth / figwidth, lheight / figheight])
        cb = plt.colorbar(im, cax=cax, orientation=lorient)

        # set colorbar label and styles
        if args.unitsize is None:
            lufs = min(float(args.label1size), float(args.label2size)) - 1
        else:
            lufs = float(args.unitsize)

        if args.unit is None:
            legend_units = ' '
        else:
            legend_units = args.unit

        cb.set_label(legend_units, rotation=lrotate, labelpad=lpad, fontproperties=font)
        if lloc == 'right' or lloc == 'left':
            cb.ax.yaxis.label.set_fontsize(lufs)
        if lloc == 'top' or lloc == 'bottom':
            cb.ax.xaxis.label.set_fontsize(lufs)

        # set colorbar tick styles: font size, family, and ticks direction
        if args.lticksize is None:
            ltfs = lufs - 1
        else:
            ltfs = float(args.lticksize)

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

        for l in cb.ax.yaxis.get_ticklabels():
            l.set_fontproperties(font)
            l.set_fontsize(ltfs)
        for l in cb.ax.xaxis.get_ticklabels():
            l.set_fontproperties(font)
            l.set_fontsize(ltfs)

        # linear data norm
        if args.norm == 'linear':

            # set colorbar major ticks
            if args.ld is None:
                ld = nice((plot_max_value - plot_min_value) / 5.0)
            else:
                ld = float(args.ld)

            if args.ltickbeg is None:
                ltickbeg = nice(plot_min_value)
                base = 0.5
                nb = 0
                while nb <= 10 and ltickbeg > plot_min_value + ld:
                    base = base / 10.0
                    ltickbeg = nice(plot_min_value, base)
                    nb = nb + 1
                if abs(ltickbeg) < abs(plot_max_value) and orderm(ltickbeg) + 2 < orderm(plot_max_value):
                    ltickbeg = 0.0
            else:
                ltickbeg = float(args.ltickbeg)
            if args.ltickend is None:
                ltickend = plot_max_value
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
            pminv = plot_min_value
            pmaxv = plot_max_value

            tbeg = max(pminv, ltickbeg)
            tend = min(pmaxv, ltickend)

            # set tick positions on colorbar
            ticks = np.asarray(
                [i for i in ticks if i >= tbeg - 1.0e-10 * abs(tbeg) and i <= tend + 1.0e-10 * abs(tend)])
            if lloc in ['left', 'right']:
                cb.ax.yaxis.set_ticks(ticks)
            else:
                cb.ax.xaxis.set_ticks(ticks)

            # set tick labels on colorbar
            tick_labels = ['' for i in range(0, len(ticks))]
            for i in range(0, len(ticks)):
                tick_labels[i] = ('%f' % (ticks[i] / cscale)).rstrip('0').rstrip('.')
            # cb.set_ticklabels(tick_labels)
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
                mticks = np.asarray([
                    i for i in mticks if i >= tbeg - 1.0e-10 * abs(tbeg) and i <= tend + 1.0e-10 * abs(tend)
                ])
                # set minor ticks
                if lloc == 'left' or lloc == 'right':
                    cb.ax.yaxis.set_ticks(mticks, minor=True)
                else:
                    cb.ax.xaxis.set_ticks(mticks, minor=True)

        # log data norm
        if args.norm == 'log':

            # set colorbar major ticks
            if args.ltickbeg is None:
                ltickbeg = np.floor(plot_min_value)
            else:
                ltickbeg = float(args.ltickbeg)
            if args.ltickend is None:
                ltickend = np.ceil(plot_max_value) + 1
            else:
                ltickend = float(args.ltickend)
            if args.ld is None:
                ld = max(1, round((ltickend - ltickbeg) / 5.0))
            else:
                ld = int(args.ld)

            ticks = np.arange(ltickbeg, ltickend + 1, ld)
            pminv = plot_min_value
            pmaxv = plot_max_value

            tbeg = max(pminv, ltickbeg)
            tend = min(pmaxv, ltickend)

            # set tick positions on colorbar
            ticks = np.asarray(
                [i for i in ticks if i >= tbeg - 1.0e-10 * abs(tbeg) and i <= tend + 1.0e-10 * abs(tend)])
            if lloc in ['left', 'right']:
                cb.ax.yaxis.set_ticks((ticks - pminv) / (pmaxv - pminv))
            else:
                cb.ax.xaxis.set_ticks((ticks - pminv) / (pmaxv - pminv))

            # set tick labels on colorbar
            tick_labels = ['' for i in range(0, len(ticks))]
            for i in range(0, len(ticks)):
                tick_labels[i] = r'$\mathregular{10^{%i}}$' % (ticks[i])
            if lloc in ['left', 'right']:
                cb.ax.yaxis.set_ticklabels(tick_labels)
            else:
                cb.ax.xaxis.set_ticklabels(tick_labels)

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
                    mticks = np.append(
                        mticks, np.log10(np.linspace(10**pticks[i], 10**pticks[i + 1], args.lmtick + 2)))
                mticks = [i for i in mticks if (i not in pticks)]
                mticks = np.asarray([
                    i for i in mticks if i >= tbeg - 1.0e-10 * abs(tbeg) and i <= tend + 1.0e-10 * abs(tend)
                ])
                # set minor ticks
                if lloc == 'left' or lloc == 'right':
                    cb.ax.yaxis.set_ticks((mticks - pminv) / (pmaxv - pminv), minor=True)
                else:
                    cb.ax.xaxis.set_ticks((mticks - pminv) / (pmaxv - pminv), minor=True)

        # make colorbar solid color continuous
        cb.solids.set_edgecolor("face")

        # colorbar reverse
        if lloc in ['left', 'right']:
            if args.lreverse == 1:
                cb.ax.invert_yaxis()
        else:
            if args.lreverse == 1:
                cb.ax.invert_xaxis()


if args.legend and args.plottype == 3:
    set_colorbar(args, sc, font, cmin, cmax, figheight, figwidth, fig)

# # background
# if args.background is not None:
    
#     bn1 = args.bn1
    
#     fsize = os.path.getsize(args.background)
#     datatype = args.datatype
#     if datatype == 'double':
#         fsize = fsize / 8
#     if datatype == 'float':
#         fsize = fsize / 4
#     if datatype == 'int':
#         fsize = fsize / 2

#     if args.bn2 is None:
#         bn2 = int(fsize * 1.0 / bn1)
#     else:
#         bn2 = args.bn2
    
#     w = np.fromfile(args.background, count=bn1*bn2, dtype=np.float32)
#     w = np.reshape(w, (bn2, bn1))
#     w = np.transpose(w)
    
#     bd1 = float(args.bd1)
#     bd2 = float(args.bd2)
    
#     backdmin = w.min()
#     backdmax = w.max()
    
#     sp1beg, sp1end, x1beg, x1end, n1beg, n1end = set_range(args.bo1, bn1, bd1, args.x1beg, args.x1end)
#     sp2beg, sp2end, x2beg, x2end, n2beg, n2end = set_range(args.bo2, bn2, bd2, args.x2beg, args.x2end)
#     w = w[n1beg:n1end, n2beg:n2end]
    
#     from module_colormap import set_colormap, set_colormap_alpha

#     backcmin, backcmax = set_clip(args, w, 'back', backdmin, backdmax)
#     colormap = set_colormap(args, 'background')
#     colormap = set_colormap_alpha(args, colormap, backcmin, backcmax, 'background')
    
#     print(ax)
    
#     if args.shading is not None:
#         from matplotlib.colors import LightSource
#         ls = LightSource(azdeg=360, altdeg=45)
#         rgb = ls.hillshade(data, vert_exag=args.shading_scale) #, blend_mode=args.shading)
#         ax.imshow(rgb, cmap=colormap, extent=[x2beg, x2end, x1beg, x1end], zorder=1)
#     else:
#         ax.imshow(data, cmap=colormap)
        

# clean
# if multiple input files then remove temporary file
if ninput > 1:
    os.remove(tempfile)

# output
output(args)

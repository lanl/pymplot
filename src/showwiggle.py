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

# this is to ignore warnings
import warnings
warnings.filterwarnings("ignore", module="matplotlib")

## read arguments
# assign description to the help doc
parser=argparse.ArgumentParser(description= \
    '''Read a 2D array from binary file and plot as wiggles, 
    written by K.G. @ 2016.05, 2016.06, 2016.08, 2016.10''')

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

# arguments -- wiggle
parser = getarg_wiggle(parser)
parser.add_argument('-plotlabel',
                    '--plotlabel',
                    type=str,
                    help='Legend labels for each plot',
                    required=False,
                    nargs='+',
                    default='')
parser.add_argument('-plotlabelloc',
                    '--plotlabelloc',
                    type=str,
                    help='Legend label location',
                    required=False,
                    default='upper_right')
parser.add_argument('-plotlabelsize',
                    '--plotlabelsize',
                    type=str,
                    help='Legend label font size',
                    required=False,
                    default='14.0')
parser.add_argument('-wiggleloc',
                    '--wiggleloc',
                    type=str,
                    help='ASCII file containing the spatial location of traces',
                    required=False,
                    default='')
parser.add_argument('-wx1beg',
                    '--wx1beg',
                    type=str,
                    help='Axis 1 begin value, in consistent with the real value of axis',
                    required=False,
                    default='')
parser.add_argument('-wx1end',
                    '--wx1end',
                    type=str,
                    help='Axis 1 end value, in consistent with the real value of axis',
                    required=False,
                    default='')
parser.add_argument('-wx2beg',
                    '--wx2beg',
                    type=str,
                    help='Axis 2 begin value, in consistent with the real value of axis',
                    required=False,
                    default='')
parser.add_argument('-wx2end',
                    '--wx2end',
                    type=str,
                    help='Axis 2 end value, in consistent with the real value of axis',
                    required=False,
                    default='')

# array for all arguments passed to script
args = parser.parse_args()

## input data
infile = args.infile[0].split(',')
nf = len(infile)

for i in range(0, nf):
    if not os.path.exists(infile[i]):
        print()
        print('input file', infile[i], 'does not exists')
        print()
        exit()

fsize = os.path.getsize(infile[0])
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

d1 = float(args.d1)
d2 = float(args.d2)

# if wiggle locations are read from a file
if len(args.wiggleloc) != 0:
    # from ASCII file
    wloc = np.loadtxt(args.wiggleloc, ndmin=2)
    [i for i in wloc if i]
    wloc = np.array(wloc)
    shape = shape(wloc)
    if args.along == 1 and n2 != shape[0]:
        print(" Error: n2 = ", n2, ' != size(wiggleloc) = ', shape[0])
        exit()
    if args.along == 2 and n1 != shape[0]:
        print(" Error: n1 = ", n1, ' != size(wiggleloc) = ', shape[0])
        exit()

    # n, o, d
    if args.along == 1:
        n2 = shape[0]
        f2 = wloc[0][0]
        d2 = (wloc.max() - wloc.min()) / (n2 - 1.0)

        n1 = n1
        f1 = float(args.f1)
        d1 = d1
    else:
        n1 = shape[0]
        f1 = wloc[0][0]
        d1 = (wloc.max() - wloc.min()) / (n1 - 1.0)

        n2 = n2
        f2 = float(args.f2)
        d2 = d2

else:
    f1 = float(args.f1)
    f2 = float(args.f2)

# print(n1, d1, f1)
# print(n2, d2, f2)
# print(wloc)

## limit of axis
from module_range import *
sp1beg, sp1end, x1beg, x1end, n1beg, n1end = set_range(f1, n1, d1, args.x1beg, args.x1end)
sp2beg, sp2end, x2beg, x2end, n2beg, n2end = set_range(f2, n2, d2, args.x2beg, args.x2end)

## plot
print()
for i in range(0, nf):
    print('input <<    ', infile[i])

# data type
from module_datatype import *
dt = set_datatype(args)

adata = np.empty([nf, n1end - n1beg, n2end - n2beg])
for i in range(0, nf):

    # read binary data
    data = fromfile(infile[i], dtype=dt, count=n1 * n2)
    if args.transpose == 0:
        data = data.reshape((n2, n1))
        data = data.transpose()
    else:
        data = data.reshape(n1, n2)

    if i == 0:
        print('shape       ', data.shape)

    # print value range of input files
    if isnan(sum(data)) == True:
        udata = data[~isnan(data)]
        if udata.shape == (0, ):
            print('error: input dataset is nan')
            exit()
        else:
            dmin = udata.min()
            dmax = udata.max()
    else:
        dmin = data.min()
        dmax = data.max()
    print('value range ', dmin, ' -- ', dmax)

    # crop data
    data = data[n1beg:n1end, n2beg:n2end]

    # assign to whole data
    adata[i, :, :] = data

# read background data
if len(args.background) != 0:

    # check if background file exists
    if not os.path.exists(args.background):
        print()
        print('input file', args.background, 'does not exists')
        print()
        exit()

    # input
    backdata = np.empty([n1, n2])
    backdata = fromfile(args.background, dtype=dt, count=n1 * n2)

    # transpose
    if args.transpose == 0:
        backdata = backdata.reshape((n2, n1))
        backdata = backdata.transpose()
    else:
        backdata = backdata.reshape(n1, n2)

    # crop
    backdata = backdata[n1beg:n1end, n2beg:n2end]

## set clip
from module_clip import *
cmin, cmax = set_clip(args, adata, 'fore')
if len(args.background) != 0:
    backcmin, backcmax = set_clip(args, backdata, 'back')

## figure size
from module_size import *
figheight, figwidth = set_size(args, n1beg, n1end, n2beg, n2end)

## begin plot
fig = plt.figure(figsize=(figwidth, figheight))
ax = fig.add_axes([0, 0, 1, 1])

## set frame
from module_frame import *
set_frame(args)

## plot image
# show image if necessary
if nf == 1 and args.overlay == 1 and len(args.background) == 0:

    # show data by imshow
    im = ax.imshow(adata[0, :, :])  #,alpha=float(args.alpha))

    # set colormap
    from module_colormap import set_colormap
    colormap = set_colormap(args)
    im.set_cmap(colormap)

    # set clip
    im.set_clim(cmin, cmax)

    # set interpolation
    im.set_interpolation(args.interp)

    #alpha = float(args.alpha)

# plot background image if necessary
if args.overlay == 0 and len(args.background) != 0:

    # beg plot
    im = ax.imshow(backdata)  #,alpha=float(args.backalpha))

    # set colormap
    from module_colormap import set_colormap
    colormap = set_colormap(args, 'background')
    im.set_cmap(colormap)

    # set clip
    im.set_clim(backcmin, backcmax)

    # set interpolation
    im.set_interpolation(args.interp)

## set font
from module_font import *
font, fontbold = set_font(args)

## set tick
from module_tick import *
set_tick(args, font, x1beg, x1end, n1beg, n1end, d1, figheight, x2beg, x2end, n2beg, n2end, d2, figwidth, extend=True)

## set grid line
from module_gridline import *
set_gridline(args)

## set title
from module_title import *
set_title(args, fontbold)

## set annotation
from module_annotation import *
set_annotation(args, font, x1beg, n1end - n1beg, d1, figheight, x2beg, n2end - n2beg, d2, figwidth)

## plot wiggles
ax = plt.gca()

# scaling factors
if len(args.wiggleloc) != 0 and args.along == 2:
    scale1 = wloc.max() - wloc.min() + d1
    scale1 = scale1 / figheight
else:
    scale1 = (n1end - n1beg) * d1 / figheight
if scale1 == 0:
    scale1 = 1.0

if len(args.wiggleloc) != 0 and args.along == 1:
    scale2 = wloc.max() - wloc.min() + d2
    scale2 = scale2 / figwidth
else:
    scale2 = (n2end - n2beg) * d2 / figwidth
if scale2 == 0:
    scale2 = 1.0

# setup line colors, if given number of colors < number of files, then cycle append colors
from itertools import cycle
# defaultcolor=cycle('brgykcp')
defaultcolor = cycle(['blue', 'red', 'green', 'yellow', 'black', 'cyan', 'magenta'])
color = args.wigglecolor[0].split(',')

nc = len(color)
if nc < nf:
    ic = 0
    for i in cycle(defaultcolor):
        color.append(i)
        ic = ic + 1
        if ic >= nf:
            break

wmin = 1.0e10
wmax = 0.0

from module_annotation import set_default

linewidth = set_default(args.wigglewidth, ',', nf, 1.0, 'float')
linestyle = set_default(args.wigglestyle, ',', nf, '-')

# plot labels if necessary
if len(args.plotlabel) != 0:
    plotlabel = args.plotlabel[0].split(',')
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

# start wiggle plot
from scipy.interpolate import InterpolatedUnivariateSpline

if args.along == 1:

    if x1beg == x1end:
        print()
        print(' Error: Wiggles are along axis 1, but have only one sample. ')
        print()
        exit()

    traceinterval = int(args.every) * abs(d2) / scale2
    adata = adata / (cmax - cmin) * traceinterval
    traces = np.arange(0, n2end - n2beg, int(args.every))
    y = (np.arange(0, n1end - n1beg) * d1) / scale1

    if len(args.interp1) != 0:
        yy = np.linspace(y.min(), y.max(), int((n1end - n1beg) * float(args.interp1) + 1))
    else:
        yy = y

    for i in traces:

        if len(args.wiggleloc) != 0:
            offset = (wloc[i][0] - wloc[0][0] + 0.5 * d2) / scale2
        else:
            offset = (i * d2 + 0.5 * d2) / scale2

        # iterate through all datasets
        for j in range(0, nf):

            # select data
            data = adata[j, :, :]

            # plot data
            x = data[:, i] + offset

            if len(args.interp1) != 0:
                spl = InterpolatedUnivariateSpline(y, x, k=3)
                xx = spl(yy)
            else:
                xx = x

            xx[where(xx >= offset + traceinterval)] = offset + traceinterval
            xx[where(xx <= offset - traceinterval)] = offset - traceinterval

            # wiggles
            if i != traces[-1]:
                plt.plot(xx, yy, color=color[j], linewidth=linewidth[j], linestyle=linestyle[j], antialiased=True)
            else:
                plt.plot(xx,
                         yy,
                         color=color[j],
                         linewidth=linewidth[j],
                         linestyle=linestyle[j],
                         antialiased=True,
                         label=plotlabel[j])

            # fill positive/negative polarity if necessary
            if args.fill == 1:
                ax.fill_betweenx(yy,
                                 offset,
                                 xx,
                                 interpolate=True,
                                 antialiased=True,
                                 lw=0,
                                 where=(xx > offset),
                                 color=color[j],
                                 edgecolor='none')
            if args.fill == -1:
                ax.fill_betweenx(yy,
                                 offset,
                                 xx,
                                 interpolate=True,
                                 antialiased=True,
                                 lw=0,
                                 where=(xx < offset),
                                 color=color[j],
                                 edgecolor='none')

            if i == traces[0]:
                wmin = min(wmin, xx.min())
            if i == traces[-1]:
                wmax = max(wmax, xx.max())

    # to ensure correct image show range
    omin = (traces[0] * d2) / scale2
    omax = (traces[-1] * d2 + d2) / scale2

    # invert y axis for make zero at the top
    if args.reverse1 == 0:
        ax.invert_yaxis()

if args.along == 2:

    if x2beg == x2end:
        print()
        print(' Error: Wiggles are along axis 2, but have only one sample. ')
        print()
        exit()

    traceinterval = int(args.every) * abs(d1) / scale1
    adata = adata / (cmax - cmin) * traceinterval
    traces = np.arange(0, n1end - n1beg, int(args.every))
    y = (np.arange(0, n2end - n2beg) * d2) / scale2

    if len(args.interp2) != 0:
        yy = np.linspace(y.min(), y.max(), int((n2end - n2beg) * float(args.interp2) + 1))
    else:
        yy = y

    for i in traces:

        if len(args.wiggleloc) != 0:
            offset = (wloc[i][0] - wloc[0][0] + 0.5 * d1) / scale1
        else:
            offset = (i * d1 + 0.5 * d1) / scale1

        # iterate through all datasets
        for j in range(0, nf):

            # select data
            data = adata[j, :, :]

            # plot data
            x = data[i, :] + offset

            if len(args.interp2) != 0:
                spl = InterpolatedUnivariateSpline(y, x, k=3)
                xx = spl(yy)
            else:
                xx = x

            xx[where(xx >= offset + traceinterval)] = offset + traceinterval
            xx[where(xx <= offset - traceinterval)] = offset - traceinterval

            # wiggles
            if i != traces[-1]:
                plt.plot(yy, xx, color=color[j], linewidth=linewidth[j], linestyle=linestyle[j], antialiased=True)
            else:
                plt.plot(yy,
                         xx,
                         color=color[j],
                         linewidth=linewidth[j],
                         linestyle=linestyle[j],
                         antialiased=True,
                         label=plotlabel[j])

            # fill positive/negative polarity if necessary
            if args.fill == 1:
                ax.fill_between(yy,
                                offset,
                                xx,
                                interpolate=True,
                                antialiased=True,
                                lw=0,
                                where=(xx > offset),
                                facecolor=color[j],
                                edgecolor='none')
            if args.fill == -1:
                ax.fill_between(yy,
                                offset,
                                xx,
                                interpolate=True,
                                antialiased=True,
                                lw=0,
                                where=(xx < offset),
                                facecolor=color[j],
                                edgecolor='none')

            if i == traces[0]:
                wmin = min(wmin, xx.min())
            if i == traces[-1]:
                wmax = max(wmax, xx.max())

    # to ensure correct image show range
    omin = (traces[0] * d1) / scale1
    omax = (traces[-1] * d1 + d1) / scale1

# add plot labels
if len(args.plotlabel) != 0:
    if args.plotlabelloc in list(locdict.keys()):
        lg = plt.legend(loc=labelloc)
    else:
        lg = plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0)
    lg.set_zorder(10)
    lg.get_frame().set_alpha(1)
    leg = ax.get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontproperties=font)
    plt.setp(ltext, fontsize=float(args.plotlabelsize))

## reset figure sizes
if (nf==1 and args.overlay==1 and len(args.background)==0) or \
    (args.overlay==0 and len(args.background)!=0):
    if args.along == 1:
        im.set_extent([omin, omax, figheight, 0])
    else:
        im.set_extent([0, figwidth, omax, omin])

apad = float(args.axispad)

if args.along == 1:
    wmin = wmin - apad
    wmax = wmax + apad
    if len(args.wx2beg) != 0:
        wmin = (float(args.wx2beg) - wloc[0][0] + 0.5 * d2) / scale2
    if len(args.wx2end) != 0:
        wmax = (float(args.wx2end) - wloc[0][0] + 0.5 * d2) / scale2
    ax.set_xlim([wmin, wmax])
    ax.set_ylim([figheight, 0])
    ax.set_aspect('auto')
if args.along == 2:
    wmin = wmin - apad
    wmax = wmax + apad
    if len(args.wx1beg) != 0:
        wmin = (float(args.wx1beg) - wloc[0][0] + 0.5 * d1) / scale1
    if len(args.wx1end) != 0:
        wmax = (float(args.wx1end) - wloc[0][0] + 0.5 * d1) / scale1
    ax.set_xlim([0, figwidth])
    if args.reverse1 == 0:
        ax.set_ylim([wmin, wmax])
    else:
        ax.set_ylim([wmax, wmin])
    ax.set_aspect('auto')

## axis invert
if args.reverse1 == 1:
    ax.invert_yaxis()
if args.reverse2 == 1:
    ax.invert_xaxis()

## set colorbar
if nf == 1 and args.overlay == 1 and args.legend == 1:
    from module_colorbar import set_colorbar
    set_colorbar(args, im, font, cmin, cmax, figheight, figwidth, fig)

## output
from module_output import *
output(args)

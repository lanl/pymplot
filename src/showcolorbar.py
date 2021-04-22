#!/usr/bin/python3

import matplotlib as mplt
import matplotlib.pyplot as plt
from module_utility import *
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
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
parser = argparse.ArgumentParser(description='''Plot a colorbar, K.G. @ 2018.12''',
                                 formatter_class=RawTextHelpFormatter)

# arguments -- general
parser.add_argument('-norm',
                    '--norm',
                    type=str,
                    help='''Norm of data values, =linear (default) or =log;
for log norm, the base is 10, and values of the data
will be converted to logarithmic norm,
the options cmin/cmax should be based on the
power of the converted data,
rather than based on the values of the raw data''',
                    required=False,
                    default='linear')
parser.add_argument('-font',
                    '--font',
                    type=str,
                    help='''Font style for the plot; valid options,
=arial (Arial) by default,
=times (Times New Roman),
=courier (Courier New)',
=helvetica (Helvetica),
=georgia (Georgia),
=consolas (Consolas) ''',
                    required=False,
                    default='arial')

# arguments -- output
parser.add_argument('-dpi',
                    '--dpi',
                    type=str,
                    help='DPI of generated figure, =300 (default)',
                    required=False,
                    default='300')
parser.add_argument('-imageonly',
                    '--imageonly',
                    type=int,
                    help='Save only plotting region, no frame or axes',
                    required=False,
                    default=0)
parser.add_argument('-out',
                    '--outfile',
                    type=str,
                    help='Output figure file',
                    nargs='+',
                    required=False,
                    default='')

# arguments -- color
parser = getarg_color(parser)

# arguments -- colorbar
parser = getarg_colorbar(parser)

# array for all arguments passed to script
args = parser.parse_args()

print()

# set font
from module_font import *
font, fontbold = set_font(args)

# default plot range
if len(args.cmin) == 0:
    plot_min_value = 0
else:
    plot_min_value = float(args.cmin)

if len(args.cmax) == 0:
    plot_max_value = 255
else:
    plot_max_value = float(args.cmax)

if len(args.lloc) == 0:
    lloc = right
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

figheight = 4
figwidth = 4
if lloc == 'left' or lloc == 'right':
    if len(args.lheight) == 0:
        lheight = figheight
    else:
        lheight = float(args.lheight)
    if len(args.lwidth) == 0:
        lwidth = 0.15
    else:
        lwidth = float(args.lwidth)
if lloc == 'top' or lloc == 'bottom':
    if len(args.lheight) == 0:
        lheight = 0.15
    else:
        lheight = float(args.lheight)
    if len(args.lwidth) == 0:
        lwidth = figwidth
    else:
        lwidth = float(args.lwidth)

if len(args.unitpad) == 0:
    if lloc == 'right':
        lpad = 10.0
    if lloc == 'bottom':
        lpad = 0.0
    if lloc == 'left':
        lpad = 0.0
    if lloc == 'top':
        lpad = 5.0
else:
    lpad = float(args.unitpad)

# set base image
if lloc == 'right' or lloc == 'left':
    data = np.empty([256, 1])
    for i in range(0, 256):
        data[i] = i
if lloc == 'top' or lloc == 'bottom':
    data = np.empty([1, 256])
    for i in range(0, 256):
        data[0][i] = i

fig = plt.figure(figsize=(lwidth, lheight))
ax = fig.add_axes([0, 0, 1, 1])
cb = ax.imshow(data, aspect='auto', cmap=args.colormap)

# set colorbar label and styles
if len(args.unitsize) == 0:
    lufs = 12
else:
    lufs = float(args.unitsize)

if args.unit is not None:
    if lloc == 'right' or lloc == 'left':
        ax.set_ylabel(legend_units, rotation=lrotate, labelpad=lpad, fontproperties=font)
        ax.yaxis.set_label_position(lloc)
        ax.yaxis.label.set_fontsize(lufs)
    if lloc == 'top' or lloc == 'bottom':
        ax.set_xlabel(legend_units, rotation=lrotate, labelpad=lpad, fontproperties=font)
        ax.xaxis.set_label_position(lloc)
        ax.xaxis.label.set_fontsize(lufs)

# set colorbar tick styles: font size, family, and ticks direction
if len(args.lticksize) == 0:
    ltfs = lufs - 1
else:
    ltfs = float(args.lticksize)

ax.tick_params(
    direction='out',
    axis='x',  # changes apply to the x-axis
    which='both',  # both major and minor ticks are affected
    top=ltop,  # ticks along the left axis
    bottom=lbottom,  # ticks along the right axis
    labeltop=ltop,  # labels along the left axis
    labelbottom=lbottom)
ax.tick_params(
    direction='out',
    axis='y',  # changes apply to the x-axis
    which='both',  # both major and minor ticks are affected
    left=lleft,  # ticks along the left axis
    right=lright,  # ticks along the right axis
    labelleft=lleft,  # labels along the left axis
    labelright=lright)

for l in ax.yaxis.get_ticklabels():
    l.set_fontproperties(font)
    l.set_fontsize(ltfs)
for l in ax.xaxis.get_ticklabels():
    l.set_fontproperties(font)
    l.set_fontsize(ltfs)

# linear data norm
if args.norm == 'linear':

    # set colorbar major ticks
    if len(args.ld) == 0:
        ld = nice((plot_max_value - plot_min_value) / 5.0)
    else:
        ld = float(args.ld)

    if len(args.ltickbeg) == 0:
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
    if len(args.ltickend) == 0:
        ltickend = plot_max_value
    else:
        ltickend = float(args.ltickend)

    # scalar
    maxtick = max(abs(ltickbeg), abs(ltickend))
    if maxtick >= 1.0e4 or maxtick <= 1.0e-3:
        scalar = int(np.floor(np.log10(maxtick)))
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
        ax.yaxis.set_ticks((ticks - pminv) / (pmaxv - pminv) * 255)
    else:
        ax.xaxis.set_ticks((ticks - pminv) / (pmaxv - pminv) * 255)

    # add power
    # if cscale != 1.0:

    #     last_tick = (ticks[-1] - pminv) / (pmaxv - pminv) * 255

    #     if lloc == 'left':
    #         p1 = -1.1
    #         p2 = max(1.01, 1.075*last_tick + 0.75 * ltfs * 0.01388888889 / lheight)
    #         ha = 'right'
    #         va = 'bottom'
    #     if lloc == 'right':
    #         p1 = 1.1
    #         p2 = max(1.01, 1.075*last_tick + 0.75 * ltfs * 0.01388888889 / lheight)
    #         ha = 'left'
    #         va = 'bottom'
    #     if lloc == 'top':
    #         p1 = max(1.01, 1.05*last_tick + 0.75 * ltfs * 0.01388888889 / lwidth)
    #         p2 = -1.6
    #         ha = 'left'
    #         va = 'center'
    #     if lloc == 'bottom':
    #         p1 = max(1.01, 1.05*last_tick + 0.75 * ltfs * 0.01388888889 / lwidth)
    #         p2 = 1.6
    #         ha = 'left'
    #         va = 'center'
    #     ax.text(p1,
    #             p2,
    #             r'$\mathregular{\times 10^{%i}}$' % scalar,
    #             fontproperties=font,
    #             size=ltfs,
    #             ha=ha,
    #             va=va)

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
        ax.set_yticklabels(tick_labels)
    else:
        if cscale != 1:
            if lloc == 'top':
                tick_labels[-1] = r'$\mathregular{\times 10^{%i}}$' % scalar + '\n' + tick_labels[-1]
            else:
                tick_labels[-1] = tick_labels[-1] + '\n' + r'$\mathregular{\times 10^{%i}}$' % scalar
        ax.set_xticklabels(tick_labels)

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
            ax.yaxis.set_ticks((mticks - pminv) / (pmaxv - pminv) * 255, minor=True)
        else:
            ax.xaxis.set_ticks((mticks - pminv) / (pmaxv - pminv) * 255, minor=True)

# log data norm
if args.norm == 'log':

    plot_min_value = np.log10(plot_min_value)
    plot_max_value = np.log10(plot_max_value)

    # set colorbar major ticks
    if len(args.ltickbeg) == 0:
        ltickbeg = np.floor(plot_min_value)
    else:
        ltickbeg = float(args.ltickbeg)
    if len(args.ltickend) == 0:
        ltickend = np.ceil(plot_max_value) + 1
    else:
        ltickend = float(args.ltickend)
    if len(args.ld) == 0:
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
        ax.yaxis.set_ticks((ticks - pminv) / (pmaxv - pminv) * 255)
    else:
        ax.xaxis.set_ticks((ticks - pminv) / (pmaxv - pminv) * 255)

    # set tick labels on colorbar
    tick_labels = ['' for i in range(0, len(ticks))]
    for i in range(0, len(ticks)):
        tick_labels[i] = '$10^{%i}$' % (ticks[i])
    if lloc in ['left', 'right']:
        ax.yaxis.set_ticklabels(tick_labels)
    else:
        ax.xaxis.set_ticklabels(tick_labels)

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
            mticks = np.append(mticks, np.log10(np.linspace(10**pticks[i], 10**pticks[i + 1],
                                                            args.lmtick + 2)))
        mticks = [i for i in mticks if (i not in pticks)]
        mticks = np.asarray(
            [i for i in mticks if i >= tbeg - 1.0e-10 * abs(tbeg) and i <= tend + 1.0e-10 * abs(tend)])
        # set minor ticks
        if lloc == 'left' or lloc == 'right':
            ax.yaxis.set_ticks((mticks - pminv) / (pmaxv - pminv) * 255, minor=True)
        else:
            ax.xaxis.set_ticks((mticks - pminv) / (pmaxv - pminv) * 255, minor=True)

# colorbar reverse
if lloc in ['left', 'right']:
    if args.lreverse is None or args.lreverse == 1:
        ax.invert_yaxis()
else:
    if args.lreverse == 1:
        ax.invert_xaxis()

## output
from module_output import *
output(args)

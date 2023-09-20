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

# ignore warnings
import warnings
warnings.filterwarnings("ignore", module="matplotlib")

# tag
program = 'colorbar'
print()

# read arguments
parser = argparse.ArgumentParser(description='''
                                purpose:
                                    Plot a horizontal or vertical colorbar.
                                ''',
                                 formatter_class=RawTextHelpFormatter)
parser = getarg(parser, program)
args = parser.parse_args()

# set font
from module_font import *
font, fontbold = set_font(args)

# default plot range
if args.cmin is None:
    plot_min_value = 0
else:
    plot_min_value = float(args.cmin)

if args.cmax is None:
    plot_max_value = 255
else:
    plot_max_value = float(args.cmax)

if args.lloc is None:
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

if args.unitpad is None:
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
cb = ax.imshow(data, aspect='auto')

from module_clip import *
cmin, cmax = np.min(data), np.max(data)

from module_colormap import set_colormap, set_colormap_alpha
colormap = set_colormap(args)
colormap = set_colormap_alpha(args, colormap, cmin, cmax)
cb.set_cmap(colormap)

# set colorbar label and styles
if args.unitsize is None:
    lufs = 12
else:
    lufs = float(args.unitsize)

if args.unit is not None:
    if lloc == 'right' or lloc == 'left':
        ax.set_ylabel(args.unit, rotation=lrotate, labelpad=lpad, fontproperties=font)
        ax.yaxis.set_label_position(lloc)
        ax.yaxis.label.set_fontsize(lufs)
    if lloc == 'top' or lloc == 'bottom':
        ax.set_xlabel(args.unit, rotation=lrotate, labelpad=lpad, fontproperties=font)
        ax.xaxis.set_label_position(lloc)
        ax.xaxis.label.set_fontsize(lufs)

# set colorbar tick styles: font size, family, and ticks direction
if args.lticksize is None:
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
    ax.invert_yaxis()
    if args.lreverse is None or args.lreverse:
        ax.invert_yaxis()
else:
    if args.lreverse:
        ax.invert_xaxis()

## output
from module_io import *
output(args)

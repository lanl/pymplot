#!/usr/bin/python3

# dependencies
from module_io import *
from module_message import fatal
from module_annotation import *
from module_title import *
from module_gridline import *
from module_tick import *
from module_frame import *
from module_font import *
from module_clip import *
from module_size import *
from module_range import *
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
import argparse
from module_getarg import getarg
from module_utility import *
from module_colormap import set_colormap
from argparse import RawTextHelpFormatter

import warnings
warnings.filterwarnings("ignore", module="matplotlib")

program = 'contour'


def run(args, ax=None):

    # input
    data, n1, n2, dmin, dmax = read_array(args, which='fore', dim=2)
    backdata, _, _, backdmin, backdmax = read_array(args, which='back', dim=2)
    mask, _, _, _, _ = read_array(args, which='mask', dim=2)
    if mask is not None:
        data = data * mask
        if backdata is not None:
            backdata = backdata * mask

    # dimension information
    d1 = float(args.d1)
    d2 = float(args.d2)

    # limit of axis
    sp1beg, sp1end, x1beg, x1end, n1beg, n1end = set_range(args.o1, n1, d1, args.x1beg, args.x1end)
    sp2beg, sp2end, x2beg, x2end, n2beg, n2end = set_range(args.o2, n2, d2, args.x2beg, args.x2end)
    data = data[n1beg:n1end, n2beg:n2end]
    if args.background is not None:
        backdata = backdata[n1beg:n1end, n2beg:n2end]

    # figure size
    figheight, figwidth = set_size(args, n1beg, n1end, n2beg, n2end)

    # set clip
    cmin, cmax = set_clip(args, data, 'fore', dmin, dmax)
    if args.norm == 'log':
        if cmin > np.floor(cmax) or cmax < np.ceil(cmin):
            fatal("log colorbar requires values spanning more than one decade; use norm='linear' or adjust cmin/cmax")

    # set font
    font, fontbold = set_font(args)

    # initiate plot
    _subplot = ax is not None
    if not _subplot:
        fig = plt.figure(figsize=(figwidth, figheight))
        ax = fig.add_axes([0, 0, 1, 1])
    else:
        fig = ax.figure
        pos = ax.get_position()
        figwidth  = pos.width  * fig.get_figwidth()
        figheight = pos.height * fig.get_figheight()
    plt.sca(ax)

    # set frame
    set_frame(args)

    # if automatically determine contour levels
    mcontour = int(args.mcontour) + 1
    if mcontour <= 0:
        fatal('mcontour must be >= 0')

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

    # plot contours
    x = np.linspace(0, figwidth, n2end - n2beg)
    y = np.linspace(0, figheight, n1end - n1beg)
    xx, yy = np.meshgrid(x, y)

    # show filled contours if necessary
    if args.contourfill and not args.overlay and args.background is None:
        # set colormap
        colormap = set_colormap(args)
        # plot contour face
        cf = plt.contourf(xx,
                          yy,
                          data,
                          levels[0:size(levels) - 1],
                          cmap=colormap,
                          extend=args.contourextend,
                          antialiased=True,
                          edgecolor='face',
                          linewidth=0.025)

    # show ordinary contours by default
    cs = plt.contour(xx,
                     yy,
                     data,
                     levels[0:size(levels) - 1:args.contourevery],
                     colors=lc,
                     linewidths=lw,
                     linestyles=ls,
                     antialiased=True)

    # this must be placed here, before clabel!
    ax.invert_yaxis()

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
                    clabels[i] = (
                        '%f' %
                        (lvl[i] / cscale)).rstrip('0').rstrip('.') + r'$\mathregular{\times 10^{%i}}$' % scalar
                else:
                    clabels[i] = ('%f' % (lvl[i])).rstrip('0').rstrip('.')

        if args.norm == 'log':
            for i in range(0, size(lvl)):
                clabels[i] = r'$\mathregular{10^{%i}}$' % (lvl[i])

        fmt = {}
        for l, s in zip(cs.levels[::mcontour], clabels):
            fmt[l] = s

        # place contour labels
        clabels = ax.clabel(cs, cs.levels[::mcontour], fmt=fmt, fontsize=clabelsize)
        for txt in clabels:
            txt.set_fontproperties(font)
            txt.set_fontsize(clabelsize)
            txt.set_color(args.clabelcolor)
            if args.clabelbackcolor is not None:
                txt.set_backgroundcolor(args.clabelbackcolor)

    # show original image if necessary
    if args.overlay and args.background is None:

        # begin plot
        im = ax.imshow(data)

        # set colormap
        colormap = set_colormap(args)
        im.set_cmap(colormap)

        # set clip
        im.set_clim(cmin, cmax)

        # set interpolation
        im.set_interpolation(args.interp)

        # set figure sizes based
        im.set_extent([0, figwidth, figheight, 0])

    if not args.overlay and args.background is not None:

        # beg plot
        im = ax.imshow(backdata)

        # set colormap
        colormap = set_colormap(args, 'background')
        im.set_cmap(colormap)

        # set clip
        backcmin, backcmax = set_clip(args, backdata, 'back', backdmin, backdmax)
        im.set_clim(backcmin, backcmax)

        # set interpolation
        im.set_interpolation(args.interp)

        # set figure sizes based
        im.set_extent([0, figwidth, figheight, 0])

    else:
        backcmin = None
        backcmax = None

    # set tick
    set_tick(args, font, x1beg, x1end, n1beg, n1end, d1, figheight, x2beg, x2end, n2beg, n2end, d2, figwidth)

    # set grid line
    set_gridline(args)

    # set title
    set_title(args, fontbold)

    # set annotation
    set_annotation(args, font, x1beg, n1end - n1beg, d1, figheight, x2beg, n2end - n2beg, d2, figwidth)

    # set colorbar
    if args.legend or args.backlegend:

        from module_colorbar import set_colorbar, set_colorbar_contour

        if args.overlay and args.background is None:
            set_colorbar(args, im, font, cmin, cmax, figheight, figwidth, fig)

        if args.backlegend:
            set_colorbar(args, im, font, backcmin, backcmax, figheight, figwidth, fig)

        if args.contourfill and not args.overlay and args.background is None:
            set_colorbar_contour(args, cs, cf, font, cmin, cmax, figheight, figwidth, fig)

    # axis invert
    if args.reverse1:
        ax.invert_yaxis()
    if args.reverse2:
        ax.invert_xaxis()

    # output
    if not _subplot:
        output(args)

    return fig


if __name__ == '__main__':
    print()
    parser = argparse.ArgumentParser(description='''
                                    purpose:
                                        Plot a 2D array as contours.
                                    ''',
                                     formatter_class=RawTextHelpFormatter)
    parser = getarg(parser, program)
    args = parser.parse_args()
    run(args)

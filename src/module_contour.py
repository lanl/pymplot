'''
    Module:
        Set contours
'''
import matplotlib.pyplot as plt
from module_utility import *
import numpy as np
from module_clip import *


def add_contour(args, figwidth, figheight, n1beg, n1end, n2beg, n2end, ax, data, backdata, levels, lc, lw, ls,
                clabelsize, mcontour, font, cmin, cmax, backcmin, backcmax):

    # plot contours
    x = np.linspace(0, figwidth, n2end - n2beg)
    y = np.linspace(0, figheight, n1end - n1beg)
    xx, yy = np.meshgrid(x, y)

    # show filled contours if necessary
    if args.contourfill and not args.overlay and args.background is None:
        # set colormap
        from module_colormap import set_colormap
        colormap = set_colormap(args)
        # plot contour face
        cf = ax.contourf(xx,
                         yy,
                         data,
                         levels[0:len(levels) - 1],
                         cmap=colormap,
                         extend='both',
                         antialiased=True)
        for l in cf.collections:
            l.set_edgecolor('face')
            l.set_linewidth(0.025)

    # show ordinary contours by default
    cs = ax.contour(xx,
                    yy,
                    data,
                    levels[0:len(levels) - 1],
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
        clabels = ['' for i in range(0, len(lvl))]
        if args.norm == 'linear':
            for i in range(0, len(lvl)):
                if lvl[i] != 0 and (abs(lvl[i]) < 1.0e-3 or abs(lvl[i]) > 1.0e3):
                    scalar = int(np.floor(np.log10(abs(lvl[i]))))
                    cscale = 10**scalar #pow(10, scalar)
                    clabels[i] = ('%f' % (lvl[i] / cscale)).rstrip('0').rstrip('.') + r'$\mathregular{\times 10^{%i}}$' % scalar
                else:
                    clabels[i] = ('%f' % (lvl[i])).rstrip('0').rstrip('.')

        if args.norm == 'log':
            for i in range(0, len(lvl)):
                clabels[i] = r'$\mathregular{10^{%i}}$' % (lvl[i])

        fmt = {}
        for l, s in zip(cs.levels[::mcontour], clabels):
            fmt[l] = s

        # place contour labels
        clabels = ax.clabel(cs, cs.levels[::mcontour], fmt=fmt,
                            fontsize=clabelsize)
        for txt in clabels:
            txt.set_fontproperties(font)
            txt.set_fontsize(clabelsize)
            txt.set_color(args.clabelcolor)
            if args.clabelbackcolor is not None:
                txt.set_backgroundcolor(args.clabelbackcolor)

    ## show original image if necessary
    if args.overlay and args.background is None:

        # begin plot
        im = ax.imshow(data)

        # set colormap
        from module_colormap import set_colormap_alpha
        colormap = set_colormap_alpha(args, args.colormap, cmin, cmax, 'foreground')
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

        # set clip
        im.set_clim(backcmin, backcmax)

        # set colormap
        from module_colormap import set_colormap_alpha
        colormap = set_colormap_alpha(args, args.backcolormap, backcmin, backcmax, 'background')
        im.set_cmap(colormap)

        # set interpolation
        im.set_interpolation(args.interp)

        # set figure sizes based
        im.set_extent([0, figwidth, figheight, 0])

    # return plots
    if args.contourfill and not args.overlay and args.background is None:
        return cf

    if args.overlay or args.background is not None:
        return im

#!/usr/bin/python3

## dependencies
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
from module_getarg import getarg
from argparse import RawTextHelpFormatter
from module_io import *
from module_message import fatal, msg_input

import warnings
warnings.filterwarnings("ignore", module="matplotlib")

program = 'wiggle'


def run(args, ax=None):

    ## input data
    infile = args.infile[0].split(',')
    nf = len(infile)

    for i in range(0, nf):
        if not os.path.exists(infile[i]):
            fatal(f'file not found: {infile[i]}')

    is_su = [os.path.splitext(i)[1].lower() == '.su' for i in infile]
    if any(is_su) and not all(is_su):
        fatal('cannot mix SU and raw binary input files')

    su_data = None
    if all(is_su):
        su_input = [read_su(i, args.endian) for i in infile]
        su_n1 = [i[1] for i in su_input]
        su_n2 = [i[2] for i in su_input]
        su_d1 = [i[3] for i in su_input]
        if len(set(su_n2)) != 1:
            fatal('multiple SU files must have the same number of traces')

        requested_d1 = float(args.d1)
        target_d1 = min(su_d1) if requested_d1 == 1.0 else requested_d1
        max_time = min([(su_n1[i] - 1) * su_d1[i] for i in range(0, nf)])
        target_n1 = args.n1 if args.n1 is not None else int(np.floor(max_time / target_d1)) + 1
        if (target_n1 - 1) * target_d1 > max_time:
            fatal('requested n1/d1 exceeds the shortest SU file time range')

        su_data = []
        for data, _, _, d1 in su_input:
            su_data.append(resample_su_data(data, d1, target_n1, target_d1))

        n1 = target_n1
        n2 = su_n2[0]
        d1 = target_d1
    else:
        fsize = os.path.getsize(infile[0])
        datatype = args.datatype
        if datatype == 'double':
            fsize = fsize / 8
        if datatype == 'float':
            fsize = fsize / 4
        if datatype == 'int':
            fsize = fsize / 2

        n1 = args.n1
        if args.n2 is None:
            n2 = int(fsize * 1.0 / n1)
        else:
            n2 = args.n2

        d1 = float(args.d1)
    d2 = float(args.d2)

    # if wiggle locations are read from a file
    if args.wiggleloc is not None:

        # from ASCII file
        wloc = np.loadtxt(args.wiggleloc, ndmin=2)
        [i for i in wloc if i]
        wloc = np.array(wloc)
        wloc_shape = wloc.shape
        if args.along == 1 and n2 != wloc_shape[0]:
            fatal(f'wiggleloc has {wloc_shape[0]} entries but n2={n2}')
        if args.along == 2 and n1 != wloc_shape[0]:
            fatal(f'wiggleloc has {wloc_shape[0]} entries but n1={n1}')

        # n, o, d
        if args.along == 1:
            n2 = wloc_shape[0]
            f2 = wloc[0][0]
            d2 = (wloc.max() - wloc.min()) / (n2 - 1.0)

            n1 = n1
            f1 = float(args.o1)
            d1 = d1
        else:
            n1 = wloc_shape[0]
            f1 = wloc[0][0]
            d1 = (wloc.max() - wloc.min()) / (n1 - 1.0)

            n2 = n2
            f2 = float(args.o2)
            d2 = d2

    else:
        f1 = float(args.o1)
        f2 = float(args.o2)

        wloc = np.zeros((2, 1))
        wloc[0, 0] = f1
        wloc[1, 0] = f2

    ## limit of axis
    from module_range import set_range
    sp1beg, sp1end, x1beg, x1end, n1beg, n1end = set_range(f1, n1, d1, args.x1beg, args.x1end)
    sp2beg, sp2end, x2beg, x2end, n2beg, n2end = set_range(f2, n2, d2, args.x2beg, args.x2end)

    # data type
    from module_datatype import set_datatype
    dt = set_datatype(args)

    if args.scaling is None:
        scaling = [1.0 for i in range(0, nf)]
    else:
        scaling = args.scaling
        if isinstance(scaling, (list, tuple, np.ndarray)):
            scaling = ','.join([str(i) for i in scaling])
        try:
            scaling = [float(i) for i in scaling.split(',')]
        except ValueError:
            fatal('scaling must be a comma-separated list of numbers')
        if len(scaling) == 1 and nf > 1:
            scaling = [scaling[0] for i in range(0, nf)]
        if len(scaling) > nf:
            fatal(f'scaling has {len(scaling)} values but only {nf} input file(s)')
        if len(scaling) < nf:
            scaling.extend([1.0 for i in range(len(scaling), nf)])

    adata = np.empty([nf, n1end - n1beg, n2end - n2beg])
    dmins = []
    dmaxs = []
    for i in range(0, nf):

        if su_data is not None:
            data = su_data[i]
        else:
            # read binary data
            data = fromfile(infile[i], dtype=dt, count=n1 * n2)
            if not args.transpose:
                data = data.reshape((n2, n1))
                data = data.transpose()
            else:
                data = data.reshape(n1, n2)

        data = data * scaling[i]

        # print value range of input files
        if isnan(sum(data)) == True:
            udata = data[~isnan(data)]
            if udata.shape == (0, ):
                fatal(f'dataset {infile[i]} is entirely NaN')
            else:
                dmin = udata.min()
                dmax = udata.max()
        else:
            dmin = data.min()
            dmax = data.max()
        msg_input(infile[i], data.shape, dmin, dmax)
        dmins.append(dmin)
        dmaxs.append(dmax)

        # crop data
        data = data[n1beg:n1end, n2beg:n2end]

        # assign to whole data
        adata[i, :, :] = data

    dmin = min(dmins)
    dmax = max(dmaxs)

    # read background data
    if args.background is not None:

        # check if background file exists
        if not os.path.exists(args.background):
            fatal(f'file not found: {args.background}')

        # input
        backdata = np.empty([n1, n2])
        backdata = fromfile(args.background, dtype=dt, count=n1 * n2)

        # transpose
        if not args.transpose:
            backdata = backdata.reshape((n2, n1))
            backdata = backdata.transpose()
        else:
            backdata = backdata.reshape(n1, n2)

        # crop
        backdata = backdata[n1beg:n1end, n2beg:n2end]
        backdmin = backdata.min()
        backdmax = backdata.max()

    ## set clip
    from module_clip import set_clip
    cmin, cmax = set_clip(args, adata, 'fore', dmin, dmax)
    if args.background is not None:
        backcmin, backcmax = set_clip(args, backdata, 'back', backdmin, backdmax)

    ## figure size
    from module_size import set_size
    figheight, figwidth = set_size(args, n1beg, n1end, n2beg, n2end)

    ## begin plot
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

    ## set frame
    from module_frame import set_frame
    set_frame(args)

    ## plot image
    # show image if necessary
    if nf == 1 and args.overlay and args.background is None:

        # show data by imshow
        im = ax.imshow(adata[0, :, :])

        # set colormap
        from module_colormap import set_colormap
        colormap = set_colormap(args)
        im.set_cmap(colormap)

        # set clip
        im.set_clim(cmin, cmax)

        # set interpolation
        im.set_interpolation(args.interp)

    # plot background image if necessary
    if not args.overlay and args.background is not None:

        # beg plot
        im = ax.imshow(backdata)

        # set colormap
        from module_colormap import set_colormap
        colormap = set_colormap(args, 'background')
        im.set_cmap(colormap)

        # set clip
        im.set_clim(backcmin, backcmax)

        # set interpolation
        im.set_interpolation(args.interp)

    ## set font
    from module_font import set_font
    font, fontbold = set_font(args)

    ## set tick
    from module_tick import set_tick
    set_tick(args,
             font,
             x1beg,
             x1end,
             n1beg,
             n1end,
             d1,
             figheight,
             x2beg,
             x2end,
             n2beg,
             n2end,
             d2,
             figwidth,
             extend=True)

    ## set grid line
    from module_gridline import set_gridline
    set_gridline(args)

    ## set title
    from module_title import set_title
    set_title(args, fontbold)

    ## set annotation
    from module_annotation import set_annotation, set_default
    set_annotation(args, font, x1beg, n1end - n1beg, d1, figheight, x2beg, n2end - n2beg, d2, figwidth)

    ## plot wiggles
    ax = plt.gca()

    # scaling factors
    if args.wiggleloc is not None and args.along == 2:
        scale1 = wloc.max() - wloc.min() + d1
        scale1 = scale1 / figheight
    else:
        scale1 = (n1end - n1beg) * d1 / figheight
    if scale1 == 0:
        scale1 = 1.0

    if args.wiggleloc is not None and args.along == 1:
        scale2 = wloc.max() - wloc.min() + d2
        scale2 = scale2 / figwidth
    else:
        scale2 = (n2end - n2beg) * d2 / figwidth
    if scale2 == 0:
        scale2 = 1.0

    # setup line colors
    from itertools import cycle
    defaultcolor = cycle(['blue', 'red', 'green', 'yellow', 'black', 'cyan', 'magenta'])
    color = args.wigglecolor[0].split(',')
    if args.tracecolor is not None:
        tracecolor = args.tracecolor[0].split(',')

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

    linewidth = set_default(args.wigglewidth, ',', nf, 1.0, 'float')
    linestyle = set_default(args.wigglestyle, ',', nf, '-')

    # plot labels if necessary
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

    # start wiggle plot
    from scipy.interpolate import InterpolatedUnivariateSpline

    if args.along == 1:

        if x1beg == x1end:
            fatal('wiggle axis 1 contains only one sample — check n1 or the crop range')

        traceinterval = int(args.every) * abs(d2) / scale2
        adata = adata / (cmax - cmin) * traceinterval
        traces = np.arange(0, n2end - n2beg, int(args.every))
        y = (np.arange(0, n1end - n1beg) * d1) / scale1

        if args.interp1 is not None:
            yy = np.linspace(y.min(), y.max(), int((n1end - n1beg) * float(args.interp1) + 1))
        else:
            yy = y

        # iterate through all datasets
        for j in range(0, nf):

            # select data
            data = adata[j, :, :]

            for i in traces:

                if args.wiggleloc is not None:
                    offset = (wloc[i][0] - wloc[0][0] + 0.5 * d2) / scale2
                else:
                    offset = (i * d2 + 0.5 * d2) / scale2

                if args.tracecolor is not None:
                    cc = tracecolor[i]
                else:
                    cc = color[j]

                # plot data
                x = data[:, i] + offset

                if args.interp1 is not None:
                    spl = InterpolatedUnivariateSpline(y, x, k=3)
                    xx = spl(yy)
                else:
                    xx = x

                xx[where(xx >= offset + traceinterval)] = offset + traceinterval
                xx[where(xx <= offset - traceinterval)] = offset - traceinterval

                # wiggles
                if i != traces[-1]:
                    plt.plot(xx,
                             yy,
                             color=cc,
                             linewidth=linewidth[j],
                             linestyle=linestyle[j],
                             antialiased=True)
                else:
                    plt.plot(xx,
                             yy,
                             color=cc,
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
                                     color=cc,
                                     edgecolor='none')
                if args.fill == -1:
                    ax.fill_betweenx(yy,
                                     offset,
                                     xx,
                                     interpolate=True,
                                     antialiased=True,
                                     lw=0,
                                     where=(xx < offset),
                                     color=cc,
                                     edgecolor='none')

                if i == traces[0]:
                    wmin = min(wmin, xx.min())
                if i == traces[-1]:
                    wmax = max(wmax, xx.max())

        # to ensure correct image show range
        omin = (traces[0] * d2) / scale2
        omax = (traces[-1] * d2 + d2) / scale2

        # invert y axis to make zero at the top
        if not args.reverse1:
            ax.invert_yaxis()

    if args.along == 2:

        if x2beg == x2end:
            fatal('wiggle axis 2 contains only one sample — check n2 or the crop range')

        traceinterval = int(args.every) * abs(d1) / scale1
        adata = adata / (cmax - cmin) * traceinterval
        traces = np.arange(0, n1end - n1beg, int(args.every))
        y = (np.arange(0, n2end - n2beg) * d2) / scale2

        if args.interp2 is not None:
            yy = np.linspace(y.min(), y.max(), int((n2end - n2beg) * float(args.interp2) + 1))
        else:
            yy = y

        # iterate through all datasets
        for j in range(0, nf):

            # select data
            data = adata[j, :, :]

            for i in traces:

                if args.wiggleloc is not None:
                    offset = (wloc[i][0] - wloc[0][0] + 0.5 * d1) / scale1
                else:
                    offset = (i * d1 + 0.5 * d1) / scale1

                # plot data
                x = data[i, :] + offset

                if args.tracecolor is not None:
                    cc = tracecolor[i]
                else:
                    cc = color[j]

                if args.interp2 is not None:
                    spl = InterpolatedUnivariateSpline(y, x, k=3)
                    xx = spl(yy)
                else:
                    xx = x

                xx[where(xx >= offset + traceinterval)] = offset + traceinterval
                xx[where(xx <= offset - traceinterval)] = offset - traceinterval

                # wiggles
                if i != traces[-1]:
                    plt.plot(yy,
                             xx,
                             color=cc,
                             linewidth=linewidth[j],
                             linestyle=linestyle[j],
                             antialiased=True)
                else:
                    plt.plot(yy,
                             xx,
                             color=cc,
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
                                    facecolor=cc,
                                    edgecolor='none')
                if args.fill == -1:
                    ax.fill_between(yy,
                                    offset,
                                    xx,
                                    interpolate=True,
                                    antialiased=True,
                                    lw=0,
                                    where=(xx < offset),
                                    facecolor=cc,
                                    edgecolor='none')

                if i == traces[0]:
                    wmin = min(wmin, xx.min())
                if i == traces[-1]:
                    wmax = max(wmax, xx.max())

        # to ensure correct image show range
        omin = (traces[0] * d1) / scale1
        omax = (traces[-1] * d1 + d1) / scale1

    # add plot labels
    if args.plotlabel is not None:
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
    if (nf == 1 and args.overlay and args.background is None) or (not args.overlay
                                                                  and args.background is not None):
        if args.along == 1:
            im.set_extent([omin, omax, figheight, 0])
        else:
            im.set_extent([0, figwidth, omax, omin])

    apad = float(args.axispad)

    if args.along == 1:
        wmin = wmin - apad
        wmax = wmax + apad
        if args.wx2beg is not None:
            wmin = (float(args.wx2beg) - wloc[0][0] + 0.5 * d2) / scale2
        if args.wx2end is not None:
            wmax = (float(args.wx2end) - wloc[0][0] + 0.5 * d2) / scale2
        ax.set_xlim([wmin, wmax])
        ax.set_ylim([figheight, 0])
        ax.set_aspect('auto')
    if args.along == 2:
        wmin = wmin - apad
        wmax = wmax + apad
        if args.wx1beg is not None:
            wmin = (float(args.wx1beg) - wloc[0][0] + 0.5 * d1) / scale1
        if args.wx1end is not None:
            wmax = (float(args.wx1end) - wloc[0][0] + 0.5 * d1) / scale1
        ax.set_xlim([0, figwidth])
        if not args.reverse1:
            ax.set_ylim([wmin, wmax])
        else:
            ax.set_ylim([wmax, wmin])
        ax.set_aspect('auto')

    ## axis invert
    if args.reverse1:
        ax.invert_yaxis()
    if args.reverse2:
        ax.invert_xaxis()

    ## set colorbar
    if nf == 1 and args.overlay and args.legend:
        from module_colorbar import set_colorbar
        set_colorbar(args, im, font, cmin, cmax, figheight, figwidth, fig)

    ## output
    if not _subplot:
        output(args)

    return fig


if __name__ == '__main__':
    print()
    parser = argparse.ArgumentParser(description='''
                                    purpose:
                                        Plot a 2D array as wiggles along the 1st or the 2nd dimension.
                                    ''',
                                     formatter_class=RawTextHelpFormatter)
    parser = getarg(parser, program)
    args = parser.parse_args()
    run(args)

'''
    Module:
        Set regular or irregular axis ticks for a plot.
'''
from module_utility import *
import numpy as np
import matplotlib.pyplot as plt

# ticks : contains irregular ticks locations
# tickbeg : regular major ticks begin location
# tickend : regular major ticks end location
# tickd : regular major ticks interval
# mtick : number of minor tick intervals betwen two major ticks
# xbeg : axis begin location
# xend : axis end location
# ns : number of points to plot
# d : interval between two points
# axislen : apparent axis length
def define_tick(ticks, tickbeg, tickend, tickd, mtick, xbeg, xend, ns, d, axislen, format, extend=False):

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

        # some ticks might out of axis range, therefore remove them if strict
        if not extend:
            if d > 0:
                tick = np.asarray([i for i in tick if i >= xbeg and i <= xend])
                minor_tick = np.asarray(
                    [i for i in minor_tick if i >= xbeg and i <= xend and (not i in tick)])
            if d < 0:
                tick = np.asarray([i for i in tick if i <= xbeg and i >= xend])
                minor_tick = np.asarray(
                    [i for i in minor_tick if i <= xbeg and i >= xend and (not i in tick)])

        # linearly scale the ticks to figure canvas
        if ns == 1:
            # if only one sample point, then tick location is 0.5
            tick_location = np.asarray([0.5])
            ntick = 1
        else:
            # if multiple sample points, then scale to apparent axis length
            tick_location = [(i - xbeg + 0.5 * d) / ((ns - 1) * d) * axislen for i in tick]
            minor_tick_location = [(i - xbeg + 0.5 * d) / ((ns - 1) * d) * axislen for i in minor_tick]
            t = tick_location

        # set major tick location and labels, note some major ticks might be out of axis range
        tl = []
        tick_label = []
        for i in range(0, len(tick)):
            if extend or ((not extend) and tick_location[i] >= 0 and tick_location[i] <= axislen + 1.0e-10):
                tl.append(tick_location[i])
                if format == 'sci' or format == 'plain':
                    tick_label.append(('%f' % tick[i]).rstrip('0').rstrip('.'))
                else:
                    tick_label.append((format % tick[i]))
        tick_location = tl

    # irregular ticks
    else:
        
        # get contents from user-specified ticks
        ticks = ticks[0].split(',')
        location = [0 for i in range(0, len(ticks))]
        label = ['' for i in range(0, len(ticks))]
        
        # set tick locations
        for i in range(0, len(ticks)):
            t = ticks[i].split(':')
            location[i] = (float(t[0]) + 0.5 * d) / ((ns - 1) * d) * axislen
            label[i] = t[1]

        # sort according to tick location
        yx = list(zip(location, label))
        yx.sort()
        tick_location = [location for location, label in yx]
        tick_label = [label for location, label in yx]

        # minor ticks
        if mtick != 0:
            mtick = mtick + 1
            minor_tick_location = np.linspace(tick_location[0], tick_location[1], mtick + 1)
            minor_tick_location = minor_tick_location[1:mtick]
            for i in range(1, len(tick_location) - 1):
                t = np.linspace(tick_location[i], tick_location[i + 1], mtick + 1)
                minor_tick_location = np.append(minor_tick_location, t[1:mtick])
        else:
            minor_tick_location = []

    # return major tick location, major tick label and minor tick location
    return tick_location, tick_label, minor_tick_location


def set_tick(args,
             font,
             x1beg,
             x1end,
             n1beg,
             n1end,
             d1,
             axis1len,
             x2beg,
             x2end,
             n2beg,
             n2end,
             d2,
             axis2len,
             extend=False):

    ax = plt.gca()

    label_1_size = float(args.label1size)
    label_2_size = float(args.label2size)

    xlabel = ax.set_xlabel(args.label2, fontsize=label_2_size, labelpad=float(args.label2pad)*72*2)
    ylabel = ax.set_ylabel(args.label1, fontsize=label_1_size, labelpad=float(args.label1pad)*72*2)

    l = ax.yaxis.get_label()
    l.set_fontproperties(font)
    l.set_fontsize(label_1_size)
    l = ax.xaxis.get_label()
    l.set_fontproperties(font)
    l.set_fontsize(label_2_size)

    if args.label2loc is not None:
        ax.xaxis.set_label_position(args.label2loc)
    else:
    	if args.ticktop:
    		ax.xaxis.set_label_position('top')
    	else:
    		ax.xaxis.set_label_position('bottom')
    if args.label1loc is not None:
    	ax.yaxis.set_label_position(args.label1loc)
    else:
    	if args.tickleft:
    		ax.yaxis.set_label_position('left')
    	else:
    		ax.yaxis.set_label_position('right')
    		ylabel.set_rotation(270)

    # ticks on/off
    ax.get_yaxis().set_tick_params(which='both', direction='out')
    ax.get_xaxis().set_tick_params(which='both', direction='out')
    plt.tick_params(
        axis='x',  # changes apply to the x1-axis
        which='both',  # both major and minor ticks are affected
        bottom=args.tickbottom,  # ticks along the bottom axis
        top=args.ticktop,  # ticks along the top axis
        labelbottom=args.tickbottom,  # labels along the bottom axis
        labeltop=args.ticktop)  # labels along the top axis
    plt.tick_params(
        axis='y',  # changes apply to the x2-axis
        which='both',  # both major and minor ticks are affected
        left=args.tickleft,  # ticks along the left axis
        right=args.tickright,  # ticks along the right axis
        labelleft=args.tickleft,  # labels along the left axis
        labelright=args.tickright)  # labels along the right axis

    # if tick font size and family not speciefied, then inherit from axis labels
    if args.tick1size is None:
        tick_1_font_size = label_1_size - 2
    else:
        tick_1_font_size = float(args.tick1size)

    if args.tick2size is None:
        tick_2_font_size = label_2_size - 2
    else:
        tick_2_font_size = float(args.tick2size)

    # axis 1
    tick_1_location, tick_1_label, tick_1_minor = define_tick(args.ticks1, args.tick1beg, args.tick1end,
                                                              args.tick1d, args.mtick1, x1beg, x1end,
                                                              n1end - n1beg + 1, d1, axis1len,
                                                              args.tick1format, extend)
    plt.yticks(tick_1_location, tick_1_label, fontsize=tick_1_font_size, rotation=float(args.tick1rot))
    if not args.tick1label:
        ax.yaxis.set_ticklabels([])

    # axis 2
    tick_2_location, tick_2_label, tick_2_minor = define_tick(args.ticks2, args.tick2beg, args.tick2end,
                                                              args.tick2d, args.mtick2, x2beg, x2end,
                                                              n2end - n2beg + 1, d2, axis2len,
                                                              args.tick2format, extend)
    plt.xticks(tick_2_location, tick_2_label, fontsize=tick_2_font_size, rotation=float(args.tick2rot))
    if not args.tick2label:
        ax.xaxis.set_ticklabels([])

    # major and minor ticks sytle
    ax.tick_params('both', length=float(args.tickmajorlen), width=float(args.tickmajorwid), which='major')

    # minor tick positions
    ax.set_yticks(tick_1_minor, minor=True)
    ax.set_xticks(tick_2_minor, minor=True)

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

    for l in ax.yaxis.get_ticklabels():
        l.set_fontproperties(font)
        l.set_fontsize(tick_1_font_size)
    for l in ax.xaxis.get_ticklabels():
        l.set_fontproperties(font)
        l.set_fontsize(tick_2_font_size)


# make tick labels rigid
def rigid_tick_label(tick_label):

    ndec = 0
    for i in tick_label:
        dec = i.split('.')
        if len(dec) == 2:
            ll = len(dec[1])
            if ll > ndec:
                ndec = ll

    for i in range(0, len(tick_label)):
        dec = tick_label[i].split('.')
        if len(dec) == 2:
            ll = len(dec[1])
            if ll < ndec:
                for k in range(0, ndec - ll):
                    tick_label[i] = tick_label[i] + '0'
        if len(dec) == 1 and ndec != 0:
            tick_label[i] = tick_label[i] + '.'
            for k in range(0, ndec):
                tick_label[i] = tick_label[i] + '0'

    return tick_label

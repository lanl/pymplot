import numpy as np
from matplotlib.collections import LineCollection, PolyCollection
from module_utility import *
import math

ipp = 0.0138889

# projection operation
# axonometric projection
def isometric_projection(x, y, z, angles):
    alpha = angles[0]
    beta = angles[1]
    bx = x * np.cos(alpha * np.pi / 180.0) - y * np.cos(beta * np.pi / 180.0)
    by = z + x * np.sin(alpha * np.pi / 180.0) + y * np.sin(beta * np.pi / 180.0)
    return bx, by


# plot image and axis on deformed rectangular mesh
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


def project_axis(ax, p1x, p1y, p2x, p2y, ticks, tickbeg, tickend, tickd, mtick, xbeg, xend, ns, d, font, tick_major_len,
                 tick_major_width, tick_minor_len, tick_minor_width, axislen, tick_orient, ticklabel_orient, tick_size,
                 label, label_orient, label_size, label_pad, tick_format):

    # regular ticks
    if len(ticks) == 0:

        # major tick interval
        if len(tickd) == 0:
            tick_interval = nice((xend - xbeg) / 5.0)
            if tick_interval == 0:
                tick_interval = 1.0e10
        else:
            tick_interval = float(tickd)

        # tick begin location
        if len(tickbeg) == 0:
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
        if len(tickend) == 0:
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
            minor_tick = np.asarray([i for i in minor_tick if i >= xbeg and i <= xend and (i not in tick)])
        if d < 0:
            tick = np.asarray([i for i in tick if i <= xbeg and i >= xend])
            minor_tick = np.asarray([i for i in minor_tick if i <= xbeg and i >= xend and (i not in tick)])

        # linearly scale the ticks to figure canvas
        if ns == 1:
            # if only one sample point, then tick location is 0.5
            tick_location = np.asarray([0.5])
            ntick = 1
        else:
            # if multiple sample points, then scale to apparent axis length
            tick_location = (tick - xbeg + 0.5 * d) / ((ns - 1) * d) * axislen
            minor_tick_location = (minor_tick - xbeg + 0.5 * d) / ((ns - 1) * d) * axislen
            t = tick_location

        # set major tick location and labels, note some major ticks might be
        # out of axis range
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
            location[i] = (float(t[0]) + 0.5 * d) / ((ns - 1) * d) * axislen
            label[i] = t[1]

        # sort according to tick location
        yx = sorted(zip(location, label))
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
        epx = tx[i] + tx2[i] / np.sqrt(tx2[i]**2 + ty2[i]**2) * ticklen
        epy = ty[i] + ty2[i] / np.sqrt(tx2[i]**2 + ty2[i]**2) * ticklen
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
        epx = tx[i] + txx[i] / np.sqrt(txx[i]**2 + tyy[i]**2) * pad
        epy = ty[i] + tyy[i] / np.sqrt(txx[i]**2 + tyy[i]**2) * pad
        ax.text(epx, epy, tick_label[i], ha=ha, va=va, rotation=ticklabel_angle, fontproperties=font, size=tick_size)

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
    labelx = hx + hx2 / np.sqrt(hx2**2 + hy2**2) * pad
    labely = hy + hy2 / np.sqrt(hx2**2 + hy2**2) * pad
    if len(label) != 0:
        ax.text(labelx, labely, label, ha=ha, va=va, rotation=label_angle, fontproperties=font, size=label_size)

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
        epx = tx[i] + tx2[i] / np.sqrt(tx2[i]**2 + ty2[i]**2) * ticklen
        epy = ty[i] + ty2[i] / np.sqrt(tx2[i]**2 + ty2[i]**2) * ticklen
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
'''
    Module:
        Set colorbar
'''
import matplotlib as mplt
import matplotlib.pyplot as plt
from module_utility import *
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from math import *
from module_projection import *
from module_colormap import *
import numpy as np


def set_colorbar(args, im, font, plot_min_value, plot_max_value, figheight, figwidth, fig):

    if (args.legend or args.backlegend) and plot_min_value != plot_max_value:

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
        cb = fig.colorbar(im, cax=cax, orientation=lorient)

        # set colorbar label and styles
        if args.unitsize is None:
            lufs = min(float(args.label1size), float(args.label2size)) - 1
        else:
            lufs = float(args.unitsize)

        if args.unit is not None:
            cb.set_label(args.unit, rotation=lrotate, labelpad=lpad, fontproperties=font)

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
            top=ltop,  # ticks along the left axis
            bottom=lbottom,  # ticks along the right axis
            labeltop=ltop,  # labels along the left axis
            labelbottom=lbottom)
        cb.ax.tick_params(
            direction='out',
            axis='y',  # changes apply to the x-axis
            which='both',  # both major and minor ticks are affected
            left=lleft,  # ticks along the left axis
            right=lright,  # ticks along the right axis
            labelleft=lleft,  # labels along the left axis
            labelright=lright)

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
                cb.ax.yaxis.set_ticks(ticks)
            else:
                cb.ax.xaxis.set_ticks(ticks)

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
                    cb.ax.yaxis.set_ticks(mticks, minor=True)
                else:
                    cb.ax.xaxis.set_ticks(mticks, minor=True)

        # make colorbar solid color continuous
        cb.solids.set_edgecolor("face")

        # colorbar reverse
        if lloc in ['left', 'right']:
            if args.lreverse:
                cb.ax.invert_yaxis()
        else:
            if args.lreverse:
                cb.ax.invert_xaxis()
                
                
def set_colorbar_contour(args, cs, cf, font, plot_min_value, plot_max_value, figheight, figwidth, fig):

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
        from matplotlib.cm import ScalarMappable
        cax = fig.add_axes([cbx, cby, lwidth / figwidth, lheight / figheight])
        cb = fig.colorbar(cf, cax=cax, orientation=lorient)
        cb.add_lines(cs)

        # set colorbar label and styles
        if args.unitsize is None:
            lufs = min(float(args.label1size), float(args.label2size)) - 1
        else:
            lufs = float(args.unitsize)

        if args.unit is not None:
            cb.set_label(args.unit, rotation=lrotate, labelpad=lpad, fontproperties=font)

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
            top=ltop,  # ticks along the left axis
            bottom=lbottom,  # ticks along the right axis
            labeltop=ltop,  # labels along the left axis
            labelbottom=lbottom)
        cb.ax.tick_params(
            direction='out',
            axis='y',  # changes apply to the x-axis
            which='both',  # both major and minor ticks are affected
            left=lleft,  # ticks along the left axis
            right=lright,  # ticks along the right axis
            labelleft=lleft,  # labels along the left axis
            labelright=lright)

        for l in cb.ax.yaxis.get_ticklabels():
            l.set_fontproperties(font)
            l.set_fontsize(ltfs)
        for l in cb.ax.xaxis.get_ticklabels():
            l.set_fontproperties(font)
            l.set_fontsize(ltfs)

        # colorbar reverse
        if lloc in ['left', 'right']:
            if args.lreverse:
                cb.ax.invert_yaxis()
        else:
            if args.lreverse:
                cb.ax.invert_xaxis()


# # self-created colorbar
# def custom_colorbar(args, im, font, cmin, cmax, figheight, figwidth, fig):

#     lloc = args.lloc

#     # colorbar values
#     cp = 512
#     cinterval = (cmax - cmin) / (cp - 1)
#     temp = np.linspace(cmin, cmax, cp)
#     cval = np.zeros([cp, 1])
#     for i in range(0, cp):
#         cval[i] = temp[i]

#     if lloc in ['left', 'right']:

#         if args.lheight is None:
#             lheight = figheight
#         else:
#             lheight = float(args.lheight)
#             if lheight > figheight:
#                 lheight = figheight

#         if args.lwidth is None:
#             lwidth = 0.2
#         else:
#             lwidth = float(args.lwidth)

#         if lloc == 'right':

#             # colorbar location
#             pad = 0.1
#             if args.lpad is None:
#                 cbpad = 0.0
#             else:
#                 cbpad = float(args.lpad)
#             pad = pad + cbpad
#             p5x = figwidth
#             p3y = 0
#             p8y = lheight
#             cx = [p5x + pad, p5x + pad, p5x + pad + lwidth, p5x + pad + lwidth]
#             dl = (figheight - lheight) / 2.0
#             cy = [p3y + dl, p8y - dl, p8y - dl, p3y + dl]

#     if args.lloc in ['top', 'bottom']:

#         if args.lheight is None:
#             lheight = 0.2
#         else:
#             lheight = float(args.lheight)

#         if args.lwidth is None:
#             lwidth = figwidth
#         else:
#             lwidth = float(args.lwidth)
#             if lwidth > figwidth:
#                 lwidth = figwidth

#         if lloc == 'bottom':

#             # colorbar location
#             if args.axis2loc == 'bottom' or args.axis2loc == 'both' or args.axis3loc == 'bottom' or args.axis3loc == 'both':
#                 pad = max(tick_2_size, tick_3_size) * ipp + float(args.tickmajorlen) * ipp + 0.2 + max(
#                     label_2_size, label_3_size) * ipp + 0.1 * abs(cos(angle1 * np.pi / 180.0))
#             else:
#                 pad = 0.2
#             if args.lpad is None:
#                 cbpad = 0.0
#             else:
#                 cbpad = float(args.lpad)
#             pad = pad + cbpad
#             dl = (figwidth - lwidth) / 2.0
#             cx = [p1x + dl, p5x - dl, p5x - dl, p1x + dl]
#             cy = [p3y - pad, p3y - pad, p3y - pad - lheight, p3y - pad - lheight]

#     # create colorbar
#     # ax = fig.gca()
#     ax = fig.add_axes([1 + pad/figwidth, 0, lwidth/figwidth, lheight/figheight])
#     colormap = set_colormap(args)
#     cb = project_image(ax, cval, cx, cy, colormap, cmin, cmax)
#     ax.tick_params(
#             direction='out',
#             axis='x',  # changes apply to the x-axis
#             which='both',  # both major and minor ticks are affected
#             top=False,  # ticks along the left axis
#             bottom=False,  # ticks along the right axis
#             labeltop=False,  # labels along the left axis
#             labelbottom=False)
#     ax.tick_params(
#             direction='out',
#             axis='y',  # changes apply to the x-axis
#             which='both',  # both major and minor ticks are affected
#             left=False,  # ticks along the left axis
#             right=False,  # ticks along the right axis
#             labelleft=False,  # labels along the left axis
#             labelright=False)

#     # crate colorbar frame
#     ipp = 0.0138889
#     line = [[(cx[0] - ipp, cy[0] - 0.5*ipp),
#              (cx[1] - ipp, cy[1] + 0.5*ipp),
#              (cx[2], cy[2] + 0.5*ipp),
#              (cx[3], cy[3] - 0.5*ipp),
#              (cx[0] - ipp, cy[0] - 0.5*ipp)]]
#     line = LineCollection(line, linewidth=1.0, color='k')
#     ax.add_collection(line)

#     if args.unitsize is None:
#         lufs = min(float(args.label1size), float(args.label2size)) - 1
#     else:
#         lufs = float(args.unitsize)

#     # tick font size
#     if args.lticksize is None:
#         ltfs = lufs - 1
#     else:
#         ltfs = float(args.lticksize)

#     if args.norm == 'linear':

#         # set colorbar major ticks
#         if args.ld is None:
#             ld = nice((cmax - cmin) / 5.0)
#         else:
#             ld = float(args.ld)

#         if args.ltickbeg is None:
#             ltickbeg = nice(cmin, 0.5)
#             base = 0.5
#             nb = 0
#             while nb <= 10 and ltickbeg > cmin + ld:
#                 base = base / 10.0
#                 ltickbeg = nice(cmin, base)
#                 nb = nb + 1
#             if abs(ltickbeg) < abs(cmax) and orderm(ltickbeg) + 2 < orderm(cmax):
#                 ltickbeg = 0.0
#         else:
#             ltickbeg = float(args.ltickbeg)
#         if args.ltickend is None:
#             ltickend = cmax
#         else:
#             ltickend = float(args.ltickend)

#         # scalar
#         maxtick = max(abs(ltickbeg), abs(ltickend))
#         if maxtick >= 1.0e4 or maxtick <= 1.0e-3:
#             scalar = int(floor(log10(maxtick)))
#             cscale = pow(10, scalar)
#         else:
#             cscale = 1.0

#         # set ticks
#         ticks = np.arange(ltickbeg, ltickend + ld, ld)
#         tbeg = max(cmin, ltickbeg)
#         tend = min(cmax, ltickend)

#         # set tick positions on colorbar
#         ticks = np.asarray([i for i in ticks if i >= tbeg - 1.0e-10 * abs(tbeg) and i <= tend + 1.0e-10 * abs(tend)])
#         tick_labels = ['' for i in range(0, len(ticks))]
#         for i in range(0, len(ticks)):
#             tick_labels[i] = ('%f' % (ticks[i] / cscale)).rstrip('0').rstrip('.')

#         # set minor ticks
#         if args.lmtick != 0:
#             # extend tail and head
#             pticks = np.append(ticks, ticks[0] - ld)
#             pticks = np.append(pticks, ticks[-1] + ld)
#             # sort all major ticks
#             pticks = np.sort(pticks)
#             # get pseudo-location of minor ticks
#             nt = len(pticks)
#             mticks = []
#             for i in range(0, nt - 1):
#                 mticks = np.append(mticks, np.linspace(pticks[i], pticks[i + 1], args.lmtick + 2))
#             mticks = [i for i in mticks if (i not in pticks)]
#             mticks = np.asarray(
#                 [i for i in mticks if i >= tbeg - 1.0e-10 * abs(tbeg) and i <= tend + 1.0e-10 * abs(tend)])

#     if args.norm == 'log':

#         # set colorbar major ticks
#         if args.ltickbeg is None:
#             ltickbeg = np.floor(cmin)
#         else:
#             ltickbeg = float(args.ltickbeg)
#         if args.ltickend is None:
#             ltickend = np.ceil(cmax)
#         else:
#             ltickend = float(args.ltickend)
#         if args.ld is None:
#             ld = max(1, round((ltickend - ltickbeg) / 5.0))
#         else:
#             ld = int(args.ld)

#         ticks = np.arange(ltickbeg, ltickend + 1, ld)
#         tbeg = max(cmin, ltickbeg)
#         tend = min(cmax, ltickend)

#         # set tick positions on colorbar
#         ticks = np.asarray([i for i in ticks if i >= tbeg - 1.0e-10 * abs(tbeg) and i <= tend + 1.0e-10 * abs(tend)])
#         tick_labels = ['' for i in range(0, len(ticks))]
#         for i in range(0, len(ticks)):
#             tick_labels[i] = r'$\mathregular{10^{%i}}$' % (ticks[i])

#         # colorbar minor ticks
#         if args.lmtick != 0:
#             # extend tail and head
#             pticks = np.append(ticks, ticks[0] - ld)
#             pticks = np.append(pticks, ticks[-1] + ld)
#             # sort all major ticks
#             pticks = np.sort(pticks)
#             # get pseudo-location of minor ticks
#             nt = len(pticks)
#             mticks = []
#             for i in range(0, nt - 1):
#                 mticks = np.append(mticks, np.log10(np.linspace(10**pticks[i], 10**pticks[i + 1], args.lmtick + 2)))
#             mticks = np.asarray(
#                 [i for i in mticks if i >= tbeg - 1.0e-10 * abs(tbeg) and i <= tend + 1.0e-10 * abs(tend)])

#     # inch per point
#     ipp = 0.0138889

#     # add ticks by drawing lines
#     if lloc == 'right':

#         # add major ticks
#         cbeg = cy[3]
#         cend = cy[2]
#         ticks = [(i - cmin) / (cmax - cmin) * (cend - cbeg) + cbeg for i in ticks]

#         tx = [cx[2] for i in range(0, len(ticks))]
#         ty = [i for i in ticks]
#         last_tick = ty[-1]
#         ticklen = float(args.tickmajorlen) * ipp * 2
#         tx2 = [i + ticklen for i in tx]
#         ty2 = [i for i in ticks]
#         ttx = [i + ticklen + 0.05 for i in tx]
#         tty = [i for i in ticks]

#         majortick = []
#         e1 = list(zip(tx, ty))
#         e2 = list(zip(tx2, ty2))
#         for i in range(0, len(ticks)):
#             majortick.append([e1[i], e2[i]])

#         tick_major_width = 1
#         majortick = LineCollection(majortick, linewidths=tick_major_width, colors='k')
#         ax.add_collection(majortick)
#         ax.autoscale()
#         ax.axis('off')

#         # add tick labels
#         for i in range(0, len(ticks)):
#             ax.text(ttx[i], tty[i], tick_labels[i], fontproperties=font, size=ltfs, ha='left', va='center')

#         # add minor ticks
#         if args.lmtick != 0:
#             mticks = [(i - cmin) / (cmax - cmin) * (cend - cbeg) + cbeg for i in mticks]

#             tx = [cx[2] for i in range(0, len(mticks))]
#             ty = [i for i in mticks]
#             ticklen = float(args.tickmajorlen) * ipp
#             tx2 = [i + ticklen for i in tx]
#             ty2 = [i for i in mticks]
#             ttx = [i + ticklen + 0.05 for i in tx]
#             tty = [i for i in mticks]

#             minortick = []
#             e1 = list(zip(tx, ty))
#             e2 = list(zip(tx2, ty2))
#             for i in range(0, len(mticks)):
#                 minortick.append([e1[i], e2[i]])

#             tick_minor_width = 0.75 * tick_major_width
#             minortick = LineCollection(minortick, linewidths=tick_minor_width, colors='k')
#             ax.add_collection(minortick)

#         # add power
#         if args.norm == 'linear' and cscale != 1.0:
#             p1 = cx[2] + 0.01
#             p2 = max(cend + 0.01, last_tick + 0.75 * ltfs * ipp)
#             ha = 'left'
#             va = 'bottom'
#             ct = ax.text(p1,
#                          p2,
#                          r'$\mathregular{\times 10^{%i}}$' % scalar,
#                          size=ltfs,
#                          fontproperties=font,
#                          ha=ha,
#                          va=va)
#             ct.size_size(ltfs)

#         # set unit
#         if args.unit is not None:
#             if args.unitpad is None:
#                 upad = 0.05
#             else:
#                 upad = float(args.unitpad)
#             if args.norm == 'linear':
#                 maxlen = max([len(i) for i in tick_labels]) * 0.75
#             if args.norm == 'log':
#                 maxlen = 3.5
#             ux = cx[2] + ticklen * 2.0 + 0.025 + maxlen * ltfs * ipp + upad
#             uy = 0.5 * (cbeg + cend)
#             ct = ax.text(ux, uy, args.unit, size=lufs, fontproperties=font, rotation=270, ha='left', va='center')
#             ct.set_size(lufs)

#     # add ticks by drawing lines
#     if lloc == 'bottom':

#         # add major ticks
#         cbeg = cx[0]
#         cend = cx[1]
#         ticks = [(i - cmin) / (cmax - cmin) * (cend - cbeg) + cbeg for i in ticks]

#         ty = [cy[3] for i in range(0, len(ticks))]
#         tx = [i for i in ticks]
#         last_tick = tx[-1]
#         ticklen = float(args.tickmajorlen) * ipp
#         ty2 = [i - ticklen for i in ty]
#         tx2 = [i for i in ticks]
#         tty = [i - ticklen - 0.05 for i in ty]
#         ttx = [i for i in ticks]

#         majortick = []
#         e1 = list(zip(tx, ty))
#         e2 = list(zip(tx2, ty2))
#         for i in range(0, len(ticks)):
#             majortick.append([e1[i], e2[i]])

#         tick_major_width = 1
#         majortick = LineCollection(majortick, linewidths=tick_major_width, colors='k')
#         ax.add_collection(majortick)

#         # add tick labels
#         for i in range(0, len(ticks)):
#             ct = ax.text(ttx[i], tty[i], tick_labels[i], fontproperties=font, size=ltfs, ha='center', va='top')
#             ct.set_size(ltfs)

#         # add minor ticks
#         if args.lmtick != 0:
#             mticks = [(i - cmin) / (cmax - cmin) * (cend - cbeg) + cbeg for i in mticks]

#             ty = [cy[3] for i in range(0, len(mticks))]
#             tx = [i for i in mticks]
#             ticklen = float(args.tickmajorlen) * ipp * 0.5
#             ty2 = [i - ticklen for i in ty]
#             tx2 = [i for i in mticks]
#             tty = [i - ticklen - 0.05 for i in ty]
#             ttx = [i for i in mticks]

#             minortick = []
#             e1 = list(zip(tx, ty))
#             e2 = list(zip(tx2, ty2))
#             for i in range(0, len(mticks)):
#                 minortick.append([e1[i], e2[i]])

#             tick_minor_width = tick_major_width * 0.75
#             minortick = LineCollection(minortick, linewidths=tick_minor_width, colors='k')
#             ax.add_collection(minortick)

#         # add power
#         if args.norm == 'linear' and cscale != 1.0:
#             p1 = cx[2] + 0.025
#             p2 = cy[3]
#             ha = 'left'
#             va = 'center'
#             ct = ax.text(p1,
#                          p2,
#                          r'$\mathregular{\times 10^{%i}}$' % scalar,
#                          size=ltfs,
#                          fontproperties=font,
#                          ha=ha,
#                          va=va)
#             ct.set_size(ltfs)

#         # set unit
#         if args.unit is not None:
#             if args.unitpad is None:
#                 upad = 0.05
#             else:
#                 upad = float(args.unitpad)
#             if args.norm == 'linear':
#                 maxlen = 1.50
#             if args.norm == 'log':
#                 maxlen = 1.75
#             ux = 0.5 * (cbeg + cend)
#             uy = cy[3] - ticklen - 0.025 - maxlen * ltfs * ipp - upad
#             ct = ax.text(ux, uy, args.unit, size=lufs, fontproperties=font, ha='center', va='top')
#             ct.set_size(lufs)

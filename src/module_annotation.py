'''
    Module:
        Add various types of annotations to a plot. 
'''
from pylab import *
import numpy as np
from matplotlib.patches import *
import matplotlib.pyplot as plt


# set default value
def set_default(arguments, separator, length, default, conversion='none', consistent=False):
    if arguments is not None:
        vals = arguments[0].split(separator)
        if len(vals) < length:
            l = len(vals)
            if consistent:
                # consistent then use the last value of vals
                avals = [vals[-1] for i in range(l, length)]
            else:
                # inconsistent then use given default value
                avals = [default for i in range(l, length)]
            # extend given values
            vals.extend(avals)
    else:
        vals = [default for i in range(0, length)]

    # conversion
    if conversion == 'float':
        vals = [float(i) for i in vals]
    if conversion == 'int':
        vals = [int(i) for i in vals]
    if conversion == 'str':
        vals = [str(i) for i in vals]

    return vals


# self-defined curve
# text
# circle/ellipse
# polygon
# arrow/line
def set_annotation(args, font, x1beg, n1, d1, axis1len, x2beg, n2, d2, axis2len):

    # get axis of current plot
    ax = plt.gca()

    # plot curve
    if args.curve is not None:

        curvefile = args.curve[0].split(",")
        nf = len(curvefile)

        curvestyle = set_default(args.curvestyle, ',', nf, 'scatter.')
        curvecolor = set_default(args.curvecolor, ',', nf, 'k')
        curvefacecolor = set_default(args.curvefacecolor, ',', nf, 'k')
        curveedgecolor = set_default(args.curveedgecolor, ',', nf, 'none')
        curvesize = set_default(args.curvesize, ',', nf, 1.0, 'float')
        curvewidth = set_default(args.curvewidth, ',', nf, 1.0, 'float')
        curveorder = set_default(args.curveorder, ',', nf, 9, 'int')

        for i in range(0, nf):

            curve = np.loadtxt(curvefile[i], ndmin=2)  # using ndmin=2 to ensure read as 2d array
            nsp = len(curve)
            x1 = curve[0:nsp, 0]
            x2 = curve[0:nsp, 1]
            px1 = (x1 - x1beg + 0.5 * d1) / (n1 * d1) * axis1len
            px2 = (x2 - x2beg + 0.5 * d2) / (n2 * d2) * axis2len
            curve[0:nsp, 0] = px1
            curve[0:nsp, 1] = px2

            # plot scatter points on the figure
            if 'scatter' in curvestyle[i]:
                ax.scatter(px2,
                           px1,
                           marker=curvestyle[i][7:],
                           facecolor=curvefacecolor[i],
                           zorder=curveorder[i],
                           edgecolor=curveedgecolor[i],
                           s=curvesize[i])
                ax.set_xlim(0, axis2len)
                ax.set_ylim(axis1len, 0)
            # plot line on the figure
            if 'line' in curvestyle[i]:
                extra = Line2D(px2,
                               px1,
                               linestyle=curvestyle[i][4:],
                               zorder=curveorder[i],
                               color=curvecolor[i],
                               linewidth=curvewidth[i])
                ax.add_artist(extra)
            # plot polygon on the figure
            if 'polygon' in curvestyle[i]:
                curve[0:nsp, 0] = px2
                curve[0:nsp, 1] = px1
                extra = Polygon(curve,
                                fill=False,
                                zorder=curveorder[i],
                                edgecolor=curvecolor[i],
                                linewidth=curvewidth[i])
                ax.add_artist(extra)

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
        center1 = x1beg + ((n1 - 1) * d1) / 2.0
        center2 = x2beg + ((n2 - 1) * d2) / 2.0
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
            x1 = float(tloc[0])
            x2 = float(tloc[1])
            px1 = (x1 - x1beg + 0.5 * d1) / (n1 * d1) * axis1len
            px2 = (x2 - x2beg + 0.5 * d2) / (n2 * d2) * axis2len

            if int(textbox[i]) == 1:
                # if textbox is required
                t = ax.text(px2,
                            px1,
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
                t.set_fontproperties(font)
                t.set_size(textsize[i])
                bb = t.get_bbox_patch()
                bb.set_boxstyle(pad=float(textboxpad[i]))

            else:
                # plain text only
                t = ax.text(px2,
                            px1,
                            text[i],
                            zorder=textorder[i],
                            color=textcolor[i],
                            fontproperties=font,
                            fontsize=textsize[i],
                            horizontalalignment='center',
                            verticalalignment='center',
                            rotation=textrotation[i])
                t.set_fontproperties(font)
                t.set_size(textsize[i])

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

            tail1 = (float(arrowloc[0]) - x1beg + 0.5 * d1) / (n1 * d1) * axis1len
            tail2 = (float(arrowloc[1]) - x2beg + 0.5 * d2) / (n2 * d2) * axis2len
            head1 = (float(arrowloc[2]) - x1beg + 0.5 * d1) / (n1 * d1) * axis1len
            head2 = (float(arrowloc[3]) - x2beg + 0.5 * d2) / (n2 * d2) * axis2len
            if arrowstyle[i] != '-':
                ax.annotate('',
                            xytext=(tail2, tail1),
                            xy=(head2, head1),
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
                r1 = scircle[2]
                r2 = scircle[3]
                angle = scircle[4]
            else:
                print('circle specification error')
                exit()

            extra = Ellipse(xy=[(float(c2) - x2beg + 0.5 * d2) / (n2 * d2) * axis2len,
                                (float(c1) - x1beg + 0.5 * d1) / (n1 * d1) * axis1len],
                            fill=True,
                            facecolor=circlefacecolor[i],
                            edgecolor=circleedgecolor[i],
                            alpha=circlealpha[i],
                            width=(float(r1) - x2beg + 0.5 * d2) / (n2 * d2) * axis2len,
                            height=(float(r2) - x1beg + 0.5 * d1) / (n1 * d1) * axis1len,
                            linestyle=circlelinestyle[i],
                            linewidth=circlelinewidth[i],
                            angle=float(angle),
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
                px1 = (float(x1[j]) - x1beg + 0.5 * d1) / (n1 * d1) * axis1len
                px2 = (float(x2[j]) - x2beg + 0.5 * d2) / (n2 * d2) * axis2len
                spolygon[j, 0] = px2
                spolygon[j, 1] = px1

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

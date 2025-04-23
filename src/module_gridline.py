'''
    Module:
        Set grid lines
'''
import matplotlib.pyplot as plt


def set_gridline(args):

    ax = plt.gca()

    if args.grid1:
        # grid line width
        if args.grid1width is None:
            grid1width = float(args.tickmajorwid)
        else:
            grid1width = float(args.grid1width)
        # add grid
        ax.grid(which='major',
                axis='y',
                linestyle=args.grid1style,
                color=args.grid1color,
                linewidth=grid1width,
                zorder=10)
        
    if args.grid2:
        # grid line width
        if args.grid2width is None:
            grid2width = float(args.tickmajorwid)
        else:
            grid2width = float(args.grid2width)
        # add grid
        ax.grid(which='major',
                axis='x',
                linestyle=args.grid2style,
                color=args.grid2color,
                linewidth=grid2width,
                zorder=10)

'''
    Module:
        Clip the input data
'''
import numpy as np
from module_message import msg_clip


def set_clip(args, data, which='fore', dmin=0, dmax=1):

    # data value range
    dlen = dmax - dmin

    if which == 'fore':

        pmin = dmin + (1.0 - float(args.cperc) / 100.0) * 0.5 * dlen
        pmax = dmax - (1.0 - float(args.cperc) / 100.0) * 0.5 * dlen

        if args.cmin is None:
            plot_min_value = -float(args.clip) if args.clip is not None else pmin
        else:
            plot_min_value = float(args.cmin)

        if args.cmax is None:
            plot_max_value = float(args.clip) if args.clip is not None else pmax
        else:
            plot_max_value = float(args.cmax)

    if which == 'back':

        pmin = dmin + (1.0 - float(args.backcperc) / 100.0) * 0.5 * dlen
        pmax = dmax - (1.0 - float(args.backcperc) / 100.0) * 0.5 * dlen

        if args.backcmin is None:
            plot_min_value = -float(args.backclip) if args.backclip is not None else pmin
        else:
            plot_min_value = float(args.backcmin)

        if args.backcmax is None:
            plot_max_value = float(args.backclip) if args.backclip is not None else pmax
        else:
            plot_max_value = float(args.backcmax)

    msg_clip(plot_min_value, plot_max_value)

    return plot_min_value, plot_max_value

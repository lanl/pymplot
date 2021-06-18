'''
    Module:
        Clip the input data
'''
import numpy as np

def set_clip(args, data, which='fore', dmin=0, dmax=1):

    # data value range
    dlen = dmax - dmin

    if which == 'fore':

        pmin = dmin + (1.0 - float(args.cperc) / 100.0) * 0.5 * dlen
        pmax = dmax - (1.0 - float(args.cperc) / 100.0) * 0.5 * dlen

        # minimum plot value
        if args.cmin is None:
            if args.clip is None:
                plot_min_value = pmin
            else:
                plot_min_value = -float(args.clip)
        else:
            plot_min_value = float(args.cmin)

        # maximum plot value
        if args.cmax is None:
            if args.clip is None:
                plot_max_value = pmax
            else:
                plot_max_value = float(args.clip)
        else:
            plot_max_value = float(args.cmax)

        # print clip information
        print('plot range  ', "{:e}".format(plot_min_value), ' -- ', "{:e}".format(plot_max_value))

    if which == 'back':

        pmin = dmin + (1.0 - float(args.backcperc) / 100.0) * 0.5 * dlen
        pmax = dmax - (1.0 - float(args.backcperc) / 100.0) * 0.5 * dlen

        # minimum plot value
        if args.backcmin is None:
            if args.backclip is None:
                plot_min_value = pmin
            else:
                plot_min_value = -float(args.backclip)
        else:
            plot_min_value = float(args.backcmin)

        # maximum plot value
        if args.backcmax is None:
            if args.backclip is None:
                plot_max_value = pmax
            else:
                plot_max_value = float(args.backclip)
        else:
            plot_max_value = float(args.backcmax)

        # print clip information
        print('plot range  ', "{:e}".format(plot_min_value), ' -- ', "{:e}".format(plot_max_value))

    return plot_min_value, plot_max_value

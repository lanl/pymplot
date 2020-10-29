## set plot minimum and maximum values

from numpy import *


def set_clip(args, data, fore_or_back='fore'):

    # data min and max
    if isnan(sum(data)) == True:
        udata = data[~isnan(data)]
        dmin = udata.min()
        dmax = udata.max()
    else:
        dmin = data.min()
        dmax = data.max()
    dlen = dmax - dmin

    if fore_or_back == 'fore':

        pmin = dmin + (1.0 - float(args.cperc) / 100.0) * 0.5 * dlen
        pmax = dmax - (1.0 - float(args.cperc) / 100.0) * 0.5 * dlen

        # minimum plot value
        if len(args.cmin) == 0:
            if len(args.clip) == 0:
                plot_min_value = pmin
            else:
                plot_min_value = -float(args.clip)
        else:
            plot_min_value = float(args.cmin)

        # maximum plot value
        if len(args.cmax) == 0:
            if len(args.clip) == 0:
                plot_max_value = pmax
            else:
                plot_max_value = float(args.clip)
        else:
            plot_max_value = float(args.cmax)

        # print clip information
        print('plot range  ', plot_min_value, ' -- ', plot_max_value)

    if fore_or_back == 'back':

        pmin = dmin + (1.0 - float(args.backcperc) / 100.0) * 0.5 * dlen
        pmax = dmax - (1.0 - float(args.backcperc) / 100.0) * 0.5 * dlen

        # minimum plot value
        if len(args.backcmin) == 0:
            if len(args.backclip) == 0:
                plot_min_value = pmin
            else:
                plot_min_value = -float(args.backclip)
        else:
            plot_min_value = float(args.backcmin)

        # maximum plot value
        if len(args.backcmax) == 0:
            if len(args.backclip) == 0:
                plot_max_value = pmax
            else:
                plot_max_value = float(args.backclip)
        else:
            plot_max_value = float(args.backcmax)

        # print clip information
        print('plot range  ', plot_min_value, ' -- ', plot_max_value)

    return plot_min_value, plot_max_value

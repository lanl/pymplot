'''
    Module:
        Output figures
'''
import matplotlib.pyplot as plt
import os
from datetime import datetime
import matplotlib as mplt
import numpy as np
from module_datatype import *
from module_message import msg_input, msg_output, msg_done, fatal, warn, error

# Global store for injected numpy arrays (used by function API)
_data_store = {}


def inject_data(which, data):
    """Store a numpy array to be used by read_array instead of a file."""
    _data_store[which] = data


def clear_injected():
    """Remove all injected arrays (call after run() returns)."""
    _data_store.clear()


# read array
def read_array(args, which='fore', dim=2):

    # Check if data was injected directly (function API path)
    if which in _data_store:
        w = _data_store[which].copy().astype(np.float64)
        if w.ndim == 2:
            n1, n2 = w.shape
            n3 = None
        else:
            n1, n2, n3 = w.shape

        if np.isnan(np.sum(w)):
            u = w[~np.isnan(w)]
            dmin = u.min() if u.size > 0 else 0.0
            dmax = u.max() if u.size > 0 else 1.0
        else:
            dmin = w.min()
            dmax = w.max()

        if which == 'mask':
            w[w == 0] = np.NaN

        if which == 'fore' and args.colorscale is not None:
            w = w * float(args.colorscale)
        if which == 'back' and args.backcolorscale is not None:
            w = w * float(args.backcolorscale)

        if which == 'fore' and args.norm == 'log':
            w = np.log(w) / np.log(float(args.base))

        if dim == 2:
            return w, n1, n2, dmin, dmax
        else:
            return w, n1, n2, n3, dmin, dmax

    # foreground, background, or mask
    if which == 'fore':
        infile = args.infile[0]
    if which == 'back':
        infile = args.background
    if which == 'mask':
        infile = args.mask

    if infile is None:
        if dim == 2:
            return None, None, None, None, None
        else:
            return None, None, None, None, None, None

    # check file existence
    if not os.path.exists(infile):
        fatal(f'file not found: {infile}')

    # size
    fsize = os.path.getsize(infile)
    datatype = args.datatype
    if datatype == 'double':
        fsize = fsize / 8
    if datatype == 'float':
        fsize = fsize / 4
    if datatype == 'int':
        fsize = fsize / 2

    if args.n1 is None:
        fatal('n1 must be specified (number of grid points along axis 1)')

    n1 = args.n1
    if args.n2 is None:
        n2 = int(fsize * 1.0 / n1)
    else:
        n2 = args.n2

    if dim == 3:
        if args.n3 is None:
            n3 = int(fsize * 1.0 / n1 / n2)
        else:
            n3 = args.n3

    # dimensions
    if dim == 2:
        shape = (n1, n2)
    else:
        shape = (n1, n2, n3)

    # data type
    dt = set_datatype(args)

    # read
    w = np.fromfile(infile, count=np.prod(shape), dtype=dt)

    # reshape
    if args.transpose:
        w = np.reshape(w, shape)
    else:
        w = np.reshape(w, shape[::-1])
        w = np.transpose(w)

    # flip
    if args.flip1:
        w = np.flip(w, axis=0)
    if args.flip2:
        w = np.flip(w, axis=1)
    if dim == 3 and args.flip3:
        w = np.flip(w, axis=2)

    # data min and max values
    if np.isnan(np.sum(w)):
        u = w[~np.isnan(w)]
        if u.shape == (0, ):
            fatal(f'dataset {infile} is entirely NaN')
        else:
            dmin = u.min()
            dmax = u.max()
    else:
        dmin = w.min()
        dmax = w.max()

    msg_input(infile, shape, dmin, dmax)

    # nan mask
    if which == 'mask':
        w[w == 0] = np.NaN

    # constant scaling
    if which == 'fore' and args.colorscale is not None:
        w = w * float(args.colorscale)
    if which == 'back' and args.backcolorscale is not None:
        w = w * float(args.backcolorscale)

    # log data
    if which == 'fore' and args.norm == 'log':
        w = np.log(w) / np.log(float(args.base))

    if dim == 2:
        return w, n1, n2, dmin, dmax
    else:
        return w, n1, n2, n3, dmin, dmax


def output(args):

    # In API mode: save if outfile given, otherwise do nothing (caller has figure)
    api_mode = getattr(args, '_api_mode', False)

    if args.outfile is None:
        plt.show()
        msg_done()
        return

    if args.outfile[0] == '...':
        import string, random
        tmpfile = "".join([random.choice(string.ascii_letters + string.digits) for n in range(10)])
        tmpfile = '/tmp/scishow_' + tmpfile + '.pdf'
        plt.savefig(tmpfile, dpi=float(args.dpi), bbox_inches='tight', pad_inches=2.0 / 72.0)
        msg_output(tmpfile)
        msg_done()
        if not api_mode:
            os.system('evince ' + tmpfile + ' &')
            exit()
        plt.close()
        return

    outfile = args.outfile[0].split(',')
    for i in outfile:
        extension = os.path.splitext(i)[1]
        if extension.lower() in {'.pdf', '.jpg', '.jpeg', '.png', '.ps', '.eps', '.tiff'}:
            if not args.imageonly:
                plt.savefig(i, dpi=float(args.dpi), bbox_inches='tight', pad_inches=2.0 / 72.0)
            else:
                plt.savefig(i, dpi=float(args.dpi), pad_inches=0.0)
            msg_output(i)
        else:
            error(f'unsupported output format "{extension}", skipping {i}')

    msg_done()
    plt.close()

    if not api_mode:
        exit()

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

# read array
def read_array(args, which='fore', dim=2):
    
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
        print()
        print(' Error: File', infile, 'does not exists; Exit. ')
        print()
        exit()
    
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
        print()
        print(' Error: n1 should > 0; Exit. ')
        print()
        exit()

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
    w = np.fromfile(infile, count=np.product(shape), dtype=dt)
    
    # reshape
    if args.transpose:
        # if the order is (n1, n2, n3) but the data is actually stored (n3, n2, n1)
        w = np.reshape(w, shape)
    else:
        # python storage is C-like by default, therefore read in reverse order
        # and transpose
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
            print(' Error: Input dataset' + infile + ' is all NaN; Exit. ')
            exit()
        else:
            dmin = u.min()
            dmax = u.max()
    else:
        dmin = w.min()
        dmax = w.max()
    
    # print value range    
    print('input <<    ', infile)
    print('shape       ', shape)
    print('value range ', "{:e}".format(dmin), ' -- ', "{:e}".format(dmax))
    
    # nan mask
    if which == 'mask':
        w[w == 0] = np.NAN

    # constant scaling
    if which == 'fore' and args.colorscale is not None:
        w = w*float(args.colorscale)
    if which == 'back' and args.backcolorscale is not None:
        w = w*float(args.backcolorscale)

    # log data
    if which == 'fore' and args.norm == 'log': 
        w = np.log(w) /  np.log(float(args.base))
    
    if dim == 2:
        return w, n1, n2, dmin, dmax
    else:
        return w, n1, n2, n3, dmin, dmax

def output(args):

    if args.outfile is None:
        # if output file is not specified, then show with native window

        plt.show()

        # time information
        print(' @ ', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        print()

    else:
        if args.outfile[0] == '...':
            # if output is ..., then export to /tmp/*random*.pdf and open

            import string, random
            tempfile = "".join([random.choice(string.ascii_letters + string.digits) for n in range(10)])
            tempfile = '/tmp/scishow_' + tempfile + '.pdf'
            plt.savefig(tempfile, dpi=float(args.dpi), bbox_inches='tight', pad_inches=2.0 / 72.0)
            print('output >>   ', tempfile)
            print()
            os.system('evince ' + tempfile + ' &')
            exit()

        else:
            # if output is specified, then save to file(s)

            outfile = args.outfile[0].split(',')

            for i in outfile:
                extension = os.path.splitext(i)[1]
                if extension.lower() in {'.pdf', '.jpg', '.jpeg', '.png', '.ps', '.eps', '.tiff'}:
                    if not args.imageonly:
                        plt.savefig(i, dpi=float(args.dpi), bbox_inches='tight', pad_inches=2.0 / 72.0)
                    else:
                        plt.savefig(i, dpi=float(args.dpi), pad_inches=0.0)
                    print('output >>   ', i)
                else:
                    print('error: unsupported output figure format ' + extension + ', skip ' + i)

        # time information
        print(' @ ', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        print()

        exit()

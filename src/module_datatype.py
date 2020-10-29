## set data type for data input
import numpy as np


def set_datatype(args):

    # double-precision
    if args.datatype == 'double':
        if args.endian == 'little':
            dt = np.dtype('<f8')
        if args.endian == 'big':
            dt = np.dtype('>f8')

    # single-precision
    if args.datatype == 'float':
        if args.endian == 'little':
            dt = np.dtype('<f4')
        if args.endian == 'big':
            dt = np.dtype('>f4')

    # integer
    if args.datatype == 'int':
        if args.endian == 'little':
            dt = np.dtype('<i2')
        if args.endian == 'big':
            dt = np.dtype('>i2')

    return dt

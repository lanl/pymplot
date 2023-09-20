'''
    Module:
        Tool functions
'''
import numpy as np
import argparse


## round number x to base with specified precision
def rounddecimalbase(x, base, prec=2):
    return round(base * round(float(x) / float(base)), prec)


## order of magnitude of number x
def orderm(x):
    if x == 0:
        return 0
    else:
        return int(np.floor(np.log10(float(abs(x)))))


## nice representation of number x: round to base 0.5
def nice(x, base=0.5):
    nm = orderm(x)
    xnice = x / 10**nm
    xnice = rounddecimalbase(xnice, base)
    xnice = xnice * 10**nm
    return xnice


## convert string to bool
def str2bool(v):

    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 'on', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'off', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Error: Boolean value expected. ')


## convert strings to bool
def strs2bool(v):

    v = v.split(',')
    n = len(v)
    r = np.zeros(n, dtype=bool)

    for i in range(n):

        if isinstance(v[i], bool):
            r[i] = v[i]
        if v[i].lower() in ('yes', 'true', 'on', 't', 'y', '1'):
            r[i] = True
        elif v[i].lower() in ('no', 'false', 'off', 'f', 'n', '0'):
            r[i] = False
        else:
            raise argparse.ArgumentTypeError('Error: Boolean value expected.')

    return r

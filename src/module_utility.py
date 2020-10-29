## several useful functions

import numpy as np


# round number x to base with specified precision
def rounddecimalbase(x, base, prec=2):
    return round(base * round(float(x) / float(base)), prec)


# order of magnitude of number x
def orderm(x):
    if x == 0:
        return 0
    else:
        return int(np.floor(np.log10(float(abs(x)))))


# nice representation of number x: round to base 0.5
def nice(x, base=0.5):
    nm = orderm(x)
    xnice = x / 10**nm
    xnice = rounddecimalbase(xnice, base)
    xnice = xnice * 10**nm
    return xnice

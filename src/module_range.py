'''
    Module:
        Set plotting range for axis
'''

def set_range(f, n, d, xbeg, xend):
    
    # sampling point begin and end
    sample_beg = float(f)
    sample_end = sample_beg + (n - 1) * d
    
    # limit of axis
    if xbeg is None:
        xbeg = sample_beg
    else:
        xbeg = float(xbeg)
    if xend is None:
        xend = sample_end
    else:
        xend = float(xend)

    if xbeg < sample_beg:
        nbeg = 0
        xbeg = sample_beg
    else:
        nbeg = int(abs(round((xbeg - sample_beg) / d)))
    if abs(xend) > abs(sample_end):
        nend = n
        xend = sample_end
    else:
        nend = int(n - abs(round((sample_end - xend) / d)))

    # check range
    if nend <= nbeg:
        print(' Error: Axis range specification error. ')
        exit()

    # return range
    return sample_beg, sample_end, xbeg, xend, nbeg, nend

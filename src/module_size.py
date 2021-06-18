'''
    Module:
        Set figure size
'''
def set_size(args, n1beg, n1end, n2beg, n2end):

    # base figure size
    figbase = 5.0

    # max/min ratio limit
    limit = 6.0
    
    golden_ratio = 1.0 / 1.61803398875
    nmax = max(n1end - n1beg + 1, n2end - n2beg + 1)
    
    if args.size1 is None:
        ratio = float(n1end - n1beg + 1) / nmax
        if ratio < 1.0 / limit:
            ratio = golden_ratio
        size1 = figbase * ratio
    else:
        size1 = float(args.size1)
    
    if args.size2 is None:
        ratio = float(n2end - n2beg + 1) / nmax
        if ratio < 1.0 / limit:
            ratio = golden_ratio
        size2 = figbase * ratio
    else:
        size2 = float(args.size2)
    
    return size1, size2

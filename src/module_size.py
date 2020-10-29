## set figure size


def set_size(args, n1beg, n1end, n2beg, n2end):

    # base figure size
    figbase = 5.0

    # golden ratio
    golden_ratio = 1.61803398875
    limit = 6.0

    if len(args.size1) == 0:
        figheight = figbase
    else:
        figheight = float(args.size1)

    if len(args.size2) == 0:
        figwidth = figbase
    else:
        figwidth = float(args.size2)

    ratio = (n1end - n1beg + 1) * 1.0 / (n2end - n2beg + 1)
    if len(args.size1) == 0:

        if len(args.size2) != 0:
            figwidth = float(args.size2)
        else:
            figwidth = figbase

        if ratio < 1.0 / limit:
            figheight = figwidth / golden_ratio
        if 1.0 / limit <= ratio and ratio < limit:
            figheight = figwidth * ratio
        if ratio > limit:
            figheight = figwidth * golden_ratio

    ratio = (n2end - n2beg + 1) * 1.0 / (n1end - n1beg + 1)
    if len(args.size1) != 0 and len(args.size2) == 0:

        figheight = float(args.size1)

        if ratio < 1.0 / limit:
            figwidth = figheight / golden_ratio
        if 1.0 / limit <= ratio and ratio < limit:
            figwidth = figheight * ratio
        if ratio > limit:
            figwidth = figheight * golden_ratio

    return figheight, figwidth

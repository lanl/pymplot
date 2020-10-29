## set frame styles
import matplotlib.pyplot as plt


def set_frame(args):

    # get current axis object
    ax = plt.gca()

    # top frame
    if args.topframe == 'off':
        ax.spines['top'].set_visible(False)

    # bottom frame
    if args.bottomframe == 'off':
        ax.spines['bottom'].set_visible(False)

    # left frame
    if args.leftframe == 'off':
        ax.spines['left'].set_visible(False)

    # right frame
    if args.rightframe == 'off':
        ax.spines['right'].set_visible(False)

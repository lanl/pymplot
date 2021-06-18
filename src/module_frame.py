'''
    Module:
        Set frame properties
'''
import matplotlib.pyplot as plt


def set_frame(args):

    # get current axis object
    ax = plt.gca()

    # top frame
    if not args.frametop:
        ax.spines['top'].set_visible(False)

    # bottom frame
    if not args.framebottom:
        ax.spines['bottom'].set_visible(False)

    # left frame
    if not args.frameleft:
        ax.spines['left'].set_visible(False)

    # right frame
    if not args.frameright:
        ax.spines['right'].set_visible(False)

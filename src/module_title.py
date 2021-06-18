'''
    Module:
        Set title
'''
import matplotlib.pyplot as plt


def set_title(args, font):

    # figure is 1x1 in size
    if args.title is not None:

        if args.titlesize is None:
            title_font_size = max(float(args.label1size), float(args.label2size)) + 2
        else:
            title_font_size = float(args.titlesize)

        # title position
        title_x = float(args.titlex)
        title_y = float(args.titley)

        # place title
        ax = plt.gca()
        plt.text(title_x,
                 title_y,
                 args.title,
                 horizontalalignment='center',
                 fontproperties=font,
                 fontsize=title_font_size,
                 transform=ax.transAxes)

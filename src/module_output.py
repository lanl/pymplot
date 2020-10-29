## output files
import matplotlib.pyplot as plt
import os
from datetime import datetime
import matplotlib as mplt

def output(args):

    if len(args.outfile) == 0:
        # if output file is not specified, then show with native window

        plt.show()

        # time information
        print(' @ ', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        print()

    else:
        if args.outfile[0] == '...':
            # if output is ..., then export to /tmp/*random*.pdf and open

            import string, random
            tempfile = "".join([
                random.choice(string.ascii_letters + string.digits)
                for n in range(10)
            ])
            tempfile = '/tmp/scishow_' + tempfile + '.pdf'
            plt.savefig(
                tempfile,
                dpi=float(args.dpi),
                bbox_inches='tight',
                pad_inches=2.0 / 72.0)
            print('output >>   ', tempfile)
            print()
            os.system('evince ' + tempfile + ' &')
            exit()

        else:
            # if output is specified, then save to file(s)

            outfile = args.outfile[0].split(',')

            for i in outfile:
                extension = os.path.splitext(i)[1]
                if extension.lower() in {
                        '.pdf', '.jpg', '.jpeg', '.png', '.ps', '.eps', '.tiff'
                }:
                    if args.imageonly == 0:
                        plt.savefig(
                            i,
                            dpi=float(args.dpi),
                            bbox_inches='tight',
                            pad_inches=2.0 / 72.0)
                    else:
                        plt.savefig(i, dpi=float(args.dpi), pad_inches=0.0)
                    print('output >>   ', i)
                else:
                    print('error: unsupported output figure format ' + extension + ', skip ' + i)

        # time information
        print(' @ ', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        print()

        exit()

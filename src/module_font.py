'''
    Module:
        Set font style for a plot.
'''
import matplotlib.font_manager as fm
import matplotlib as mplt
import inspect
import os
from matplotlib import rcParams
#from matplotlib import type1font


def set_font(args):

    # sans-serif fonts
    if args.font == 'arial':
        basefamily = 'sans-serif'
        basefont = 'Arial'
        fontset = 'custom'

    if args.font == 'helvetica':
        basefamily = 'sans-serif'
        basefont = 'Helvetica'
        fontset = 'custom'

    if args.font == 'plex':
        basefamily = 'sans-serif'
        basefont = 'IBMPlexSans'
        fontset = 'custom'
        
    if args.font == 'noto-sans':
        basefamily = 'sans-serif'
        basefont = 'NotoSans'
        fontset = 'custom'

    # serif fonts
    if args.font == 'times':
        basefamily = 'serif'
        basefont = 'Times New Roman'
        fontset = 'custom'

    if args.font == 'georgia':
        basefamily = 'serif'
        basefont = 'Georgia'
        fontset = 'custom'

    if args.font == 'cm':
        basefamily = 'serif'
        basefont = 'CMU Serif'
        fontset = 'cm'
        
    if args.font == 'noto-serif':
        basefamily = 'serif'
        basefont = 'NotoSerif'
        fontset = 'custom'

    # monospace fonts
    if args.font == 'courier':
        basefamily = 'monospace'
        basefont = 'Courier Prime'
        fontset = 'custom'

    if args.font == 'consolas':
        basefamily = 'monospace'
        basefont = 'Consolas'
        fontset = 'custom'
        
    if args.font == 'noto-mono':
        basefamily = 'monospace'
        basefont = 'NotoSansMono'
        fontset = 'custom'

    # set fonts
    rcParams['font.family'] = basefamily
    rcParams['font.' + basefamily] = basefont
    mplt.rcParams['mathtext.fontset'] = fontset
    mplt.rcParams['mathtext.rm'] = basefont
    mplt.rcParams['mathtext.sf'] = basefont
    mplt.rcParams['mathtext.it'] = basefont + ':italic'
    mplt.rcParams['mathtext.bf'] = basefont + ':bold'

    # enforce type 1 font in output: usually result in larger files than type 3 font
    if not args.type3font:
        mplt.rcParams['pdf.fonttype'] = 42
        mplt.rcParams['ps.fonttype'] = 42
        
    # tex math expressions
    mplt.rcParams['mathtext.default'] = 'it' #'regular'
    #mplt.rcParams["text.usetex"] = True

    # font paths
    srcdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    if not args.type3font:
        font = fm.FontProperties(fname=srcdir + '/fonts_subset/' + basefont + '.ttf')
        fontbold = fm.FontProperties(fname=srcdir + '/fonts_subset/' + basefont + 'Bold.ttf')
    else:
        font = fm.FontProperties(fname=srcdir + '/fonts/' + basefont + '.ttf')
        fontbold = fm.FontProperties(fname=srcdir + '/fonts/' + basefont + 'Bold.ttf')

    return font, fontbold

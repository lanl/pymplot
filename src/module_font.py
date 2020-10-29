## set font style

import matplotlib.font_manager as fm
import matplotlib as mplt
import inspect
import os
from matplotlib import rcParams
from matplotlib import type1font


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

    # monospace fonts
    if args.font == 'courier':
        basefamily = 'monospace'
        basefont = 'Courier Prime'
        fontset = 'custom'

    if args.font == 'consolas':
        basefamily = 'monospace'
        basefont = 'Consolas'
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
    # however, type 3 font is sometimes not accepted for publications
    mplt.rcParams['pdf.fonttype'] = 42
    mplt.rcParams['ps.fonttype'] = 42
    #mplt.rcParams['text.antialiased']=True
    #mplt.rcParams['text.usetex']=True

    # font paths
    srcdir = os.path.dirname(
        os.path.abspath(inspect.getfile(inspect.currentframe())))
    font = fm.FontProperties(
        fname=srcdir + '/fonts_subset/' + basefont + '.ttf')
    fontbold = fm.FontProperties(
        fname=srcdir + '/fonts_subset/' + basefont + 'Bold.ttf')
    #font=fm.FontProperties(fname=srcdir+'/fonts/'+basefont+'.ttf')
    #fontbold=fm.FontProperties(fname=srcdir+'/fonts/'+basefont+'Bold.ttf')

    return font, fontbold

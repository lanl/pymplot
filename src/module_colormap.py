'''
    Module:
        Set colormap for a plot.
'''
import numpy as np
import matplotlib as mpl
import matplotlib.pylab as pl
import matplotlib.pyplot as plt
from pylab import *
from matplotlib.colors import ListedColormap
import matplotlib.colors as colors


## Custom colormaps for matplotlib, adopted from Chris Slocum
'''
	Purpose:
		This program shows how to implement make_cmap which is a function that
		generates a colorbar. If you want to look at different color schemes,
		check out https://kuler.adobe.com/create. 
		The function takes a list of tuples which contain RGB values. The RGB
		values may either be in 8-bit [0 to 255] (in which bit must be set to
		True when called) or arithmetic [0 to 1] (default). make_cmap returns
		a cmap with equally spaced colors. 
		Arrange your tuples so that the first color is the lowest value for the
		colorbar and the last is the highest.
		position contains values from 0 to 1 to dictate the location of each color.
		Chris Slocum
	History:
		20130411 -- Initial version created
		20140313 -- Small changes made and code posted online
		20140320 -- Added the ability to set the position of each color
'''
def make_cmap(colors, position=None, bit=False, ncolor=256):

    bit_rgb = np.linspace(0, 1, 256)
    if position == None:
        position = np.linspace(0, 1, len(colors))
    else:
        if len(position) != len(colors):
            sys.exit("position length must be the same as colors")
        elif position[0] != 0 or position[-1] != 1:
            sys.exit("position must start with 0 and end with 1")
    if bit:
        for i in range(len(colors)):
            colors[i] = (bit_rgb[colors[i][0]], bit_rgb[colors[i][1]], bit_rgb[colors[i][2]])
    cdict = {'red': [], 'green': [], 'blue': []}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))

    cmap = mpl.colors.LinearSegmentedColormap('my_colormap', cdict, ncolor)
    return cmap


## truncate colormaps
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=256):

    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))

    return new_cmap


## set colormap with a uniform alpha
def set_colormap(args, which='foreground'):

    # find corresponding colormap if in matplotlib colormap pool
    if which == 'foreground':
        argcolor = args.colormap
    if which == 'background':
        argcolor = args.backcolormap

    colormap = [color for color in dir(cm) if color == argcolor]

    if colormap == []:
        # default colormap is jet of matplotlib
        colormap = cm.jet
        if argcolor == 'mild':
            colormap = cm.RdYlBu_r
        elif argcolor == 'mild_r':
            colormap = cm.RdYlBu
        elif argcolor == 'spectral':
            colormap = cm.Spectral
        elif argcolor == 'spectral_r':
            colormap = cm.Spectral_r
        elif argcolor == 'ncar':
            colormap = cm.gist_ncar
        elif argcolor == 'ncar_r':
            colormap = cm.gist_ncar_r
        elif argcolor == 'nspec':
            colormap = cm.nipy_spectral
        elif argcolor == 'nspec_r':
            colormap = cm.nipy_spectral_r
    else:
        colormap = colormap[0]
            
    cmap = plt.get_cmap(colormap, args.ncolor)

    colors = []

    # several self-defined colormaps
    if argcolor == 'warmcold':
        colors = [(255, 0, 0), (255, 255, 0), (255, 255, 255), (0, 255, 0), (0, 0, 255)]
        position = [0, 0.25, 0.5, 0.75, 1]
    if argcolor == 'coldwarm':
        colors = [(0, 0, 255), (0, 255, 0), (255, 255, 255), (255, 255, 0), (255, 0, 0)]
        position = [0, 0.25, 0.5, 0.75, 1]

    if argcolor == 'rainbow256':
        colors = [(0, 0, 255), (0, 255, 255), (0, 255, 0), (255, 255, 0), (255, 0, 0)]
        position = [0, 0.25, 0.5, 0.75, 1]
    if argcolor == 'rainbow256_r':
        colors = [(0, 0, 255), (0, 255, 255), (0, 255, 0), (255, 255, 0), (255, 0, 0)]
        colors = list(reversed(colors))
        position = [0, 0.25, 0.5, 0.75, 1]

    if argcolor == 'honeypot':
        colors = [(16, 91, 99), (255, 250, 213), (255, 211, 78), (219, 158, 54), (189, 73, 50)]
        position = [0, 0.25, 0.5, 0.75, 1]
    if argcolor == 'honeypot_r':
        colors = [(16, 91, 99), (255, 250, 213), (255, 211, 78), (219, 158, 54), (189, 73, 50)]
        colors = list(reversed(colors))
        position = [0, 0.25, 0.5, 0.75, 1]

    if argcolor == 'dbwr':
        colors = [(0, 0, 215), (5, 75, 255), (255, 255, 255), (215, 35, 15), (185, 0, 15)]
        position = [0, 0.25, 0.5, 0.75, 1]
    if argcolor == 'dbwr_r':
        colors = [(0, 0, 215), (5, 75, 255), (255, 255, 255), (215, 35, 15), (185, 0, 15)]
        colors = list(reversed(colors))
        position = [0, 0.25, 0.5, 0.75, 1]

    if argcolor == 'kwr':
        colors = [(0, 0, 0), (255, 255, 255), (255, 0, 0)]
        position = [0, 0.5, 1]
    if argcolor == 'kwr_r':
        colors = [(0, 0, 0), (255, 255, 255), (255, 0, 0)]
        colors = list(reversed(colors))
        position = [0, 0.5, 1]

    if argcolor == 'kwyr':
        colors = [(0, 0, 0), (128, 128, 128), (255, 255, 255), (255, 255, 0), (255, 0, 0)]
        position = [0, 0.25, 0.5, 0.75, 1]
    if argcolor == 'kwyr_r':
        colors = [(0, 0, 0), (128, 128, 128), (255, 255, 255), (255, 255, 0), (255, 0, 0)]
        colors = list(reversed(colors))
        position = [0, 0.25, 0.5, 0.75, 1]

    if argcolor == 'ywr':
        colors = [(255, 255, 0), (255, 255, 255), (255, 0, 0)]
        position = [0, 0.5, 1]
    if argcolor == 'ywr_r':
        colors = [(255, 255, 0), (255, 255, 255), (255, 0, 0)]
        colors = list(reversed(colors))
        position = [0, 0.5, 1]

    if argcolor == 'ryw':
        colors = [(255, 0, 0), (255, 255, 0), (255, 255, 255)]
        position = [0, 0.5, 1]
    if argcolor == 'ryw_r':
        colors = [(255, 0, 0), (255, 255, 0), (255, 255, 255)]
        colors = list(reversed(colors))
        position = [0, 0.5, 1]

    if argcolor == 'ry':
        colors = [(255, 0, 0), (255, 255, 0)]
        position = [0, 1]
    if argcolor == 'ry_r':
        colors = [(255, 0, 0), (255, 255, 0)]
        colors = list(reversed(colors))
        position = [0, 1]

    if argcolor == 'by':
        colors = [(0, 0, 255), (255, 255, 0)]
        position = [0, 1]
    if argcolor == 'by_r':
        colors = [(0, 0, 255), (255, 255, 0)]
        colors = list(reversed(colors))
        position = [0, 1]

    if argcolor == 'by':
        colors = [(0, 0, 255), (255, 255, 0)]
        position = [0, 1]
    if argcolor == 'by_r':
        colors = [(0, 0, 255), (255, 255, 0)]
        colors = list(reversed(colors))
        position = [0, 1]
        
    if argcolor == 'rainbowcmyk':
        colors = [(int(0.00*255), int(0.00*255), int(0.51*255)), 
                  (int(0.00*255), int(0.81*255), int(1.00*255)), 
                  (int(0.87*255), int(1.00*255), int(0.12*255)), 
                  (int(1.00*255), int(0.20*255), int(0.00*255)), 
                  (int(0.51*255), int(0.00*255), int(0.00*255))]
        position = [0.00, 0.34, 0.61, 0.84, 1.00]
    if argcolor == 'rainbowcmyk_r':
        colors = [(int(0.00*255), int(0.00*255), int(0.51*255)), 
                  (int(0.00*255), int(0.81*255), int(1.00*255)), 
                  (int(0.87*255), int(1.00*255), int(0.12*255)), 
                  (int(1.00*255), int(0.20*255), int(0.00*255)), 
                  (int(0.51*255), int(0.00*255), int(0.00*255))]
        colors = list(reversed(colors))
        position = [0.00, 0.34, 0.61, 0.84, 1.00]

    # make self-defined colormap
    if colors != []:
        cmap = make_cmap(colors, position, bit=True)

    # truncate colormap if necessary
    if which == 'foreground':
        colormap = truncate_colormap(cmap, float(args.ctruncbeg), float(args.ctruncend))
    if which == 'background':
        colormap = truncate_colormap(cmap, float(args.backctruncbeg), float(args.backctruncend))

    # make nan values transparent
    colormap.set_bad('white', 0.0)
    
    # resampled
    if args.ncolor != 256:
    	colormap = colormap.resampled(args.ncolor)
    
    return colormap

## set colormap with non-uniform alphas
def set_colormap_alpha(args, colormap, cmin, cmax, which='foreground'):

    if which == 'foreground':
        if args.alphas is not None:
            n = colormap.N
            crange = np.linspace(cmin, cmax, n)
            alphas = args.alphas.split(",")
            l = len(alphas)
            v = np.zeros(l)
            a = np.zeros(l)
            for i in range(0, l):
                d = alphas[i].split(':')
                v[i] = d[0]
                a[i] = d[1]
            alpha_colormap = colormap(np.arange(n))
            alpha_colormap[:, -1] = np.interp(crange, v, a)
            return ListedColormap(alpha_colormap)
        else:
            return colormap

    if which == 'background':
        if args.backalphas is not None:
            n = colormap.N
            crange = np.linspace(cmin, cmax, n)
            alphas = args.backalphas.split(",")
            l = len(alphas)
            v = np.zeros(l)
            a = np.zeros(l)
            for i in range(0, l):
                d = alphas[i].split(':')
                v[i] = d[0]
                a[i] = d[1]
            alpha_colormap = colormap(np.arange(n))
            alpha_colormap[:, -1] = np.interp(crange, v, a)
            return ListedColormap(alpha_colormap)
        else:
            return colormap

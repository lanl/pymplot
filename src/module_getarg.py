# get arguments
import argparse


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


## arguments -- general
def getarg_general(parser, dim=2):

    flags = parser.add_argument_group('required arguments')

    flags.add_argument('-in', '--infile', type=str, help='Input raw binary file', nargs='+', required=True)
    parser.add_argument('-out',
                        '--outfile',
                        type=str,
                        help='Output figure file',
                        nargs='+',
                        required=False,
                        default='')
    flags.add_argument('-n1', '--n1', type=int, help='Number of points in first dimension', required=True)

    parser.add_argument('-imageonly',
                        '--imageonly',
                        type=int,
                        help='Save only plotting region, no frame or axes',
                        required=False,
                        default=0)
    if dim == 2:
        parser.add_argument('-n2',
                            '--n2',
                            type=int,
                            help='Number of points in second dimension',
                            required=False,
                            default=0)
    parser.add_argument('-dtype',
                        '--datatype',
                        type=str,
                        help='Input file data type, =double, =float (default) or =int',
                        required=False,
                        default='float')
    parser.add_argument('-endian',
                        '--endian',
                        type=str,
                        help='Endianness of data, =big or =little (default)',
                        required=False,
                        default='little')
    parser.add_argument('-transpose',
                        '--transpose',
                        type=int,
                        help='Plot transposed data, =0 by default or =1',
                        required=False,
                        default=0)
    parser.add_argument('-size1', '--size1', type=str, help='Axis 1 size in inch', required=False, default='')
    parser.add_argument('-size2', '--size2', type=str, help='Axis 2 size in inch', required=False, default='')
    parser.add_argument('-mask',
                        '--mask',
                        type=str,
                        help='Input mask file, =0 no draw, =1 draw',
                        required=False,
                        default='')
    parser.add_argument('-flip1', '--flip1', type=int, help='Flip first dimension', required=False, default=0)
    parser.add_argument('-flip2',
                        '--flip2',
                        type=int,
                        help='Flip second dimension',
                        required=False,
                        default=0)

    if dim == 3:
        flags.add_argument('-n2',
                           '--n2',
                           type=int,
                           help='Number of points in second dimension',
                           required=False)
        parser.add_argument('-n3',
                            '--n3',
                            type=int,
                            help='Number of points in the third dimension',
                            required=False,
                            default=0)
        parser.add_argument('-slice1',
                            '--slice1',
                            type=str,
                            help='Position of slice on axis 1, in consistent with actual value of axis',
                            required=False,
                            default='')
        parser.add_argument('-slice2',
                            '--slice2',
                            type=str,
                            help='Position of slice on axis 2, in consistent with actual value of axis',
                            required=False,
                            default='')
        parser.add_argument('-slice3',
                            '--slice3',
                            type=str,
                            help='Position of slice on axis 3, in consistent with actual value of axis',
                            required=False,
                            default='')
        parser.add_argument('-size3', '--size3', type=str, help='Axis 3 in inch', required=False, default='')
        parser.add_argument('-flip3',
                            '--flip3',
                            type=int,
                            help='Flip third dimension',
                            required=False,
                            default=0)

    # common options
    parser.add_argument('-dpi',
                        '--dpi',
                        type=str,
                        help='DPI of generated figure, =300 (default)',
                        required=False,
                        default='300')
    parser.add_argument('-norm',
                        '--norm',
                        type=str,
                        help='''Norm of data values, =linear (default) or =log;
for log norm, the base is 10, and values of the data
will be converted to logarithmic norm,
the options cmin/cmax should be based on the
power of the converted data,
rather than based on the values of the raw data''',
                        required=False,
                        default='linear')
    parser.add_argument('-background',
                        '--background',
                        type=str,
                        help='''File name of the raw binary file that is expected
to be plotted as background in the figure; the background
raw binary file must have the same dimension (i.e., number
of points in each dimesnion), the data ordering (transpoe or
no transpose) and data type with the raw binary file
specified in --infile option; for options to specify the
background plot, please refer to relavent options''',
                        required=False,
                        default=None)
    parser.add_argument('-font',
                        '--font',
                        type=str,
                        help='''Font style for the plot, valid options include:
    =arial (Arial),
    =helvetica (Helvetica) (default, usually smaller file size),
    =times (Times New Roman),
    =courier (Courier New),
    =consolas (Consolas)''',
                        required=False,
                        default='arial')
    #    =georgia (Georgia),
    return parser


## arguments -- axis
def getarg_axis(parser, dim):

    parser.add_argument('-d1',
                        '--d1',
                        type=str,
                        help='''Sample interval in the first dimension, =1 by default;
if <0, then axis tick values descending''',
                        required=False,
                        default='1.0')
    parser.add_argument('-d2',
                        '--d2',
                        type=str,
                        help='''Sample interval in the second dimension, =1 by default;
if <0, then axis tick values descending''',
                        required=False,
                        default='1.0')
    parser.add_argument('-f1',
                        '--f1',
                        type=str,
                        help='First sample value in the first dimension, =0 by default',
                        required=False,
                        default='0.0')
    parser.add_argument('-f2',
                        '--f2',
                        type=str,
                        help='First sample value in the second dimension, =0 by default',
                        required=False,
                        default='0.0')
    parser.add_argument('-label1',
                        '--label1',
                        type=str,
                        help='Axis label of the first dimension; Use $$ to render TeX symbols',
                        required=False,
                        default='Axis 1')
    parser.add_argument('-label2',
                        '--label2',
                        type=str,
                        help='Axis label of the second dimension; Use $$ to render TeX symbols',
                        required=False,
                        default='Axis 2')
    parser.add_argument('-label1loc',
                        '--label1loc',
                        type=str,
                        help='Axis 1 label location, =left (default) or =right',
                        required=False,
                        default='left')
    parser.add_argument('-label2loc',
                        '--label2loc',
                        type=str,
                        help='Axis 2 label location, =top (default) or =bottom',
                        required=False,
                        default='top')
    parser.add_argument('-label1pad',
                        '--label1pad',
                        type=str,
                        help='Axis 1 label padding size from tick labels, unit is point',
                        required=False,
                        default='7.5')
    parser.add_argument('-label2pad',
                        '--label2pad',
                        type=str,
                        help='Axis 2 label padding size from tick labels, unit is point',
                        required=False,
                        default='7.5')
    parser.add_argument('-label1size',
                        '--label1size',
                        type=str,
                        help='Axis 1 label font size, =16 by default, unit is point',
                        required=False,
                        default='16.0')
    parser.add_argument('-label2size',
                        '--label2size',
                        type=str,
                        help='Axis 2 label font size, =16 by default, unit is point',
                        required=False,
                        default='16.0')
    parser.add_argument('-x1beg',
                        '--x1beg',
                        type=str,
                        help='Axis 1 begin value, in consistent with the real value of axis',
                        required=False,
                        default='')
    parser.add_argument('-x1end',
                        '--x1end',
                        type=str,
                        help='Axis 1 end value, in consistent with the real value of axis',
                        required=False,
                        default='')
    parser.add_argument('-x2beg',
                        '--x2beg',
                        type=str,
                        help='Axis 2 begin value, in consistent with the real value of axis',
                        required=False,
                        default='')
    parser.add_argument('-x2end',
                        '--x2end',
                        type=str,
                        help='Axis 2 end value, in consistent with the real value of axis',
                        required=False,
                        default='')
    parser.add_argument('-reverse1',
                        '--reverse1',
                        type=str2bool,
                        help='''Axis 1 tick values descend from top to bottom,
=0 by default (values ascend from top to bottom) or =1;
only valid in showmatrix, showcontour, showwiggle and showgraph''',
                        required=False,
                        default='n')
    parser.add_argument('-reverse2',
                        '--reverse2',
                        type=str2bool,
                        help='''Axis 2 tick values descend from left to right,
=0 by default (values ascend from left to right) or =1;
only valid in showmatrix, showcontour, showwiggle and showgraph''',
                        required=False,
                        default='n')

    if dim == 3:

        parser.add_argument('-d3',
                            '--d3',
                            type=str,
                            help='''Sample interval in the third dimension, =1 by default;
if <0, then axis tick values descending''',
                            required=False,
                            default='1.0')
        parser.add_argument('-f3',
                            '--f3',
                            type=str,
                            help='First sample value in the third dimension, =0 by default',
                            required=False,
                            default='0.0')
        parser.add_argument('-label3',
                            '--label3',
                            type=str,
                            help='Axis label of the third dimension; Use $$ to render TeX symbols',
                            required=False,
                            default='Axis 3')
        parser.add_argument('-label3pad',
                            '--label3pad',
                            type=str,
                            help='Axis 3 label padding size from tick labels, unit is point',
                            required=False,
                            default='7.5')
        parser.add_argument('-label3size',
                            '--label3size',
                            type=str,
                            help='Axis 3 label font size, =16 by default, unit is point',
                            required=False,
                            default='16.0')
        parser.add_argument('-x3beg',
                            '--x3beg',
                            type=str,
                            help='Axis 3 begin value, in consistent with the real value of axis',
                            required=False,
                            default='')
        parser.add_argument('-x3end',
                            '--x3end',
                            type=str,
                            help='Axis 3 end value, in consistent with the real value of axis',
                            required=False,
                            default='')

    return parser


## arguments -- color
def getarg_color(parser):

    parser.add_argument('-cscale',
                        '--colorscale',
                        type=str,
                        help='Scale the value of data by multiplying the scalar',
                        required=False,
                        default='1.0')
    parser.add_argument(
        '-clip',
        '--clip',
        type=str,
        help='Set the color rendering range to +-clip; Default is the whole array value range',
        required=False,
        default='')
    parser.add_argument('-cmin',
                        '--cmin',
                        type=str,
                        help='Color rendering range lower bound; Default is min(data)',
                        required=False,
                        default='')
    parser.add_argument('-cmax',
                        '--cmax',
                        type=str,
                        help='Color randering range upper bound; Default is max(data)',
                        required=False,
                        default='')
    parser.add_argument('-cperc',
                        '--cperc',
                        type=str,
                        help='''Percentage of color rendering, the final color rendering range is
cperc*[min(data),max(data)]; Valid value is [0,100]; =100 by default''',
                        required=False,
                        default='100.0')
    parser.add_argument('-color',
                        '--colormap',
                        type=str,
                        help='''Figure colormap;
Available colormaps:
    >> Uniform sequential --
        viridis, inferno, plasma, magma,
    >> Sequential --
        Blues, BuGn, BuPu,
        GnBu, Greens, Greys, Oranges, OrRd,
        PuBu, PuBuGn, PuRd, Purples, RdPu,
        Reds, YlGn, YlGnBu, YlOrBr, YlOrRd,
    >> Sequential (2) --
        afmhot, autumn, bone, cool, copper,
        gist_heat, gray, hot, pink, spring, summer, winter,
    >> Diverging --
        BrBG, bwr, coolwarm, PiYG, PRGn, PuOr,
        RdBu, RdGy, RdYlBu, RdYlGn, Spectral, seismic,
    >> Qualitative --
        Accent, Dark2, Paired, Pastel1,
        Pastel2, Set1, Set2, Set3,
    >> Miscellaneous --
        gist_earth,terrain, ocean, gist_stern,
        brg, CMRmap, cubehelix, gnuplot, gnuplot2, gist_ncar,
        nipy_spectral, jet, rainbow, gist_rainbow,hsv, flag, prism.
    >> Self-defined --
        rainbow256 is a linearly-changing colormap ranging from blue to red;
        kgbwr is a self-defined blue-white-red colormap;
        kwr is a linearly-changing colormap ranging from black, white to red;
        kwyr is a linearly-changing colormap ranging from black, white, yellow to red;
        coldwarm and warmcold is a diverging colormap that changes
        from cold colors to warm colors or vice versa;
Default colormap is jet''',
                        required=False,
                        default='jet')
    parser.add_argument('-ncolor',
                        '--ncolor',
                        type=int,
                        help='Number of discrete colors',
                        required=False,
                        default=256)
    parser.add_argument('-ctruncbeg',
                        '--ctruncbeg',
                        type=str,
                        help='Figure colormap truncation begin; Valid value is [0,1], =0 by default',
                        required=False,
                        default='0.0')
    parser.add_argument('-ctruncend',
                        '--ctruncend',
                        type=str,
                        help='Figure colormap truncation end; Valid value is [0,1], =1.0 by default',
                        required=False,
                        default='1.0')
    parser.add_argument('-backctruncbeg',
                        '--backctruncbeg',
                        type=str,
                        help='Figure colormap truncation begin; Valid value is [0,1], =0 by default',
                        required=False,
                        default='0.0')
    parser.add_argument('-backctruncend',
                        '--backctruncend',
                        type=str,
                        help='Figure colormap truncation end; Valid value is [0,1], =1.0 by default',
                        required=False,
                        default='1.0')
    #    parser.add_argument(
    #        '-alpha',
    #        '--alpha',
    #        type=str,
    #        help='''Image transparency in the range [0,1], =1.0 default;
    # To plot background image, set this value <1 so that background image can be shown''',
    #        required=False,
    #        default='1.0')
    #
    #    parser.add_argument(
    #        '-alphaclipmin',
    #        '--alphaclipmin',
    #        type=str,
    #        help='''Below this value alpha = 0, for foreground image''',
    #        required=False,
    #        default='')
    #    parser.add_argument(
    #        '-alphaclipmax',
    #        '--alphaclipmax',
    #        type=str,
    #        help='''Above this value alpha = 0, for foreground image''',
    #        required=False,
    #        default='')
    #    parser.add_argument(
    #        '-backalphaclipmin',
    #        '--backalphaclipmin',
    #        type=str,
    #        help='''Below this value alpha = 0, for background image''',
    #        required=False,
    #        default='')
    #    parser.add_argument(
    #        '-backalphaclipmax',
    #        '--backalphaclipmax',
    #        type=str,
    #        help='''Above this value alpha = 0, for background image''',
    #        required=False,
    #        default='')
    parser.add_argument(
        '-interp',
        '--interp',
        type=str,
        help='''Image interpolation method, =nearest by default, none, Gaussian, bicubic, etc.;
Not applicable for showvolume''',
        required=False,
        default='none')
    parser.add_argument(
        '-backinterp',
        '--backinterp',
        type=str,
        help='''Background image interpolation method, =nearest by default, none, Gaussian, bicubic, etc.;
Not applicable for showvolume''',
        required=False,
        default='none')
    #    parser.add_argument(
    #        '-backalpha',
    #        '--backalpha',
    #        type=str,
    #        help=
    #        '''Background image transparency in the range [0,1], =1.0 by default;
    # Not applicable for showvolume''',
    #        required=False,
    #        default='1.0')
    parser.add_argument('-backcolor',
                        '--backcolormap',
                        type=str,
                        help='''Background image colormap, =hsv by default;
Not applicable for showvolume (the following 'back' options as well)''',
                        required=False,
                        default='hsv')
    parser.add_argument('-backclip',
                        '--backclip',
                        type=str,
                        help='''Set the color rendering range to +-clip for background image;
Default is the whole array value range''',
                        required=False,
                        default='')
    parser.add_argument('-backcmin',
                        '--backcmin',
                        type=str,
                        help='Color rendering range lower bound; Default is min(backdata)',
                        required=False,
                        default='')
    parser.add_argument('-backcmax',
                        '--backcmax',
                        type=str,
                        help='Color rendering range upper bound; Default is max(backdata)',
                        required=False,
                        default='')
    parser.add_argument('-backcperc',
                        '--backcperc',
                        type=str,
                        help='''Percentage of color rendering, the final color rendering range is
backcperc*[min(backdata),max(backdata)]; Valid value is [0,100]; =100 by default''',
                        required=False,
                        default='100.0')
    parser.add_argument('-nlcolor',
                        '--nlcolor',
                        type=str,
                        help='''Nonlinear colormap by specifying levels''',
                        required=False,
                        default='')

    parser.add_argument('-alphas',
                        '--alphas',
                        type=str,
                        help='''Specify alphas of foreground colormap for different values;
        The default alpha for all values within the colormap is 1; To specify,
        use the following format: -alphas=v1:a1,v2:a2,...; The alpha curve
        will be interpolated from the given vi:ai pairs''',
                        required=False,
                        default='')
    parser.add_argument('-backalphas',
                        '--backalphas',
                        type=str,
                        help='''Specify alphas of background colormap for different values;
        The default alpha for all values within the colormap is 1; To specify,
        use the following format: -backalphas=v1:a1,v2:a2,...; The alpha curve
        will be interpolated from the given vi:ai pairs''',
                        required=False,
                        default='')

    return parser


## arguments -- colorbar
def getarg_colorbar(parser):

    parser.add_argument('-legend',
                        '--legend',
                        type=str2bool,
                        help='Show colorbar, =n by default (no colorbar), =y with colorbar',
                        required=False,
                        default='n')
    parser.add_argument('-lloc',
                        '--lloc',
                        type=str,
                        help='Collorbar location, = left, right (default),top or bottom',
                        required=False,
                        default='right')
    parser.add_argument('-lpad',
                        '--lpad',
                        type=str,
                        help='Colorbar pad away from plot',
                        required=False,
                        default='')
    parser.add_argument('-unit', '--unit', type=str, help='Colorbar unit', required=False, default=None)
    parser.add_argument('-unitsize',
                        '--unitsize',
                        type=str,
                        help='Colorbar unit font size',
                        required=False,
                        default='')
    parser.add_argument('-ltickbeg',
                        '--ltickbeg',
                        type=str,
                        help='Colorbar tick begin',
                        required=False,
                        default='')
    parser.add_argument('-ltickend',
                        '--ltickend',
                        type=str,
                        help='Colorbar tick end',
                        required=False,
                        default='')
    parser.add_argument('-lheight',
                        '--lheight',
                        type=str,
                        help='Colorbar height in inches',
                        required=False,
                        default='')
    parser.add_argument('-lwidth',
                        '--lwidth',
                        type=str,
                        help='Colorbar width in inches',
                        required=False,
                        default='')
    parser.add_argument('-lmtick',
                        '--lmtick',
                        type=int,
                        help='Number of minor ticks between two major colorbar ticks, =0 by default',
                        required=False,
                        default=0)
    parser.add_argument('-ld', '--ld', type=str, help='Colorbar tick interval', required=False, default='')
    parser.add_argument('-unitpad',
                        '--unitpad',
                        type=str,
                        help='Colorbar unit pad distance',
                        required=False,
                        default='')
    parser.add_argument('-lticksize',
                        '--lticksize',
                        type=str,
                        help='Colorbar tick label font size',
                        required=False,
                        default='')
    parser.add_argument('-lreverse',
                        '--lreverse',
                        type=str2bool,
                        help='Reverse colorbar (=y) or not (=n, default)',
                        required=False,
                        default='n')
    parser.add_argument(
        '-backlegend',
        '--backlegend',
        type=str2bool,
        help='Show colorbar of background image, =n by default (no colorbar), =y with colorbar',
        required=False,
        default='n')

    return parser


## arguments -- tick
def getarg_tick(parser, dim):

    parser.add_argument('-ticks1',
                        '--ticks1',
                        type=str,
                        help='''Mannualy set ticks along axis 1; e.g., 100,'Planet1':300,'Planet2' ''',
                        required=False,
                        nargs='+',
                        default='')
    parser.add_argument('-ticks2',
                        '--ticks2',
                        type=str,
                        help='''Mannualy set ticks along axis 2; format same as ticks1 ''',
                        required=False,
                        nargs='+',
                        default='')
    parser.add_argument('-tick1beg',
                        '--tick1beg',
                        type=str,
                        help='First tick along axis 1',
                        required=False,
                        default='')
    parser.add_argument('-tick2beg',
                        '--tick2beg',
                        type=str,
                        help='First tick along axis 2',
                        required=False,
                        default='')
    parser.add_argument('-tick1end',
                        '--tick1end',
                        type=str,
                        help='Last tick along axis 1',
                        required=False,
                        default='')
    parser.add_argument('-tick2end',
                        '--tick2end',
                        type=str,
                        help='Last tick along axis 2',
                        required=False,
                        default='')
    parser.add_argument('-tick1d',
                        '--tick1d',
                        type=str,
                        help='Tick interval along axis 1',
                        required=False,
                        default='')
    parser.add_argument('-tick2d',
                        '--tick2d',
                        type=str,
                        help='Tick interval along axis 2',
                        required=False,
                        default='')
    parser.add_argument('-mtick1',
                        '--mtick1',
                        type=int,
                        help='Number of minor ticks between two major ticks along axis 1',
                        required=False,
                        default=0)
    parser.add_argument('-mtick2',
                        '--mtick2',
                        type=int,
                        help='Number of minor ticks between two major ticks along axis 2',
                        required=False,
                        default=0)
    parser.add_argument('-tick1size',
                        '--tick1size',
                        type=str,
                        help='Font size of axis 1 tick labels',
                        required=False,
                        default='')
    parser.add_argument('-tick2size',
                        '--tick2size',
                        type=str,
                        help='Font size of axis 2 tick labels',
                        required=False,
                        default='')
    parser.add_argument('-tickmajorlen',
                        '--tickmajorlen',
                        type=str,
                        help='Lenght of major ticks, =5.0 by default, unit is pt. ',
                        required=False,
                        default='5.0')
    parser.add_argument('-tickminorlen',
                        '--tickminorlen',
                        type=str,
                        help='Length of minor ticks, =0.5*tickmajorlen by default, unit is pt. ',
                        required=False,
                        default='')
    parser.add_argument('-tickmajorwid',
                        '--tickmajorwid',
                        type=str,
                        help='Width of major ticks, =1.0 by default, unit is pt. ',
                        required=False,
                        default='1.0')
    parser.add_argument('-tickminorwid',
                        '--tickminorwid',
                        type=str,
                        help='Width of minor ticks, =0.75*tickmajorwid by default, unit is pt. ',
                        required=False,
                        default='')
    parser.add_argument('-tickleft',
                        '--tickleft',
                        type=str2bool,
                        help='Ticks on the left axis, =y (default) or =n',
                        required=False,
                        default='y')
    parser.add_argument('-tickright',
                        '--tickright',
                        type=str2bool,
                        help='Ticks on the right axis, =y or =n (default)',
                        required=False,
                        default='n')
    parser.add_argument('-ticktop',
                        '--ticktop',
                        type=str2bool,
                        help='Ticks on the top axis, =y (default) or =n',
                        required=False,
                        default='y')
    parser.add_argument('-tickbottom',
                        '--tickbottom',
                        type=str2bool,
                        help='Ticks on the bottom axis, =y or =n (default)',
                        required=False,
                        default='n')
    parser.add_argument('-tick1label',
                        '--tick1label',
                        type=str2bool,
                        help='Axis 1 tick labels, =y (default) or n',
                        required=False,
                        default='y')
    parser.add_argument('-tick2label',
                        '--tick2label',
                        type=str2bool,
                        help='Axis 2 tick labels, =y (default) or n',
                        required=False,
                        default='y')
    parser.add_argument('-tick1format',
                        '--tick1format',
                        type=str,
                        help='''Axis 1 tick label format, =plain or =sci by default, or any legal format''',
                        required=False,
                        default='sci')
    parser.add_argument('-tick2format',
                        '--tick2format',
                        type=str,
                        help='''Axis 2 tick label format, =plain or =sci by default, or any legal format''',
                        required=False,
                        default='sci')
    parser.add_argument('-grid1',
                        '--grid1',
                        type=str2bool,
                        help='Grid lines along axis 1, =y or =n (default)',
                        required=False,
                        default='n')
    parser.add_argument('-grid2',
                        '--grid2',
                        type=str2bool,
                        help='Grid lines along axis 2, =y or =n (default)',
                        required=False,
                        default='n')
    parser.add_argument('-grid1width',
                        '--grid1width',
                        type=str,
                        help='Width of grid lines along axis 1',
                        required=False,
                        default='')
    parser.add_argument('-grid2width',
                        '--grid2width',
                        type=str,
                        help='Width of grid lines along axis 2',
                        required=False,
                        default='')
    parser.add_argument('-grid1color',
                        '--grid1color',
                        type=str,
                        help='Color of grid lines along axis 1, =black by default',
                        required=False,
                        default='k')
    parser.add_argument('-grid2color',
                        '--grid2color',
                        type=str,
                        help='Color of grid lines along axis 2, =black by default',
                        required=False,
                        default='k')
    parser.add_argument(
        '-grid1style',
        '--grid1style',
        type=str,
        help='''Style of grid lines along axis 1, available styles: solid or - (default), dashed or --,
dashdot or -., dotted or :''',
        required=False,
        default='-')
    parser.add_argument(
        '-grid2style',
        '--grid2style',
        type=str,
        help='''Style of grid lines along axis 2, available styles: solid or - (default), dashed or --,
dashdot or -., dotted or :''',
        required=False,
        default='-')
    parser.add_argument('-topframe',
                        '--topframe',
                        type=str2bool,
                        help='Show top frame or hide, =y (default) or =n',
                        required=False,
                        default='y')
    parser.add_argument('-bottomframe',
                        '--bottomframe',
                        type=str2bool,
                        help='Show bottom frame or hide, =y (default) or =n',
                        required=False,
                        default='y')
    parser.add_argument('-leftframe',
                        '--leftframe',
                        type=str2bool,
                        help='Show left frame or hide, =y (default) or =n',
                        required=False,
                        default='y')
    parser.add_argument('-rightframe',
                        '--rightframe',
                        type=str2bool,
                        help='Show right frame or hide, =y (default) or =n',
                        required=False,
                        default='y')
    parser.add_argument('-centerframe',
                        '--centerframe',
                        type=str2bool,
                        help='Show center frame or hide, =y or =n (default)',
                        required=False,
                        default='y')
    parser.add_argument('-tick1rot',
                        '--tick1rot',
                        type=str,
                        help='Axis 1 tick label rotation, =0 by default',
                        required=False,
                        default='0.0')
    parser.add_argument('-tick2rot',
                        '--tick2rot',
                        type=str,
                        help='Axis 2 tick label rotation, =0 by default',
                        required=False,
                        default='0.0')

    if dim == 3:

        parser.add_argument('-ticks3',
                            '--ticks3',
                            type=str,
                            help='Mannualy set ticks along axis 3; format same as ticks1 ',
                            required=False,
                            nargs='+',
                            default='')
        parser.add_argument('-tick3beg',
                            '--tick3beg',
                            type=str,
                            help='First tick along axis 3',
                            required=False,
                            default='')
        parser.add_argument('-tick3end',
                            '--tick3end',
                            type=str,
                            help='Last tick along axis 3',
                            required=False,
                            default='')
        parser.add_argument('-tick3d',
                            '--tick3d',
                            type=str,
                            help='Tick interval along axis 3',
                            required=False,
                            default='')
        parser.add_argument('-mtick3',
                            '--mtick3',
                            type=int,
                            help='Number of minor ticks between two major ticks along axis 3',
                            required=False,
                            default=0)
        parser.add_argument('-tick3size',
                            '--tick3size',
                            type=str,
                            help='Font size of axis 3 tick labels',
                            required=False,
                            default='')
        parser.add_argument('-grid3',
                            '--grid3',
                            type=str2bool,
                            help='Grid lines along axis 3, =y or =n (default)',
                            required=False,
                            default='n')
        parser.add_argument('-grid3width',
                            '--grid3width',
                            type=str,
                            help='Width of grid lines along axis 3',
                            required=False,
                            default='')
        parser.add_argument('-grid3color',
                            '--grid3color',
                            type=str,
                            help='Color of grid lines along axis 3',
                            required=False,
                            default='k')
        parser.add_argument(
            '-grid3style',
            '--grid3style',
            type=str,
            help='''Style of grid lines along axis 3, available styles: solid or - (default), dashed or --,
dashdot or -., dotted or :''',
            required=False,
            default='-')
        parser.add_argument('-sliceline1',
                            '--sliceline1',
                            type=str,
                            help='Draw line to indicate slice 1 location on axis 1',
                            required=False,
                            default='on')
        parser.add_argument('-sliceline1width',
                            '--sliceline1width',
                            type=str,
                            help='Slice 1 line width, =1.0 by default',
                            required=False,
                            default='1.0')
        parser.add_argument('-sliceline1color',
                            '--sliceline1color',
                            type=str,
                            help='Slice 1 line color, =black by default',
                            required=False,
                            default='k')
        parser.add_argument('-sliceline1style',
                            '--sliceline1style',
                            type=str,
                            help='Slice 1 line style, =solid by default',
                            required=False,
                            default='-')
        parser.add_argument('-sliceline2',
                            '--sliceline2',
                            type=str,
                            help='Draw line to indicate slice 2 location on axis 2',
                            required=False,
                            default='on')
        parser.add_argument('-sliceline2width',
                            '--sliceline2width',
                            type=str,
                            help='Slice 2 line width, =1.0  by default',
                            required=False,
                            default='1.0')
        parser.add_argument('-sliceline2color',
                            '--sliceline2color',
                            type=str,
                            help='Slice 2 line color, =black by default',
                            required=False,
                            default='k')
        parser.add_argument('-sliceline2style',
                            '--sliceline2style',
                            type=str,
                            help='Slice 2 line style, =solid by default',
                            required=False,
                            default='-')
        parser.add_argument('-sliceline3',
                            '--sliceline3',
                            type=str,
                            help='Draw line to indicate slice 3 location on axis 3',
                            required=False,
                            default='on')
        parser.add_argument('-sliceline3width',
                            '--sliceline3width',
                            type=str,
                            help='Slice 3 line width, =1.0  by default',
                            required=False,
                            default='1.0')
        parser.add_argument('-sliceline3color',
                            '--sliceline3color',
                            type=str,
                            help='Slice 3 line color, =black by default',
                            required=False,
                            default='k')
        parser.add_argument('-sliceline3style',
                            '--sliceline3style',
                            type=str,
                            help='Slice 3 line style, =solid by default',
                            required=False,
                            default='-')
        parser.add_argument('-slicegap',
                            '--slicegap',
                            type=str,
                            help='White space size between slice plots',
                            required=False,
                            default='0.1')
        parser.add_argument('-tick3rot',
                            '--tick3rot',
                            type=str,
                            help='Axis 3 tick label rotation, =0 by default',
                            required=False,
                            default='0.0')
        parser.add_argument('-tick3label',
                            '--tick3label',
                            type=str2bool,
                            help='Axis 3 tick labels, =y (default) or n',
                            required=False,
                            default='y')
        parser.add_argument(
            '-tick3format',
            '--tick3format',
            type=str,
            help='''Axis 3 tick label format, =plain or =sci by default, or any legal format''',
            required=False,
            default='sci')
        parser.add_argument('-ticks',
                            '--ticks',
                            type=str2bool,
                            help='Axis 1,2,3 tick labels, =y (default) or n',
                            required=False,
                            default='y')

    return parser


## arguments -- title
def getarg_title(parser, dim=2):

    parser.add_argument('-title', '--title', type=str, help='Figure title', required=False, default='')
    parser.add_argument('-titlesize',
                        '--titlesize',
                        type=str,
                        help='Font size of title',
                        required=False,
                        default='')

    if dim == 2:
        parser.add_argument('-titlex',
                            '--titlex',
                            type=str,
                            help='Title position x2 (in [0,1]), =0.50 default',
                            required=False,
                            default='0.50')
        parser.add_argument('-titley',
                            '--titley',
                            type=str,
                            help='Title position x1 (in [0,1]), =1.25 default',
                            required=False,
                            default='1.25')

    if dim == 3:
        parser.add_argument('-titlex',
                            '--titlex',
                            type=str,
                            help='Title position x2 (in [0,size3+size2])',
                            required=False,
                            default='')
        parser.add_argument('-titley',
                            '--titley',
                            type=str,
                            help='Title position x1 (in [0,size1+size2])',
                            required=False,
                            default='')

    return parser


## arguments -- annotation
def getarg_annotation(parser):

    # draw open curves
    parser.add_argument('-curve',
                        '--curve',
                        type=str,
                        help='''Files that specify
the coordinates of the curve points. The filenames are separated by comma, i.e.
curve1.txt,curve2.txt,...; the coordinates in each file are in the following format:
px py
px py
px py
...
all in ASCII format. ''',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-curvestyle',
                        '--curvestyle',
                        type=str,
                        help='''Curve plot
style, =scatter(+style): plot given points as scatter points; =line(+style): plot
given points as open continous curve; =polygon: plot given points as closed
continous curve (i.e., polygon). Scatter and line styles follow matplotlib standard names,
e.g., scatter., scattero, scatter^, scatter*, line-, line--, etc..
The styles are separated by comma, e.g., scatter., line-,scatter*,line--,...''',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-curvesize',
                        '--curvesize',
                        type=str,
                        help='''Curve plot size (line width or scatter plot marker area size),
=1.0 default, separated by comma, e.g., 1.0,2.0,3.0.,...''',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-curvewidth',
                        '--curvewidth',
                        type=str,
                        help='''Curve plot size (line width or scatter plot marker area size),
=1.0 default, separated by comma, e.g., 1.0,2.0,3.0.,...''',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-curvecolor',
                        '--curvecolor',
                        type=str,
                        help='''(Face) color of curves, separated by comma, e.g., k,b,w,y,...''',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-curvefacecolor',
                        '--curvefacecolor',
                        type=str,
                        help='''Face color of curves, separated by comma, e.g., k,b,w,y,...''',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-curveedgecolor',
                        '--curveedgecolor',
                        type=str,
                        help='''Edge color of curves, separated by comma, e.g., k,b,w,y,...''',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-curveorder',
                        '--curveorder',
                        type=str,
                        help='''Overlay order of curves, separated by comma, e.g.,10,9,10,6,...;
integer values only; the higher the shallower; =9 by default''',
                        nargs='+',
                        required=False,
                        default='')

    # add text
    parser.add_argument('-text',
                        '--text',
                        type=str,
                        help='Extra text on the plot, separated by comma, e.g.,text1,text2,... ',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-textstyle',
                        '--textstyle',
                        type=str,
                        help='Styles of text, separated by comma',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-textloc',
                        '--textloc',
                        type=str,
                        help='''Coordinates of text (w.r.t. real axis of the base plot),
separated by comma, e.g., x1,y1,x2,y2,... ''',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-textrotation',
                        '--textrotation',
                        type=str,
                        help='Rotation of text, separated by comma, e.g., angle1,angle2,... ',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-textsize',
                        '--textsize',
                        type=str,
                        help='Text font size, =12.0 default, separated by comma',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-textcolor',
                        '--textcolor',
                        type=str,
                        help='Text color, separated by comma, e.g., k,b,w,y',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-textbox',
                        '--textbox',
                        type=str,
                        help='''Bounding box around text, =0 no bounding box (by default) or =1,
separted by comma, e.g., 1,0,1,1,1,... ''',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-textboxedgecolor',
                        '--textboxedgecolor',
                        type=str,
                        help='Text bounding box edge color =y by default, separated by comma ',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-textboxfacecolor',
                        '--textboxfacecolor',
                        type=str,
                        help='Text bounding box face color =y by default, separated by comma ',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-textboxalpha',
                        '--textboxalpha',
                        type=str,
                        help='Text bounding box transparency =0.5 by default, separated by comma ',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-textboxpad',
                        '--textboxpad',
                        type=str,
                        help='Pad size of bounding box around text, =10 by default, separted by comma ',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-textboxstyle',
                        '--textboxstyle',
                        type=str,
                        help='''Text bounding box style, =square by default.
Availalbe choices are:
circle, darrow, larrow, rarrow, round, round4, roundtooth, sawtooth, square. ''',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-textorder',
                        '--textorder',
                        type=str,
                        help='''Overlay order of texts, separated by comma,
e.g.,10,9,10,6,...; integer values only; the higher the shallower;
=10 by default''',
                        nargs='+',
                        required=False,
                        default='')

    # draw arrows/line
    parser.add_argument('-arrow',
                        '--arrow',
                        type=str,
                        help='''The coordinates of arrow tail and head (4 numbers per arrow),
separated by semi-colon, e.g, 1,2,3,4:10,11,12,13:...;
can also be used to draw straight lines by using
arrowstyle=- (none arrow head/tail)''',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-arrowfacecolor',
                        '--arrowfacecolor',
                        type=str,
                        help='Arrow face colors, separted by comma ',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-arrowedgecolor',
                        '--arrowedgecolor',
                        type=str,
                        help='Arrow edge colors, separted by comma ',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-arrowlinestyle',
                        '--arrowlinestyle',
                        type=str,
                        help='''Arrow styles, available choices are:
solid, dashed, dashdot, dotted, separated by comma ''',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-arrowstyle',
                        '--arrowstyle',
                        type=str,
                        help='''Arrow styles, available choices are:
-, ->, -[, |-|, -|>, <-, <->, <|-, <|-|>, fancy, simple, wedge,
separated by comma ''',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-arrowconnect',
                        '--arrowconnect',
                        type=str,
                        help='''Arrow connection styles, separted by semi-colon,
see http://matplotlib.org/users/annotations_guide.html for details. ''',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-arrowwidth',
                        '--arrowwidth',
                        type=str,
                        help='Arrow widths, separted by comma ',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-arroworder',
                        '--arroworder',
                        type=str,
                        help='''Overlay order of arrows, separated by comma,
e.g.,10,9,10,6,...; integer values only; the higher the shallower;
=9 by default''',
                        nargs='+',
                        required=False,
                        default='')

    # draw closed polygon
    parser.add_argument('-polygon',
                        '--polygon',
                        type=str,
                        help='''Polygon coordinates, coordinates in each polygon
is separated by comma, different polygons are separted by semi-colon ''',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-polygonfacecolor',
                        '--polygonfacecolor',
                        type=str,
                        help='Polygon facecolors. ',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-polygonedgecolor',
                        '--polygonedgecolor',
                        type=str,
                        help='Polygon edgecolors. ',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-polygonalpha',
                        '--polygonalpha',
                        type=str,
                        help='Polygon transparencies. ',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-polygonlinestyle',
                        '--polygonlinestyle',
                        type=str,
                        help='Polygon line styles. ',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-polygonlinewidth',
                        '--polygonlinewidth',
                        type=str,
                        help='Polygon line widths. ',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-polygonorder',
                        '--polygonorder',
                        type=str,
                        help='''Overlay order of polygons, separated by comma,
e.g.,10,9,10,6,...; integer values only; the higher the shallower;
=7 by default''',
                        nargs='+',
                        required=False,
                        default='')

    # draw closed circle/ellipse
    parser.add_argument('-circle',
                        '--circle',
                        type=str,
                        help='''Circle center coodinates (x1,x2), circle radius (width and height)
and circle rotation angle, e.g., 3.0,4.0,2.0,1.2,40 specifies a ellipse with center
at (3.0,4.0), with horizontal radius 2.0 and vertical radius 1.2, and rotation angle 40 degress.
Each circle is separted by semi-colon. ''',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-circlefacecolor',
                        '--circlefacecolor',
                        type=str,
                        help='Circle facecolors.',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-circleedgecolor',
                        '--circleedgecolor',
                        type=str,
                        help='Circle edge colors.',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-circlealpha',
                        '--circlealpha',
                        type=str,
                        help='Circle transparencies.',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-circlelinestyle',
                        '--circlelinestyle',
                        type=str,
                        help='Circle line styles. ',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-circlelinewidth',
                        '--circlelinewidth',
                        type=str,
                        help='Circle line widths. ',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-circleorder',
                        '--circleorder',
                        type=str,
                        help='''Overlay order of circles, separated by comma,
e.g.,10,9,10,6,...; integer values only; the higher the shallower;
=8 by default''',
                        nargs='+',
                        required=False,
                        default='')

    return parser


## arguments -- contours
def getarg_contour(parser):

    parser.add_argument('-contours',
                        '--contours',
                        type=str,
                        help='Explicitly specify contours',
                        required=False,
                        nargs='+',
                        default='')
    parser.add_argument('-overlay',
                        '--overlay',
                        type=str2bool,
                        help='Contour overlays on the image, =n (default, no overlay) or =y (overlay)',
                        required=False,
                        default='n')
    parser.add_argument('-contourbeg',
                        '--contourbeg',
                        type=str,
                        help='Start value for contour plot',
                        required=False,
                        default='')
    parser.add_argument('-contourend',
                        '--contourend',
                        type=str,
                        help='End value for contour plot',
                        required=False,
                        default='')
    parser.add_argument('-contourcolor',
                        '--contourcolor',
                        type=str,
                        help='Color of the contour on the image, =black by default',
                        required=False,
                        nargs='+',
                        default='k')
    parser.add_argument('-contourwidth',
                        '--contourwidth',
                        type=str,
                        help='Width of the contour on the image, =1.0 by default',
                        required=False,
                        default='1.0')
    parser.add_argument('-contourstyle',
                        '--contourstyle',
                        type=str,
                        help='Style of the contour on the image, =- (solid, by default)',
                        required=False,
                        default='-')
    parser.add_argument('-contourlevel',
                        '--contourlevel',
                        type=str,
                        help='Interval of the contour on the image',
                        required=False,
                        default='')
    parser.add_argument('-clabelsize',
                        '--clabelsize',
                        type=str,
                        help='Font size of the contour label on the image',
                        required=False,
                        default='')
    parser.add_argument('-clabelcolor',
                        '--clabelcolor',
                        type=str,
                        help='Color of contour label on the image, =black by default',
                        required=False,
                        default='k')
    parser.add_argument('-clabelbackcolor',
                        '--clabelbackcolor',
                        type=str,
                        help='Background color of contour label on the image',
                        required=False,
                        default='')
    parser.add_argument('-mcontour',
                        '--mcontour',
                        type=int,
                        help='Number of minor contours between two numbered major contours',
                        required=False,
                        default=0)
    parser.add_argument('-mcontourwidth',
                        '--mcontourwidth',
                        type=str,
                        help='Minor contour width',
                        required=False,
                        default='')
    parser.add_argument('-mcontourstyle',
                        '--mcontourstyle',
                        type=str,
                        help='Style of the minor contour on the image',
                        required=False,
                        default='-')
    parser.add_argument('-contourfill',
                        '--contourfill',
                        type=str2bool,
                        help='Fill contours =n (default) or =y',
                        required=False,
                        default='n')

    return parser


## arguments -- wiggle
def getarg_wiggle(parser):

    parser.add_argument('-wigglecolor',
                        '--wigglecolor',
                        type=str,
                        help='''Colors of the wiggles, separated by comma for different input files,
e.g., b,k,r,g; if not specified, then using the following
colors in a cycle manner: k,b,r,g,yellow,p,c,b,r,g,yellow,p,c,b,r,g,yellow ''',
                        required=False,
                        nargs='+',
                        default='k,b,r,g,yellow,p,c,b,r,g,yellow,p,c,b,r,g,yellow')
    parser.add_argument('-wigglewidth',
                        '--wigglewidth',
                        type=str,
                        help='''Widths of the wiggles, separated by comma for different input files;
                        =1.0 by default, unit is pt''',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-wigglestyle',
                        '--wigglestyle',
                        type=str,
                        help='''Styles of wiggles, separated by comma for different input files;
                        =- (solid line) by default''',
                        nargs='+',
                        required=False,
                        default='')
    parser.add_argument('-axispad',
                        '--axispad',
                        type=str,
                        help='Pad space of plot =0.1 inch by default',
                        required=False,
                        default='0.1')
    parser.add_argument('-interp1',
                        '--interp1',
                        type=str,
                        help='''Interpolation along axis 1, for wiggle plot;
                        valid only when wiggle is along axis 1''',
                        required=False,
                        default='')
    parser.add_argument('-interp2',
                        '--interp2',
                        type=str,
                        help='''Interpolation along axis 2, for wiggle plot;
                        valid only when wiggle is along axis 2''',
                        required=False,
                        default='')
    parser.add_argument('-overlay',
                        '--overlay',
                        type=int,
                        help='Plot wiggles over corresponding data image',
                        required=False,
                        default=0)
    parser.add_argument('-along',
                        '--along',
                        type=int,
                        help='Plot wiggles along which axis, =1 by default or =2',
                        required=False,
                        default=1)
    parser.add_argument('-fill',
                        '--fill',
                        type=int,
                        help='''Fill the positive or negative part of the wiggles with wiggles' color,
                        =-1 (fill negative polarity), 0 (no fill) or =1 (fill positive polarity);
                        =0 by default''',
                        required=False,
                        default=0)
    parser.add_argument('-every',
                        '--every',
                        type=int,
                        help='Show every N traces to avoid crowd figure; =1 by default (show all traces)',
                        required=False,
                        default=1)

    return parser

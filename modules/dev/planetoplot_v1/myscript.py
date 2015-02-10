def getparseroptions(parser):

    ### I/O
    parser.add_option('-f', '--file',   action='append',dest='file',     type="string",  default=None,  help='[NEEDED] filename. Append: different figures. Comma-separated: same figure (+ possible --operation). Regex OK: use -f "foo*" DONT FORGET QUOTES "" !!!!')
    parser.add_option('-L', '--large',  action='store_true',dest='monster',              default=False, help='speedy version for large files (EXPERIMENTAL)')
    parser.add_option('--seevar',       action='store_true',dest='seevar',               default=False, help='display the list of variables in the file')
    parser.add_option('-t', '--target', action='store',dest='tgt',       type="string",  default=None,  help='destination folder')
    parser.add_option('-S', '--save',   action='store',dest='save',      type="string",  default="gui", help='save mode (gui,png,eps,svg,pdf,txt,html,avi) [gui]')
    parser.add_option('-d', '--display',action='store_false',dest='display',             default=True,  help='do not pop up created images')
    parser.add_option('-O','--output',  action='store',dest='out',       type="string",  default=None,  help='output file name')
    parser.add_option('--rate',         action='store'     ,dest='rate', type="int",   default=None, help='output is a movie along Time dimension [None]')
    parser.add_option('--quality',      action='store_true',dest='quality',            default=False,help='For movie mode: improves movie quality.(slower)')

    ### WHAT I PLOT
    parser.add_option('-v', '--var',    action='append',dest='var',      type="string",  default=None,  help='variable color-shaded')
    parser.add_option('-w', '--with',   action='store',dest='var2',      type="string",  default=None,  help='variable contoured')
    parser.add_option('-a', '--anomaly',action='store_true',dest='anomaly',              default=False, help='compute and plot relative anomaly in %')
    parser.add_option('--mult',         action='store',dest='mult',      type="float",   default=1.,    help='multiplicative factor to plotted field (2718=log)')
    parser.add_option('--add',          action='store',dest='add',       type="float",   default=0.,    help='additive factor to plotted field')
    parser.add_option('-m', '--min',    action='append',dest='vmin',     type="float",   default=None,  help='bounding minimum value [min]')    
    parser.add_option('-M', '--max',    action='append',dest='vmax',     type="float",   default=None,  help='bounding maximum value [max]') 
    parser.add_option('-H', '--hole',   action='store_true',dest='hole',                 default=False, help='holes above max and below min')
    parser.add_option('--nolow',        action='store_true',dest='nolow',                default=False, help='do not plot low |values| [False]')
    parser.add_option('--redope',       action='store',dest='redope',    type="string",  default=None,  help='REDuce OPErators: mint,maxt for the moment [None]')

    ### VERTICAL INTERPOLATION
    parser.add_option('-l', '--level',  action='store',dest='lvl',       type="string",  default="0",   help='level / start,stop,step (-i 2: p,Pa)(-i 3,4: z,km) [0]')
    parser.add_option('-i', '--interp', action='store',dest='itp',       type="int",     default=None,  help='interpolation (2: p, 3: z-amr, 4:z-als, -1)')
    parser.add_option('--intas',        action='store',dest='intas',     type="string",  default=None,  help='specify "mcs" or "tes" for gcm P interpolation grid')
    parser.add_option('-N', '--no-api', action='store_true',dest='nocall',               default=False, help='do not recreate api file')

    ### GENERIC GRAPHICS SETTINGS
    parser.add_option('-c', '--color',  action='store',dest='clb',       type="string",  default=None,  help='change colormapS (also: nobar,onebar,def)')
    parser.add_option('--trycol',       action='store_true',dest='trycol',               default=False, help='try 9 typical color palette')
    parser.add_option('--nocolorb',     action='store_true',dest='nocolorb',             default=False, help='no color bar please')
    parser.add_option('--div',          action='store',dest='ndiv',      type="int",     default=10,    help='number of divisions in colorbar [10]')
    parser.add_option('--title',        action='store',dest='zetitle',   type="string",  default="fill",help='customize the whole title')
    parser.add_option('-T', '--tiled',  action='store_true',dest='tile',                 default=False, help='draw a tiled plot (3D) or add crosses (1D)')
    parser.add_option('--res',          action='store',dest='res',       type="float",   default=200.,  help='resolution for png outputs. --save png needed. [200.]')
    parser.add_option('--trans',        action='store',dest='trans',     type="float",   default=1.,    help='shaded plot transparency, 0 to 1 (=opaque) [1]')
    parser.add_option('--area',         action='store',dest='area',       type="string",   default=None,  help='area on the map to be plot [None]')
    parser.add_option('--xlabel',       action='store',dest='xlab',       type="string",  default=None, help='customize the x-axis label')
    parser.add_option('--ylabel',       action='store',dest='ylab',       type="string",  default=None, help='customize the y-axis label')
    parser.add_option('--labels',       action='store',dest='labels',     type="string",  default=None, help='customize 1D curve labels. Str comma-separated. [None]')
    parser.add_option('--lstyle',       action='store',dest='linestyle',  type="string",  default=None, help='customize 1D curve linestyles. Str comma-separ. [None]')

    ### SPECIFIC FOR MAPPING [MAPMODE 1]
    parser.add_option('-p',  '--proj',     action='store',dest='proj',      type="string",  default=None,  help='projection')
    parser.add_option('-b',  '--back',     action='store',dest='back',      type="string",  default=None,  help='background image [None]')
    parser.add_option('-W',  '--winds',    action='store_true',dest='winds',                default=False, help='wind vectors [False]')
    parser.add_option('--facwind',         action='store',dest='facwind',   type="float",   default=1,     help='wind vectors magnifying factor [1]')
    parser.add_option('-s',  '--stride',   action='store',dest='ste',       type="int",     default=3,     help='stride vectors [3]')
    parser.add_option('-z',  '--zoom',     action='store',dest='zoom',      type="float",   default=None,  help='zoom factor in %')
    parser.add_option('--blat',            action='store',dest='blat',      type="int",     default=None,  help='reference lat (or bounding lat for stere) [computed]')
    parser.add_option('--blon',            action='store',dest='blon',      type="int",     default=None,  help='reference lon [computed]')
    parser.add_option('--mark',            action='append',dest='mark',  type="string",  default=None, help='superimpose a crossmark at given lon,lat [None]')
    parser.add_option('--finddevil',       action='store_true',dest='mdevil',               default=False, help='superimpose a crossmark where the steepest dust devil is [False]')

    ### SPECIFIC FOR SLICING [MAPMODE 0]
    parser.add_option('--lat',          action='append',dest='slat',   type="string",  default=None, help='slices along lat. 2 comma-separated values: averaging')
    parser.add_option('--lon',          action='append',dest='slon',   type="string",  default=None, help='slices along lon. 2 comma-separated values: averaging')
    parser.add_option('--vert',         action='append',dest='svert',  type="string",  default=None, help='slices along vert. 2 comma-separated values: averaging') 
    parser.add_option('--column',       action='store_true',dest='column',             default=False,help='changes --vert z1,z2 from MEAN to INTEGRATE along z')
    parser.add_option('--time',         action='append',dest='stime',  type="string",  default=None, help='slices along time. 2 comma-separated values: averaging') 
    parser.add_option('--xmax',         action='store',dest='xmax',    type="float",   default=None, help='max value for x-axis in contour-plots [max(xaxis)]')
    parser.add_option('--ymax',         action='store',dest='ymax',    type="float",   default=None, help='max value for y-axis in contour-plots [max(yaxis)]')
    parser.add_option('--xmin',         action='store',dest='xmin',    type="float",   default=None, help='min value for x-axis in contour-plots [min(xaxis)]')
    parser.add_option('--ymin',         action='store',dest='ymin',    type="float",   default=None, help='min value for y-axis in contour-plots [min(yaxis)]')
    parser.add_option('--inverty',      action='store_true',dest='inverty',            default=False,help='force decreasing values along y-axis (e.g. p-levels)')
    parser.add_option('--logx',         action='store_true',dest='logx',               default=False,help='set x-axis to logarithmic')
    parser.add_option('--logy',         action='store_true',dest='logy',               default=False,help='set y-axis to logarithmic') 
    parser.add_option('--axtime',       action='store',dest='axtime',  type="string",  default=None, help='choose "ls","sol","lt" for time ref (1D or --time)')

    ### OPERATIONS BETWEEN FILES
    parser.add_option('--operation',    action='store',dest='operat',  type="string",  default=None,  help='operation to perform on input files given through -f. "+" or "-" acts on each input file by adding or substracting the ref file specified through --fref. "cat" acts on all input files in-a-row. "add_var" "sub_var" "mul_var" "div_var" acts on two variables (add _only to get only operation plot). "-_histo" will add an histogram plot for the "-" operation.')
    parser.add_option('--fref',         action='store',dest='fref',    type="string",  default=None,  help='reference namefile for the --operation option.')
    parser.add_option('--mope',         action='store',dest='vminope', type="float",   default=0.,  help='bounding minimum value for inter-file operation')
    parser.add_option('--Mope',         action='store',dest='vmaxope', type="float",   default=0.,  help='bounding maximum value for inter-file operation')
    parser.add_option('--titleref',     action='store',dest='titref',  type="string",  default="fill",  help='title for the reference file. [title of fig (1)]')

    ### SPECIAL
    parser.add_option('--tsat',         action='store_true',dest='tsat',               default=False,help='convert temperature field T in Tsat-T using pressure')
    parser.add_option('--stream',       action='store_true',dest='stream',             default=False,help='plot streamlines from streamfunction.e for vert/lat slices.')
    parser.add_option('--analysis',     action='store'     ,dest='analysis',           default=None ,help='Analysis operation. histo, density (kernel distribution estimate, with gaussian kernel only for the moment (many other distributions are available and can be added)), histodensity (overplot of both density and histo), fft. (currently fft works only on the z-axis, i.e. spatial Fourier Transform, and yields spectrum amplitude. To get enough bandwith, use API with a step of 10m on your data. Note that if you use --time A,B, the result will be the mean of FT at each timestep, and not the FT of the mean. The same apply to --lon and --lat, but does not apply to histo and density (for which arrays are flattened -> no mean).) [None]')

    return parser

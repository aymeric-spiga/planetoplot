#! /usr/bin/env python
##############################################
##  A MINIMAL PP.PY SCRIPT USING PPCLASS.PY ##
##  Author: A. Spiga 03/2013                ##
##############################################
from optparse import OptionParser ### TBR by argparse
from ppclass import pp, inspect
import sys
##############################################

# NB: this is a convenient command-line script
#     ... but ppclass is more versatile 
#     ... than what is proposed here
#     ... e.g. differences between files, 
#     ... complex operations,
#     ... see sample scripts

######################################
# define parser with version and usage 
######################################
parser = OptionParser()
parser.version = \
'''**************************************************
******** PLANETOPLOT (for help: pp.py -h) ********
**************************************************'''
parser.usage = \
'''pp.py [options] netCDF file(s)
(NB: no options --> simple inspection of variables and dimensions in netCDF files)
-------------------
PLANETOPLOT
--> command line tool to make nice & quick plots from netCDF files
--> based on python + numpy + scipy + matplotlib + basemap + netCDF4
--> Author: A. Spiga (LMD/UPMC) aymeric.spiga@upmc.fr
-------------------'''
parser.print_version()

########################################
# set options for the pp.py command line
########################################
parser.add_option('--verbose',action='store_true',dest='verbose',default=False,help='make the program verbose')
# field --> lower case
parser.add_option('-v','--var',action='append',dest='var',type="string",default=None,help="'variable' or ['var1','var2',etc]")
parser.add_option('-x','--lon',action='append',dest='x',type="string",default=None,help="x axis value. one value; or val1,val2 (computations)")
parser.add_option('-y','--lat',action='append',dest='y',type="string",default=None,help="y axis value. one value; or val1,val2 (computations)")
parser.add_option('-z','--vert',action='append',dest='z',type="string",default=None,help="z axis value. one value; or val1,val2 (computations)")
parser.add_option('-t','--time',action='append',dest='t',type="string",default=None,help="t axis value. one value; or val1,val2 (computations)")
parser.add_option('-u','--compute',action='store',dest='compute',type="string",default="mean",help="computation: mean, min, max, meanarea")
parser.add_option('-c','--contour',action='store',dest='contour',type="string",default=None,help="one 'variable' for contour")
parser.add_option('-i','--vecx',action='store',dest='vecx',type="string",default=None,help="one 'variable' for wind vector x component")
parser.add_option('-j','--vecy',action='store',dest='vecy',type="string",default=None,help="one 'variable' for wind vector y component")
parser.add_option('-m','--mult',action='store',dest='mult',type="float",default=None,help="multiplicative factor on field")
parser.add_option('-a','--add',action='store',dest='add',type="float",default=None,help="additive factor on field")
parser.add_option('-o','--output',action='store',dest='filename',type="string",default="myplot",help="name of output files")
parser.add_option('-d','--directory',action='store',dest='folder',type="string",default="./",help="directory of output files")
parser.add_option('-s','--changetime',action='store',dest='changetime',type="string",default=None,\
                  help="transformation on time axis : [None] | correctls | mars_sol2ls | mars_dayini | mars_meso_ls | mars_meso_sol | mars_meso_utc | mars_meso_lt ")
parser.add_option('-p','--print',action='store_true',dest='savtxt',default=False,help="[1D] output field+coord in an ASCII file")
parser.add_option('--stridex',action='store',dest='stridex',type="int",default=1,help="Load data every stridex grid points over x dimension")
parser.add_option('--stridey',action='store',dest='stridey',type="int",default=1,help="Load data every stridex grid points over y dimension")
parser.add_option('--stridez',action='store',dest='stridez',type="int",default=1,help="Load data every stridex grid points over z dimension")
parser.add_option('--stridet',action='store',dest='stridet',type="int",default=1,help="Load data every stridex grid points over t dimension")
# plot --> upper case
# -- generic
parser.add_option('-T','--title',action='append',dest='title',type="string",default=None,help="change 'title'")
parser.add_option('-X','--xlabel',action='append',dest='xlabel',type="string",default=None,help="change 'xlabel'")
parser.add_option('-Y','--ylabel',action='append',dest='ylabel',type="string",default=None,help="change 'ylabel'")
parser.add_option('-D','--div',action='store',dest='div',type="int",default=20,help="integer for number of divisions")
parser.add_option('-H','--trans',action='store',dest='trans',type="float",default=1.0,help="float for transparency (0=transp,1=opaque)")
parser.add_option('-Z','--logy',action='store_true',dest='logy',default=False,help="set log for vertical axis")
parser.add_option('-O','--save',action='store',dest='out',type="string",default="gui",help="save mode: 'gui' 'png' 'pdf' 'eps' 'svg' 'ps'")
parser.add_option('-V','--void',action='store_true',dest='void',default=False,help="no colorbar, no title, no labels")
parser.add_option('-U','--units',action='append',dest='units',type="string",default=None,help="units for the field")
# -- 1D plot
parser.add_option('-L','--lstyle',action='append',dest='lstyle',type="string",default=None,help="[1D] linestyle: '-' '--' '.' '..'")
parser.add_option('-Q','--color',action='append',dest='color',type="string",default=None,help="[1D] color: 'b' 'g' 'r' etc")
parser.add_option('-K','--marker',action='append',dest='marker',type="string",default=None,help="[1D] marker: '' 'x' 'o' etc")
parser.add_option('-S','--superpose',action='store_true',dest='superpose',default=False,help="[1D] use same axis for all plots")
parser.add_option('-E','--label',action='append',dest='label',type="string",default=None,help="[1D] label for line")
parser.add_option('--xcoeff',action='append',dest='xcoeff',type="float",default=None,help="[1D] multiply x axis")
parser.add_option('--ycoeff',action='append',dest='ycoeff',type="float",default=None,help="[1D] multiply y axis")
parser.add_option('--xmin',action='append',dest='xmin',type="float",default=None,help="[1D] min bound x axis")
parser.add_option('--ymin',action='append',dest='ymin',type="float",default=None,help="[1D] min bound y axis")
parser.add_option('--xmax',action='append',dest='xmax',type="float",default=None,help="[1D] max bound x axis")
parser.add_option('--ymax',action='append',dest='ymax',type="float",default=None,help="[1D] max bound y axis")
parser.add_option('--modx',action='append',dest='modx',type="float",default=None,help="[1D] change xticks with a modulo")
# -- 2D plot
parser.add_option('-C','--colorb',action='append',dest='colorb',type="string",default=None,help="[2D] colormap: http://micropore.files.wordpress.com/2010/06/colormaps.png")
parser.add_option('-P','--proj',action='append',dest='proj',type="string",default=None,help="[2D] map projection: 'cyl' 'npstere' 'spstere' 'ortho' 'moll' 'robin' 'lcc' 'laea' 'merc' 'noproj'")
parser.add_option('-B','--back',action='append',dest='back',type="string",default=None,help='[2D] predefined map background (cf. set_back.txt)')
parser.add_option('-A','--area',action='append',dest='area',type="string",default=None,help='[2D] predefined region of mapping (cf. set_area.txt)')
parser.add_option('-I','--blon',action='append',dest='blon',type="float",default=None,help='[2D] float: bounding longitude for stere (or center longitude for ortho)')
parser.add_option('-J','--blat',action='append',dest='blat',type="float",default=None,help='[2D] float: bounding latitude for stere (or center latitude for ortho) ')
parser.add_option('-N','--vmin',action='append',dest='vmin',type="float",default=None,help='[2D] float: minimum value for displayed field')
parser.add_option('-M','--vmax',action='append',dest='vmax',type="float",default=None,help='[2D] float: maximum value for displayed field')
parser.add_option('-W','--wscale',action='append',dest='wscale',type="float",default=None,help='[2D] float: set size of reference wind vector')
parser.add_option('--stridevecx',action='store',dest='stridevecx',type="int",default=1,help="Define an abscissa stride on vectors only -- not on field")
parser.add_option('--stridevecy',action='store',dest='stridevecy',type="int",default=1,help="Define an ordinate stride on vectors only -- not on field")
###########################
(opt,args) = parser.parse_args()
# remains F G R  

######################################
# get arguments (one or several files)
######################################
if args is None:
    print "Stop here! I need file(s) as argument(s)!" ; exit()
else:
    files = args

#############################################
# a possibility to simply inspect the file(s)
#############################################
if opt.var is None:
    for filename in files: inspect(filename)
    exit()

######################################
# use ppclass to get field and plot it
######################################
# treat the case of additional vectors or contours (contours must be before vectors)
var = [] ; vargoal = []
for element in opt.var:
    var.append(element) ; vargoal.append("main")
    if opt.contour is not None: var.append(opt.contour) ; vargoal.append("contour")
    if opt.vecx is not None: var.append(opt.vecx) ; vargoal.append("vector")
    if opt.vecy is not None: var.append(opt.vecy) ; vargoal.append("vector")
# set pp object
user = pp()
user.file = files
user.var = var ; user.vargoal = vargoal
user.x = opt.x ; user.y = opt.y 
user.z = opt.z ; user.t = opt.t
user.verbose = opt.verbose
user.compute = opt.compute
user.changetime = opt.changetime
user.stridex = opt.stridex
user.stridey = opt.stridey
user.stridez = opt.stridez
user.stridet = opt.stridet
user.stridevecx = opt.stridevecx
user.stridevecy = opt.stridevecy
user.savtxt = opt.savtxt
# define field
user.define()
# retrieve field
user.retrieve()
# some possible operations
if opt.add is not None: user = user + opt.add
if opt.mult is not None: user = user * opt.mult
# get some options
user.superpose = opt.superpose
user.filename = opt.filename
user.folder = opt.folder
user.out = opt.out
# if noproj is given for proj, no map mode
if opt.proj is not None:
    if 'noproj' in opt.proj: 
        user.noproj = True
# if user wants to give a name, we drop the indication of date
if opt.filename != "myplot":
    user.includedate = False
# define plot
user.defineplot()
# user-defined plot settings
# ... shouldn't this be before defineplot?
user.getopt(opt)
# make plot
user.makeplot()

####################################
# save a .sh file with the command #
####################################
command = ""
for arg in sys.argv: command = command + arg + ' '
try:
    f = open(opt.folder+'/'+opt.filename+'.sh', 'w')
    f.write(command)	
except IOError:
    print "!! WARNING !! pp.py command not saved. Probably do not have permission to write here."

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
parser.add_option('-o','--output',action='store',dest='filename',type="string",default=None,help="name of output files")
parser.add_option('-d','--directory',action='store',dest='folder',type="string",default="./",help="directory of output files")
parser.add_option('-s','--changetime',action='store',dest='changetime',type="string",default=None,\
                  help="transformation on time axis : [None] | correctls | mars_sol2ls | mars_dayini | mars_meso_ls | mars_meso_sol | mars_meso_utc | mars_meso_lt ")
parser.add_option('-p','--print',action='store_true',dest='savtxt',default=False,help="[1D] output field+coord in an ASCII file")
parser.add_option('--sx',action='store',dest='sx',type="int",default=1,help="Load data every sx grid points over x dimension")
parser.add_option('--sy',action='store',dest='sy',type="int",default=1,help="Load data every sy grid points over y dimension")
parser.add_option('--sz',action='store',dest='sz',type="int",default=1,help="Load data every sz grid points over z dimension")
parser.add_option('--st',action='store',dest='st',type="int",default=1,help="Load data every st grid points over t dimension")
# plot --> upper case
# -- generic
parser.add_option('-T','--title',action='append',dest='title',type="string",default=None,help="change 'title'")
parser.add_option('-X','--xlabel',action='append',dest='xlabel',type="string",default=None,help="change 'xlabel'")
parser.add_option('-Y','--ylabel',action='append',dest='ylabel',type="string",default=None,help="change 'ylabel'")
parser.add_option('-D','--div',action='store',dest='div',type="int",default=20,help="integer for number of divisions")
parser.add_option('-H','--trans',action='store',dest='trans',type="float",default=1.0,help="float for transparency (0=transp,1=opaque)")
parser.add_option('-Z','--logy',action='store_true',dest='logy',default=False,help="set log for vertical axis")
parser.add_option('-O','--save',action='store',dest='out',type="string",default=None,help="save mode: 'gui' 'png' 'pdf' 'eps' 'svg' 'ps'")
parser.add_option('-V','--void',action='store_true',dest='void',default=False,help="no colorbar, no title, no labels")
parser.add_option('-U','--units',action='append',dest='units',type="string",default=None,help="units for the field")
parser.add_option('-F','--fmt',action='append',dest='fmt',type="string",default=None,help="values formatting. ex: '%.0f' '%3.1e'")
parser.add_option('--xcoeff',action='append',dest='xcoeff',type="float",default=None,help="multiply x axis [not for 2D map]")
parser.add_option('--ycoeff',action='append',dest='ycoeff',type="float",default=None,help="multiply y axis [not for 2D map]")
parser.add_option('--xmin',action='append',dest='xmin',type="float",default=None,help="min bound x axis [not for 2D map]")
parser.add_option('--ymin',action='append',dest='ymin',type="float",default=None,help="min bound y axis [not for 2D map]")
parser.add_option('--xmax',action='append',dest='xmax',type="float",default=None,help="max bound x axis [not for 2D map]")
parser.add_option('--ymax',action='append',dest='ymax',type="float",default=None,help="max bound y axis [not for 2D map]")
parser.add_option('--nxticks',action='append',dest='nxticks',type="float",default=None,help="ticks for x axis [not for 2D map]")
parser.add_option('--nyticks',action='append',dest='nyticks',type="float",default=None,help="ticks for y axis [not for 2D map]")
parser.add_option('--xp',action='store',dest='xp',type="int",default=None,help="x size of figure (integer)")
parser.add_option('--yp',action='store',dest='yp',type="int",default=None,help="y size of figure (integer)")
# -- 1D plot
parser.add_option('-L','--linestyle',action='append',dest='linestyle',type="string",default=None,help="[1D] linestyle: '-' '--' '.' '..'")
parser.add_option('-Q','--color',action='append',dest='color',type="string",default=None,help="[1D] color: 'b' 'g' 'r' etc")
parser.add_option('-K','--marker',action='append',dest='marker',type="string",default=None,help="[1D] marker: '' 'x' 'o' etc")
parser.add_option('-S','--superpose',action='store_true',dest='superpose',default=False,help="[1D] use same axis for all plots")
parser.add_option('-E','--legend',action='append',dest='legend',type="string",default=None,help="[1D] legend for line")
parser.add_option('--modx',action='append',dest='modx',type="float",default=None,help="[1D] change xticks with a modulo")
# -- 2D plot
parser.add_option('-C','--colorbar',action='append',dest='colorbar',type="string",default=None,help="[2D] colormap: http://micropore.files.wordpress.com/2010/06/colormaps.png")
parser.add_option('-P','--proj',action='store',dest='proj',type="string",default=None,help="[2D] map projection: 'cyl' 'npstere' 'spstere' 'ortho' 'moll' 'robin' 'lcc' 'laea' 'merc'")
parser.add_option('-B','--back',action='append',dest='back',type="string",default=None,help='[2D] predefined map background (cf. set_back.txt)')
parser.add_option('-A','--area',action='append',dest='area',type="string",default=None,help='[2D] predefined region of mapping (cf. set_area.txt)')
parser.add_option('-I','--blon',action='append',dest='blon',type="float",default=None,help='[2D] float: bounding longitude for stere (or center longitude for ortho)')
parser.add_option('-J','--blat',action='append',dest='blat',type="float",default=None,help='[2D] float: bounding latitude for stere (or center latitude for ortho) ')
parser.add_option('-N','--vmin',action='append',dest='vmin',type="float",default=None,help='[2D] float: minimum value for displayed field')
parser.add_option('-M','--vmax',action='append',dest='vmax',type="float",default=None,help='[2D] float: maximum value for displayed field')
parser.add_option('-W','--wscale',action='append',dest='wscale',type="float",default=None,help='[2D] float: set size of reference wind vector')
parser.add_option('--svx',action='store',dest='svx',type="int",default=None,help="Define an abscissa stride on vectors only -- not on field")
parser.add_option('--svy',action='store',dest='svy',type="int",default=None,help="Define an ordinate stride on vectors only -- not on field")
parser.add_option('--cbticks',action='append',dest='cbticks',type="float",default=None,help="ticks for colorbar")
###########################
(opt,args) = parser.parse_args()
# remains G R  

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
user.sx = opt.sx ; user.sy = opt.sy
user.sz = opt.sz ; user.st = opt.st
user.svx = opt.svx ; user.svy = opt.svy
user.savtxt = opt.savtxt
if opt.xp is not None: user.xp = opt.xp
if opt.yp is not None: user.yp = opt.yp
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
user.proj = opt.proj
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
if opt.filename is not None:
  try:
    f = open(opt.folder+'/'+opt.filename+'.sh', 'w')
    f.write(command)	
  except IOError:
    print "!! WARNING !! pp.py command not saved. Probably do not have permission to write here."

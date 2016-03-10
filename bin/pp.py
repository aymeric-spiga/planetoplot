#! /usr/bin/env python
##############################################
##  A MINIMAL PP.PY SCRIPT USING PPCLASS.PY ##
##  Author: A. Spiga 03/2013                ##
##############################################
from optparse import OptionParser ### TBR by argparse
from ppclass import pp, inspect
import ppplot
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
parser.add_option('--useindex',action='store_true',dest="useindex",default=False,help="Use index for arrays and not values of dimensions")
parser.add_option('--kind3d',action='store',dest='kind3d',type="string",default="tyx",help="dimensions if rank<4: tzy, tyx (default)")
# plot options --> upper case. see ppplot.
parser = ppplot.opt(parser)
parser = ppplot.opt1d(parser)
parser = ppplot.opt2d(parser)
###########################
(opt,args) = parser.parse_args()
# remains G R  
if (len(args) == 0):
    parser.print_version()

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
if not user.verbose: user.quiet = True
user.compute = opt.compute
user.changetime = opt.changetime
user.useindex = opt.useindex
user.sx = opt.sx ; user.sy = opt.sy
user.sz = opt.sz ; user.st = opt.st
user.svx = opt.svx ; user.svy = opt.svy
user.savtxt = opt.savtxt
user.kind3d = opt.kind3d
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
user.res = opt.res
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

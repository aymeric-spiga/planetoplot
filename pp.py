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

#####################################
# get arguments from the command line
#####################################
parser = OptionParser()
parser.add_option('--verbose',action='store_true',dest='verbose',default=False,help='')
# field --> lower case
parser.add_option('-f','--file',action='append',dest='file',type="string",default=None,help='')
parser.add_option('-v','--var',action='append',dest='var',type="string",default=None,help='')
parser.add_option('-x','--lon',action='append',dest='x',type="string",default=None,help='')
parser.add_option('-y','--lat',action='append',dest='y',type="string",default=None,help='')
parser.add_option('-z','--vert',action='append',dest='z',type="string",default=None,help='')
parser.add_option('-t','--time',action='append',dest='t',type="string",default=None,help='')
parser.add_option('-c','--contour',action='store',dest='contour',type="string",default=None,help='')
parser.add_option('-i','--vecx',action='store',dest='vecx',type="string",default=None,help='')
parser.add_option('-j','--vecy',action='store',dest='vecy',type="string",default=None,help='')
parser.add_option('-m','--mult',action='store',dest='mult',type="float",default=None,help='')
parser.add_option('-a','--add',action='store',dest='add',type="float",default=None,help='')
parser.add_option('-o','--output',action='store',dest='filename',type="string",default="myplot",help='')
parser.add_option('-d','--directory',action='store',dest='folder',type="string",default="./",help='')
# plot --> upper case
# -- generic
parser.add_option('-T','--title',action='append',dest='title',type="string",default=None,help='')
parser.add_option('-X','--xlabel',action='append',dest='xlabel',type="string",default=None,help='')
parser.add_option('-Y','--ylabel',action='append',dest='ylabel',type="string",default=None,help='')
parser.add_option('-D','--div',action='store',dest='div',type="int",default=20,help='')
parser.add_option('-H','--trans',action='store',dest='trans',type="float",default=1.0,help='')
parser.add_option('-Z','--logy',action='store_true',dest='logy',default=False,help='')
parser.add_option('-O','--save',action='store',dest='out',type="string",default="gui",help='')
# -- 1D plot
parser.add_option('-L','--lstyle',action='append',dest='lstyle',type="string",default=None,help='')
parser.add_option('-S','--superpose',action='store_true',dest='superpose',default=False,help='')
# -- 2D plot
parser.add_option('-C','--colorb',action='append',dest='colorb',type="string",default=None,help='')
parser.add_option('-P','--proj',action='append',dest='proj',type="string",default=None,help='')
parser.add_option('-B','--back',action='append',dest='back',type="string",default=None,help='')
parser.add_option('-A','--area',action='append',dest='area',type="string",default=None,help='')
parser.add_option('-I','--blon',action='append',dest='blon',type="float",default=None,help='')
parser.add_option('-J','--blat',action='append',dest='blat',type="float",default=None,help='')
parser.add_option('-N','--vmin',action='append',dest='vmin',type="float",default=None,help='')
parser.add_option('-M','--vmax',action='append',dest='vmax',type="float",default=None,help='')
(opt,args) = parser.parse_args()

##########################################
# a possibility to simply inspect the file
##########################################
if opt.file is None:
    print "Stop here. I need at least a file: -f FILE" ; exit()
if opt.var is None:
    for filename in opt.file: inspect(filename)
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
user = pp(file=opt.file, var=var, vargoal=vargoal, \
          x=opt.x, y=opt.y, z=opt.z, t=opt.t,\
          verbose=opt.verbose)
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
# define plot
user.defineplot()
# user-defined plot settings
user.getopt(opt)
# make plot
user.makeplot()

####################################
# save a .sh file with the command #
####################################
command = ""
for arg in sys.argv: command = command + arg + ' '
f = open(opt.filename+'.sh', 'w')
f.write(command)

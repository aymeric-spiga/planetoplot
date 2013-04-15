#! /usr/bin/env python
from api_wrapper import api_onelevel
from ppclass import separatenames,inspect
from optparse import OptionParser ### TBR by argparse
import numpy as np
# -----------------------------------------------------------------------
# A simple script to call the API program to interpolate mesoscale files
# ... NB: API_WRAPPER necessary with api compiled with f2py
# Author: A. Spiga 04/2013
# -----------------------------------------------------------------------

# default settings
# ---------------- 
defitp = 4
deflvl = "0.1"

# define parser with version and usage 
# ------------------------------------
parser = OptionParser()
parser.version = \
'''***********************************************
******** EXAPI (for help: exapi.py -h) ********
***********************************************'''
parser.usage = \
'''exapi.py [options] netCDF file(s)'''
parser.print_version()

# define options and get options+arguments
# ----------------------------------------
parser.add_option('-v','--var',action='store',dest='var',type="string",default=None,help='Variables processed. Not used means all variables.')
parser.add_option('-i','--interp',action='store',dest='itp',type="int",default=defitp,help='Mode: 2=press /  3=z-amr / 4=z-als (default)')
parser.add_option('-l','--level',action='store',dest='lvl',type="string",default=deflvl,help='Levels: start[,stop,step] (-i 2: Pa)(-i 3,4: km)')
parser.add_option('-o','--output',action='store',dest='output',type="string",default=None,help="name of output files")
parser.add_option('-d','--directory',action='store',dest='folder',type="string",default=None,help="directory of output files")
(opt,args) = parser.parse_args()

# get variables
# if None, take all!
# ------------------
if opt.var is None: 
  ze_process = "all"
  zevars = None
else: 
  ze_process = "list"
  # (we start by unravelling user input in an array)
  zevars = separatenames(opt.var)
  # (we help a little the user about naming certain variables)
  for i in range(len(zevars)):
    if zevars[i] in ['t','temp','temperature']: zevars[i] = 'tk'
    elif zevars[i] in ['T','temppot','theta']: zevars[i] = 'tpot'
    elif zevars[i] in ['u','v','U','V','Um','Vm','uv','UV','wind']: zevars[i] = 'uvmet'
  # (we recombine them all for call to api)
  list = ""
  for el in zevars: list = list + el + ","
  zevars = list

# get the kind of interpolation
# -----------------------------
inputnvert = separatenames(opt.lvl)

# prepare levels: one-level only
# ------------------------------
if len(inputnvert) == 1:
    zelevel = float(inputnvert[0])
    ze_interp_levels = [-9999.]

# prepare levels: several levels
# ------------------------------
elif len(inputnvert) > 1:
    # initialize. make number of interp levels to 20 if not given.
    zelevel = -99.
    start = float(inputnvert[0])
    stop = float(inputnvert[1])
    if len(inputnvert) == 2:  numsample = 20
    else:                     numsample = float(inputnvert[2])
    # make the interval. either normal -- or log if pressure.
    if stop > start:
        # altitude coordinates
        ze_interp_levels = np.linspace(start,stop,numsample)
    else:
        # pressure coordinates
        ze_interp_levels = np.logspace(np.log10(start),np.log10(stop),numsample)

# main api call
# -------------
for file in args:
    print "EXAPI: working on file",file 
    newname = api_onelevel ( \
        path_to_input = '', \
        path_to_output = None, \
        input_name = file, \
        output_name = opt.output, \
        fields = zevars, \
        interp_method = opt.itp, \
        interp_level = ze_interp_levels, \
        onelevel = zelevel, \
        process = ze_process )
    print "EXAPI: inspect the new file."
    inspect(newname)
    print "EXAPI: OK. done with this file."
    



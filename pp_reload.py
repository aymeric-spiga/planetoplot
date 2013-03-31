#! /usr/bin/env python
# -----------------------------------------------------------------------
# A very simple script to reload a previously saved plot array in ppclass
# ... for a demo, try "pp_reload.py demo_data/*"
# Author: A. Spiga 03/2013
# -----------------------------------------------------------------------
from optparse import OptionParser ### TBR by argparse
from ppclass import pp
# -----------------------------------------------------------------------
parser = OptionParser()
parser.add_option('-O','--out',action='store',dest='out',type="string",default="gui",help='Specify a new output format')
parser.add_option('-K','--marker',action='store',dest='marker',type="string",default=None,help="[1D] Define a new marker")
parser.add_option('-P','--proj',action='store',dest='proj',type="string",default=None,help='[2D] Define a new map projection')
parser.add_option('-C','--colorb',action='store',dest='colorb',type="string",default=None,help="[2D] Define a new colormap")
(opt,args) = parser.parse_args()
# -----------------------------------------------------------------------
for files in args:
    yeah = pp()
    yeah.defineplot(loadfile=files)
    yeah.out = opt.out
    if opt.proj is not None:
      for plot in yeah.p:
        plot.proj = opt.proj
    if opt.colorb is not None:
      for plot in yeah.p:
        plot.colorb = opt.colorb
    if opt.marker is not None:
      for plot in yeah.p:
        plot.marker = opt.marker
    yeah.makeplot()


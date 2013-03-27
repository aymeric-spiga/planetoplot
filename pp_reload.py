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
parser.add_option('-P','--proj',action='store',dest='proj',type="string",default=None,help='Define a new map projection')
(opt,args) = parser.parse_args()
# -----------------------------------------------------------------------
for files in args:
    yeah = pp()
    yeah.defineplot(loadfile=files)
    yeah.out = opt.out
    if opt.proj is not None:
      for plot in yeah.p:
        plot.proj = opt.proj
    yeah.makeplot()


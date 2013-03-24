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
parser.add_option('-o',action='store',dest='out',type="string",default="gui",help='Specify a new output format')
(opt,args) = parser.parse_args()
# -----------------------------------------------------------------------
for files in args:
    yeah = pp()
    yeah.defineplot(loadfile=files)
    yeah.out = opt.out
    yeah.makeplot()


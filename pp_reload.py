#! /usr/bin/env python
# -----------------------------------------------------------------------
# A very simple script to reload a previously saved plot array in ppclass
# ... for a demo, try "pp_reload.py demo_data/*"
# Author: A. Spiga 03/2013
# -----------------------------------------------------------------------
from optparse import OptionParser ### TBR by argparse
from ppclass import pp
from ppplot import rainbow
# -----------------------------------------------------------------------
parser = OptionParser()
parser.add_option('-O','--out',action='store',dest='out',type="string",default="gui",help='Specify a new output format')
parser.add_option('-T','--title',action='store',dest='title',type="string",default=None,help='Give a new title')
parser.add_option('-K','--marker',action='store',dest='marker',type="string",default=None,help="[1D] Define a new marker")
parser.add_option('-P','--proj',action='store',dest='proj',type="string",default=None,help='[2D] Define a new map projection')
parser.add_option('-C','--colorb',action='store',dest='colorb',type="string",default=None,help="[2D] Define a new colormap")
parser.add_option('-E','--explore',action='store_true',dest='explore',default=False,help="[2D] Try many colormaps")
(opt,args) = parser.parse_args()
# -----------------------------------------------------------------------
clb = ["Greys","Blues","YlOrRd",\
       "jet","spectral","hot",\
       "RdBu","RdYlBu","Paired",\
       "gist_ncar","gist_rainbow","gist_stern"]
# -----------------------------------------------------------------------

for files in args:

    yeah = pp()
    yeah.quiet = True
    yeah.defineplot(loadfile=files)
    yeah.out = opt.out

    if opt.proj is not None:
      for plot in yeah.p:
        plot.proj = opt.proj
    if opt.colorb is not None:
      yeah.colorb = opt.colorb
      for plot in yeah.p:
        plot.colorb = opt.colorb
    if opt.marker is not None:
      for plot in yeah.p:
        plot.marker = opt.marker
    if opt.title is not None:
      for plot in yeah.p:
        plot.title = opt.title

    if opt.explore:
      for cm in clb:
       for plot in yeah.p:
         plot.colorb = cm
         plot.title = cm
       yeah.makeplot()
    else:         
      yeah.makeplot()


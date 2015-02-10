#! /usr/bin/env python
import ppplot
import numpy as np
from optparse import OptionParser ### TBR by argparse

# inputs and options
parser = OptionParser()
parser.usage = \
'''
asciiplot.py [options] text file(s)
-- default is col2 for field and col1 for coord
-- this can be set with -x and -y options
   (or use --swap if applicable)
-- one-column files are also supported
'''
parser.add_option('-x','--colx',action='store',dest='colx',type="int",default=1,help='column for x axis')
parser.add_option('-y','--coly',action='store',dest='coly',type="int",default=2,help='column for y axis')
parser.add_option('-s','--skip',action='store',dest='skiprows',type="int",default=0,help='skip first rows in file(s)')
parser = ppplot.opt(parser) # common options for plots
parser = ppplot.opt1d(parser) # common options for plots
(opt,args) = parser.parse_args()

# plot object + options
pl = ppplot.plot1d()

# for all input files
count = 0
for fff in args:

  # transfer options to plot object
  pl.transopt(opt,num=count)

  # load data
  var = np.transpose(np.loadtxt(fff,skiprows=opt.skiprows))

  # get coord
  if len(var.shape) == 1:
    pl.f = var
    pl.x = None # important for chained plots
  else:
    pl.f = var[opt.coly-1]
    pl.x = var[opt.colx-1]

  # make plot
  pl.make()
  count = count + 1

# show plot
ppplot.show()

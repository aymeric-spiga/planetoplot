#! /usr/bin/env python
from ppclass import pp

# define object
p = pp()

# define file, variable, slice
p.file = "/home/aymeric/Big_Data/DATAPLOT/diagfired.nc"
p.var = "icetot"
p.t = 0.5

# get field. convert from kg m-2 to pr_mic
p.get()
p = 1000.*p

# customize plot
# -- map projection, background image
p.proj = "robin"
p.back = "vishires"
p.trans = 0.5
# -- title, colorbar
p.title = 'Water ice clouds on Mars'
p.units = "pr-$\mu$m"
p.colorbar = "spectral"
p.vmin = 0.
p.vmax = 8.
p.div = 32
p.fmt = '%.1f'
# -- save info
p.filename = "weblowres"
p.includedate = False
p.out = "png"
p.res = 30

# generate plot
p.plot()

# generate a higher resolution version
p.filename = "webhires"
p.res = 150
p.plot()

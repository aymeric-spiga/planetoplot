#! /usr/bin/env python
from ppclass import pp

u = pp()
u.file = "/home/aymeric/Big_Data/DATAPLOT/diagfired.nc"
u.var = ["temp","phisinit","u","v"]
u.vargoal = ["main","contour","vector","vector"]
u.t = "0.5,0.8"
u.z = "50000"
u.filename = "vector"

# stride both x and y
# this impacts field + vector
u.sx = 3
u.sy = 3
u.getplot()
u.sx = 1 # (reinitialise)
u.sy = 1 # (reinitialise)

# stride vectors only
# not field (here topography)
u.sx = 3
u.sy = 3
u.getplot()

u.z = "50"
u.filename = "myplot"
u.getplot()

u.colorbar = "jet"
u.trans = 0.0
u.back = "vis"

u.makeplot()

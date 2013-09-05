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
u.stridex = 3
u.stridey = 3
u.getplot()
u.stridex = 1 # (reinitialise)
u.stridey = 1 # (reinitialise)

# stride vectors only
# not field (here topography)
u.stridevecx = 3
u.stridevecy = 3
u.getplot()



u.z = "50"
u.filename = "myplot"
u.getplot()


u.p[0].colorb = "jet"
u.p[0].trans = 0.0
u.p[0].back = "vis"

u.makeplot()

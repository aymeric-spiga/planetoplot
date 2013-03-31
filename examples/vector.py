#! /usr/bin/env python
from ppclass import pp

u = pp()
u.file = "/home/aymeric/Big_Data/DATAPLOT/diagfired.nc"
u.var = ["temp","phisinit","u","v"]
u.vargoal = ["main","contour","vector","vector"]
u.t = "0.5,0.8"
u.z = "50000"
u.filename = "vector"
u.get()
u.defineplot()
#u.p[0].proj = "ortho"
#u.p[0].blat = -45.
u.makeplot()


u.z = "50"
u.filename = "myplot"
u.getplot()


u.p[0].colorb = "jet"
u.p[0].trans = 0.0
u.p[0].back = "vis"

u.makeplot()

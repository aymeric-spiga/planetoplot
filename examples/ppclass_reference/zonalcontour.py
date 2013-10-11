#! /usr/bin/env python
from ppclass import pp

u = pp()
u.file = "/home/aymeric/Big_Data/DATAPLOT/diagfired.nc"
u.var = ["temp","u"]
u.vargoal = ["main","contour"]
u.x = "0"
u.t = "0.5"
u.filename = "zonalcontour"
u.get()
u.defineplot()
u.div = 30.
u.colorbar = "spectral"
u.makeplot()

u.var = ["u","u"]
u.vargoal = ["main","contour"]
u.get()
u.defineplot()
u.div = 30.
u.makeplot()

#! /usr/bin/env python
from ppclass import pp

u = pp()
u.file = "/home/aymeric/Big_Data/DATAPLOT/diagfired.nc"
u.var = ["temp","u"]
u.vargoal = ["main","contour"]
u.x = "0"
u.t = "0.5"
u.get()
u.defineplot()
u.p[0].div = 30.
u.p[0].colorb = "spectral"
u.makeplot()


u.var = ["u","u"]
u.vargoal = ["main","contour"]
u.get()
u.defineplot()
u.p[0].div = 30.
u.makeplot()

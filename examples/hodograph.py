#! /usr/bin/env python
from ppclass import pp

u = pp()
u.file = "/home/aymeric/Big_Data/DATAPLOT/diagfired.nc"
u.var = "u"
u.x = 0.
u.y = 0.
u.z = 10.
u.get()

v = pp()
v << u
v.var = "v"
v.get()

# u as a function of v
hodo = u.f(v)
hodo.filename = "hodograph"
hodo.makeplot()

# v as a function of u
hodo2 = v.f(u)
hodo2.makeplot()

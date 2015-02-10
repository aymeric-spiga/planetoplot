#! /usr/bin/env python
from ppclass import pp

m = pp()
m.file = "/home/aymeric/Big_Data/POLAR_APPERE_lowres_wrfout_d01_2024-03-04_06.nc"
m.var = ["tk","tk"]

m.y = ["80","82"]
m.t = ["2","5"]

m.getplot()


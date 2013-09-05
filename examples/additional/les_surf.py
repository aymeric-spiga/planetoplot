#! /usr/bin/env python
from ppclass import pp

for t in range(3500,3600):
  les = pp()
  #les.file = "/home/aymeric/Big_Data/press_ustm_exomars.nc"
  les.file = "/home/aymeric/Big_Data/uv.nc"
  les.var = ["U","V"]
  les.t = t
  les.xcoeff = 100./1000.
  les.ycoeff = 100./1000.
  les.xlabel = "x distance (km)"
  les.ylabel = "y distance (km)"
  les.vmin = -8.
  les.vmax = 8.
  les.div = 32
  #les.out = "png"
  les.res = 75
  les.get()
  les.defineplot()
  #les.p[0].vmin = 715.5
  #les.p[0].vmax = 716.5
  #les.p[1].vmin = 0.0
  #les.p[1].vmax = 1.0
  les.p[0].title = "East-West wind"
  les.p[1].title = "North-South wind"
  les.makeplot()



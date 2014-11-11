#! /usr/bin/env python
from ppclass import pp


opac = pp()
opac.file = "../tau.nc"
opac.var = "TAU"
opac.contour = "HGT"
#opac.changetime = "mars_meso_lt"
opac.vmin = 0.
opac.vmax = 4.
opac.back = "vishires"
opac.title = ""
opac.proj = "ortho"
opac.out = "png"
opac.res = 75
opac.nopickle = True

count = 0
for ttt in range(0,135,3):
  opac.filename = "taustorm"+"%03d"%(count)
  opac.includedate = False
  opac.t = ttt
  opac.getplot()
  count = count+1

#! /usr/bin/env python
from ppclass import pp
import matplotlib.pyplot as mpl
import numpy as np
import ppplot

u = pp()
u.file = "/home/aymeric/Big_Data/ustar.101x101x201.CaseA.w30_zipbl.nc"
u.var = "USTM"
u.x = "0,1000"
u.y = "0,1000"

for yeah in ["300","200,300"]:
  u.t = yeah
  u.compute = "nothing"
  ustm = u.getf()
  u.compute = "max" ; zemax = u.getf()
  u.compute = "min" ; zemin = u.getf()
  u.compute = "mean" ; zemean = u.getf()
  ppplot.figuref(x=4,y=4)
  mpl.hist(np.ravel(ustm),np.arange(zemin,zemax,0.1),log=True)
  mpl.title("$\mu$=%.2f / m=%.2f / M=%.2f" % (zemean,zemin,zemax))
  mpl.xlabel('Friction velocity $u_{\star}$ (m s$^{-1}$)')
  ppplot.save(mode="png",filename="hist")
  ppplot.close()



u.x = None
u.y = None
u.t = "0,1000"
u.compute = "max"
u.xcoeff = 0.1
u.ycoeff = 0.1
u.xlabel = "x (km)"
u.ylabel = "y (km)"
u.title = 'Daytime maximum $u\star$'
u.vmin = 1.5
u.vmax = 2.5
u.fmt = "%.1f"
u.xp = 8
u.yp = 8
u.filename = "maxustm"
u.includedate = False
u.out = "png"
u.getplot()


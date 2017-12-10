#! /usr/bin/env python
from ppclass import pp
import matplotlib.pyplot as mpl
import numpy as np
import ppplot

u = pp()
#u.file = "/home/aymeric/Big_Data/ustar.101x101x201.CaseA.w30_zipbl.nc"
u.file = "BIGLES10m_wind5_USTM_9-11.nc"
u.var = "USTM"
u.x = "0,1000"
u.y = "0,1000"

tttall = "0,1e10"

for yeah in [tttall]:
#for yeah in ["0"]:
  u.t = yeah
  u.compute = "nothing"
  ustm = u.getf()
  u.compute = "max" ; zemax = u.getf()
  u.compute = "min" ; zemin = u.getf()
  u.compute = "mean" ; zemean = u.getf()
  ppplot.figuref(x=4,y=4)
  dval = 0.05
  bins = np.arange(zemin,zemax,dval)
  hh = mpl.hist(np.ravel(ustm),bins,log=True)
  print hh
  mpl.title("$\mu$=%.2f / m=%.2f / M=%.2f" % (zemean,zemin,zemax))
  mpl.xlabel('Friction velocity $u_{\star}$ (m s$^{-1}$)')
  ppplot.save(mode="png",filename="roughness_hist")
  ppplot.close()


u.x = None
u.y = None
u.t = tttall
u.compute = "max"
u.xcoeff = 0.01
u.ycoeff = 0.01
u.xlabel = "x (km)"
u.ylabel = "y (km)"
u.title = 'maximum $u\star$'
u.vmin = 0.4
u.vmax = 1.1
u.div = 70
u.colorbar = "gist_ncar" #"CMRmap"
u.fmt = "%.3f"
u.xp = 10
u.yp = 8
u.filename = "maxustm"
u.includedate = False
u.out = "png"
u.getplot()


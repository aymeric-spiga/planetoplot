#!/bin/env python

from ppplot import plot2d
from ppclass import pp
import numpy as np

fi = '/planeto/islmd/north_pole_mesoscale/2014/5th_nest/4th_try_several_days/wrfout_d05_2024-03-46_06:00:00_zabg'
fi = 'testisaac/wrfout_d04_2024-03-46_06:00:00_zabg'
fi = 'testisaac/wrfout_d03_2024-03-46_06:00:00_zabg'
ti = 1
ti = 10
zi = 25. # in meters
deltaT = 3.
grav = 3.72
tpottest = 195. #test potential temperature for finding height of katabatic layer

### LOAD
tpot = pp(file=fi,var='tpot',t=ti,z=zi,verbose=True).getf()
Um = pp(file=fi,var='Um',t=ti,z=zi).getf()
Vm = pp(file=fi,var='Vm',t=ti,z=zi).getf()
W = pp(file=fi,var='W',t=ti,z=zi).getf()
PT,xx,yy,zz,tt = pp(file=fi,var='tpot',t=ti).getfd()
HGT = pp(file=fi,var='HGT',t=ti).getf()
potheight = pp(file=fi,var='HGT',t=ti).getf() # just to get an array the right size

### CALCULATE
vel = Um**2 + Vm**2
vel = np.sqrt(vel)
print "searching for height of katabatic layer"
for i in range(xx.shape[1]):
 for j in range(yy.shape[0]):
  w = np.where( PT[:,j,i] > tpottest)
  potheight[j,i] = zz[w][0]
print "... done."
# from pettre and andre
denom = (deltaT/tpot)*grav*potheight
denom = np.sqrt(denom)
froude = vel / denom

### PLOT
p2d = plot2d()
p2d.f = HGT
p2d.c = HGT
p2d.x = xx
p2d.y = yy
p2d.title = "Topography"
p2d.colorbar = "gist_rainbow"
p2d.makeshow()

p2d = plot2d()
p2d.f = froude
p2d.c = HGT
p2d.x = xx
p2d.y = yy
p2d.title = "Froude number"
p2d.vmin = 0.
p2d.vmax = 4.
p2d.div = 40
p2d.colorbar = "gist_stern"
p2d.fmt = "%.1f"
p2d.makeshow()

p2d = plot2d()
p2d.f = potheight
p2d.c = HGT
p2d.x = xx
p2d.y = yy
p2d.title = "Height of the katabatic layer"
p2d.vmin = 0.
p2d.vmax = 800.
p2d.div = 40
p2d.colorbar = "spectral"
p2d.fmt = "%.0f"
p2d.makeshow()

p2d = plot2d()
p2d.f = vel
p2d.c = HGT
p2d.x = xx
p2d.y = yy
p2d.title = "Horizontal wind velocity"
p2d.vmin = 0.
p2d.vmax = 20.
p2d.div = 40 
p2d.colorbar = "jet"
p2d.fmt = "%.0f"
p2d.makeshow()

p2d = plot2d()
p2d.f = W
p2d.c = HGT
p2d.x = xx
p2d.y = yy
p2d.title = "Vertical wind velocity"
p2d.vmin = -1.
p2d.vmax = 1.
p2d.div = 20
p2d.colorbar = "RdBu_r"
p2d.fmt = "%.1f"
p2d.makeshow()


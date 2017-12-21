#! /usr/bin/env python

import ppplot
import numpy as np

lons = np.linspace(-180, 180, 100)
lats = np.linspace(-90, 90, 100)

lons, lats = np.meshgrid(lons, lats)
v10 = np.ones((lons.shape)) * 15
u10 = np.zeros((lons.shape))
u10 = v10 #diagonal
v10 = 0.*v10 # zonal


pl = ppplot.plot2d()
pl.proj = "npstere"
#pl.proj = "ortho"
pl.f = lons
pl.x = lons
pl.y = lats
pl.vx = u10
pl.vy = v10

pl.makeshow()

#! /usr/bin/env python
from ppclass import pp
from ppplot import plot2d
from ppcompute import divort
import numpy as np

####################################################
fff = "/home/aymeric/Big_Data/DATAPLOT/diagfired.nc"
rad = 3389500
ttt = 1
zzz = 10.
#####################################################

# get wind fields
u = pp(file=fff,var="u",t=ttt,z=zzz).getf()
v = pp(file=fff,var="v",t=ttt,z=zzz).getf()
s = pp(file=fff,var="phisinit",t=ttt).getf()
   
# set latitude and longitude arrays
# -- asuming regular grid!
ny,nx = u.shape
lon = np.linspace(-180.,180.,nx)
lat = np.linspace(90.,-90.,ny)
   
# compute vorticity and divergence
div, vorti = divort(u,v,lon,lat,rad)

##########
# figure #
##########
pl = plot2d()
# -- mapping stuff
pl.x = lon
pl.y = lat
pl.proj = "cyl"
pl.mapmode = True
# -- field: vorticity
pl.f = vorti*1e5
pl.trans = 0.5
pl.colorbar = "RdBu_r"
pl.vmin = -5.
pl.vmax = 5.
pl.fmt = "%.0f"
# -- contour: topography
pl.c = s
# -- vectors: winds
pl.vx = u
pl.vy = v
pl.wscale = 25.
# -- make the plot
pl.makesave(mode="jpg",filename="vortimars",includedate=False)

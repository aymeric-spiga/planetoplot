#! /usr/bin/env python
from ppclass import pp
from ppplot import plot2d
from ppcompute import divort
import numpy as np

####################################################
fff = "diagfi.nc"
uuu = "u"
vvv = "v"
ttt = 1000
zzz = 1e4
rad = 60000000
#####################################################

# get wind fields
u = pp(file=fff,var="u",t=ttt,z=zzz).getf()
v = pp(file=fff,var="v",t=ttt,z=zzz).getf()
   
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
pl.f = vorti*1e6
pl.colorbar = "RdBu_r"
pl.vmin = -3.
pl.vmax = 3.
pl.fmt = "%.0f"
# -- make the plot
pl.makesave(mode="jpg",filename="vortisaturn",includedate=False)
#pl.f = div
#pl.makeshow()




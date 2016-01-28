#! /usr/bin/env python
from ppclass import pp
from ppplot import plot2d,figuref
from ppcompute import divort
import numpy as np
import planets

####################################################
myplanet = planets.Saturn
#myplanet = planets.Mars
rad = myplanet.a
####################################################
fff = "diagfi.nc"
fff = "DRAG90days_DISSIP50000_lat10_913286052-951338051_red.nc"
uuu = "u"
vvv = "v"
ttt = 1000
ttt = 1
zzz = 1e4
zzz = 1e5
#####################################################

# get wind fields
u,lon,lat,foo,foo = pp(file=fff,var="u",t=ttt,z=zzz,compute="pert_x").getfd()
v = pp(file=fff,var="v",t=ttt,z=zzz,compute="pert_x").getf()
#s = pp(file=fff,var="phisinit",t=ttt).getf()
   
# work out the pole problem 
# -- spurious high values if extrem lati close to +/-90)
ny,nx = u.shape
u,v = u[1:ny-2,:],v[1:ny-2,:]
lon,lat = lon[1:ny-2,:],lat[1:ny-2,:]

#########################################
# compute vorticity and divergence
div, vorti = divort(u,v,lon,lat,rad)
#########################################

# for plot: get highest value & highest power
absmax = np.abs(vorti).max()
exponent=int(round(np.log10(absmax)))
norm = 10.**exponent

##########
# figure #
##########
pl = plot2d()
# -- mapping stuff
pl.x = lon
pl.y = lat
pl.proj = "cyl"
#pl.mapmode = True
# -- field: vorticity
pl.f = vorti/norm
pl.units = '$10^{'+str(exponent)+'}$ s$^{-1}$'
pl.colorbar = "RdBu_r"
pl.vmin = -absmax/norm
pl.vmax = absmax/norm
pl.fmt = "%.1f"
## -- vectors: winds
#pl.vx,pl.vy = u,v
#pl.svx,pl.svy = 2,2
#pl.wscale = 50.
# -- make the plot
pl.makesave(mode="jpg",filename="vorticity",includedate=False)

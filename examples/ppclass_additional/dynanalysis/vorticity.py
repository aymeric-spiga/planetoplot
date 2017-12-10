#! /usr/bin/env python
from ppclass import pp
from ppplot import plot2d
from ppcompute import divort
import numpy as np
import planets

####################################################
myplanet = planets.Saturn
rad = myplanet.a
####################################################
#fff = "DRAG90days_DISSIP10000_year8_1842401736_512_z5"
fff = "DRAG90days_DISSIP10000_year9_z5_256_every20d_278921160-316973159"
uuu,vvv = "u","v"
ttt = 1
zzz = 1.5e5
#####################################################

# get wind fields
u,lon,lat,foo,foo = pp(file=fff+".nc",var="u",t=ttt,z=zzz,compute="pert_x").getfd()
v = pp(file=fff+".nc",var="v",t=ttt,z=zzz,compute="pert_x").getf()
#temp = pp(file=fff+".nc",var="temp",t=ttt,z=zzz,compute="pert_x").getf()
   
# work out the pole problem 
# -- spurious high values if extrem lati close to +/-90
ny,nx = u.shape
u,v = u[1:ny-2,:],v[1:ny-2,:]
lonsave,latsave = lon,lat
lon,lat = lon[1:ny-2,:],lat[1:ny-2,:]

#########################################
# compute vorticity and divergence
div, vorti = divort(u,v,lon,lat,rad)
#########################################

vorti = vorti / myplanet.fcoriolis(lat=90.)

##########
# figure #
##########
pl = plot2d()
# -- mapping stuff
pl.x = lon
pl.y = lat
pl.proj = "cyl"
## -- field: vorticity
pl.f = vorti
pl.colorbar = "seismic"
pl.vmin = -0.16
pl.vmax = +0.16
pl.fmt = "%.2f"
pl.units = r'$\times 2 \Omega$'
# -- vectors: winds
pl.vx,pl.vy = u,v
pl.svx,pl.svy = 5,5  #10,10
pl.wscale = 40.
## -- normalize
#pl.normalize()
# -- make the plot
pl.makesave(mode="png",filename=fff+"_vorticity")

#pl.f = temp
#pl.x = lonsave
#pl.y = latsave
#pl.vmin = -5.
#pl.vmax = +5.
#pl.makesave(mode="png",filename=fff+"_temp")

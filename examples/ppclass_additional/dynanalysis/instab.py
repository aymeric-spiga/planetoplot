#! /usr/bin/env python

from ppclass import pp
import ppcompute
import ppplot
import numpy as np
import planets

########################
myp = planets.Saturn
zefile="../DRAG90days_DISSIP10000_year1-10_512_every200d_zonmean_stride4lat_precast.nc"
ttt = 0.25
ttt = 7.4 
ttt = 7.45
ttt = "7.40,7.43" # mean over 5 values ~ 1000 days
ttt = "7.43,7.45"
#ttt = "7.40,7.45"
#ttt = "7.35,7.45"
########################

## fields
u,lon,lat,z,time=pp(file=zefile,var="u",x=999,t=ttt).getfd()
t=pp(file=zefile,var="temp",x=999,t=ttt).getf()

## coordinates
p=pp(file=zefile,var="p",x=999,y=999).getf()
p2d = np.transpose(np.tile(p, (u.shape[1], 1)))
lat2d = np.tile(lat, (u.shape[0], 1))
acosphi2d = myp.acosphi(lat=lat2d)
cosphi2d = acosphi2d / myp.a
latrad,lat2drad = lat*np.pi/180.,lat2d*np.pi/180.
beta = myp.beta(lat=lat2d)
f = myp.fcoriolis(lat=lat2d)

## derivatives
du_dy,dummy = ppcompute.deriv2d(u*cosphi2d,latrad,p) / acosphi2d
d2u_dy2,dummy = ppcompute.deriv2d(du_dy*cosphi2d,latrad,p ) / acosphi2d

## computations 

## BAROTROPIC INSTABILITY
## Rayleigh-Kuo criterion dqs_dy = beta - d2u_dy2 vanishes
dqs_dy = beta - d2u_dy2
dqs_dy_dim = (1. - d2u_dy2/beta)*100.

## BAROCLINIC INSTABILITY
## -- from Holton 2004 sections 8.4.1 and 8.4.2 (Rayleigh theorem)
zeh = myp.H()
pseudoz = zeh*np.log(p[0]/p) # log-pressure coordinates
dummy,dtdz = ppcompute.deriv2d(t,latrad,pseudoz)
N2 = (myp.R/zeh)*( dtdz + ((myp.R/myp.cp)*t/zeh) ) # equation (8.45) [rather than myp.N2()]
epsilon = f*f / N2 # equation (8.46)
rho = p2d / (myp.R*t)
dummy,dudz = ppcompute.deriv2d(u,latrad,pseudoz)
interm = epsilon*rho*dudz
dummy,interm = ppcompute.deriv2d(interm,latrad,pseudoz)
baroclinic = -interm/rho # equation (8.49)
baroclinic_dim = 100.*baroclinic/beta

## plot
ppplot.figuref(x=8,y=6)
pl = ppplot.plot2d()
pl.f,pl.c = dqs_dy_dim,u
pl.x,pl.y = lat2d,p2d
pl.xlabel = "Latitude"
pl.ylabel = "Pressure (Pa)"
pl.logy,pl.invert = True,True
pl.xmin,pl.xmax = -65.,-45.
pl.ymin,pl.ymax = 5e3,2e5
pl.colorbar = "hot"
pl.colorbar = "nipy_spectral"
pl.vmin,pl.vmax,pl.div = -12.,0.,24
#pl.colorbar = "seismic"
#pl.vmin,pl.vmax,pl.div = -10.,10.,40
pl.fmt = "%0.f"
pl.units = "%"
pl.make()
#ppplot.show()
ppplot.save(mode="png",filename="dimensionless-rayleigh-kuo_section")
#
ppplot.figuref(x=8,y=6)
pl.f = dqs_dy_dim+baroclinic_dim
pl.make()
ppplot.save(mode="png",filename="dimensionless-baroclinic-rayleigh_section")
#
ppplot.figuref(x=8,y=6)
pl = ppplot.plot2d()
pl.f,pl.c = u,t
pl.x,pl.y = lat2d,p2d
pl.xlabel = "Latitude"
pl.ylabel = "Pressure (Pa)"
pl.logy,pl.invert = True,True
pl.xmin,pl.xmax = -65.,-45.
pl.ymin,pl.ymax = 5e3,2e5
pl.fmt = "%0.f"
pl.units = r'm s$^{-1}$'
pl.colorbar = "seismic"
pl.vmin,pl.vmax,pl.div = -80.,80.,40
pl.make()
ppplot.save(mode="png",filename="jet_section")

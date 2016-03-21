#! /usr/bin/env python

from ppclass import pp
import ppcompute
import ppplot
import numpy as np
import planets

########################
myp = planets.Saturn
zefile="../DRAG90days_DISSIP10000_year1-10_512_every200d_zonmean_stride4lat_precast.nc"
ttt = "7.40,7.43" # mean over 5 values ~ 1000 days
                  # not more because jet is being displaced!
########################

## fields
u,lon,lat,z,time=pp(file=zefile,var="u",x=999,t=ttt).getfd()
v=pp(file=zefile,var="v",x=999,t=ttt).getf()
t=pp(file=zefile,var="temp",x=999,t=ttt).getf()
vpup=pp(file=zefile,var="vpup",x=999,t=ttt).getf()
vptp=pp(file=zefile,var="vptp",x=999,t=ttt).getf()

## coordinates
p=pp(file=zefile,var="p",x=999,y=999).getf()
p2d = np.transpose(np.tile(p, (u.shape[1], 1)))
lat2d = np.tile(lat, (u.shape[0], 1))
acosphi2d = myp.acosphi(lat=lat2d)
cosphi2d = acosphi2d / myp.a
latrad,lat2drad = lat*np.pi/180.,lat2d*np.pi/180.

## derivatives
du_dphi,du_dp = ppcompute.deriv2d(u,latrad,p)
dt_dphi,dt_dp = ppcompute.deriv2d(t,latrad,p)

## computations [cf. Andrews et al. JAS 83]
## (2.2) function psi
rcp = myp.R / myp.cp
psi = - vptp / ( (rcp*t/p2d) - (dt_dp) ) 
## (2.1) EP flux
Fphi = acosphi2d * ( - vpup + psi*du_dp ) 
## (2.3) divergence of EP flux
divFphi,dummy = ppcompute.deriv2d(Fphi*cosphi2d,latrad,p) / acosphi2d
du_dt_EPF = divFphi / acosphi2d
## (2.6) residual mean meridional circulation
dummy,dpsi_dp = ppcompute.deriv2d(psi,latrad,p)
vstar = v - dpsi_dp

## plot
ppplot.figuref(x=8,y=6)
pl = ppplot.plot2d()
pl.f,pl.c = du_dt_EPF,u
pl.x,pl.y = lat2d,p2d
pl.xlabel = "Latitude"
pl.ylabel = "Pressure (Pa)"
pl.logy,pl.invert = True,True
pl.vmin,pl.vmax,pl.div = -15,15,20
pl.xmin,pl.xmax = -65.,-45.
pl.ymin,pl.ymax = 5e3,2e5
pl.colorbar = "RdBu_r"
pl.fmt = "%0.f"
pl.normalize()
pl.make()
ppplot.save(mode="png",filename="epflux_section")
exit()
#
ppplot.figuref(x=8,y=6)
pl = ppplot.plot2d()
pl.f = vstar 
pl.c = u
pl.x = lat2d
pl.y = p2d
pl.logy = True
pl.invert = True
pl.vmin = -2
pl.vmax = 2
pl.normalize()
pl.makeshow()

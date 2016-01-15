#! /usr/bin/env python

from ppclass import pp
import ppplot
import ppcompute
import planets
import numpy as np

###############
planet=planets.Saturn
fileAP="DRAG90days_DISSIP50000_lat10_913286052-951338051_red.nc"
verb=True
factor=1e23
###############

print "get 4D fields"
u4D,longit,latit,pniv,time=pp(file=fileAP,var="u",verbose=verb).getfd()
v4D=pp(file=fileAP,var="v",verbose=verb).getf()

#print "get 4D fields (zonal anomaly)"
#primu4D=pp(file=fileAP,var="u",verbose=verb,compute="pert_x").getf()
#primv4D=pp(file=fileAP,var="v",verbose=verb,compute="pert_x").getf()

print "get 4D fields (temporal anomaly)"
primu4D=pp(file=fileAP,var="u",verbose=verb,compute="pert_t").getf()
primv4D=pp(file=fileAP,var="v",verbose=verb,compute="pert_t").getf()

print "get 2D fields (temporal+zonal mean)"
meanu2D=pp(file=fileAP,var="u",verbose=verb,t="0,1e15",x="-180,180").getf()
meanv2D=pp(file=fileAP,var="v",verbose=verb,t="0,1e15",x="-180,180").getf()
#meanmass2D=pp(file=fileAPM,var="dmass",verbose=verb,t="0,1e15",x="-180,180").getf()

##########################################
####
lat=latit[:,0]
acoslat2D=np.empty_like(meanv2D)
for i in range(len(pniv)):
  acoslat2D[i,:]=planet.acosphi(lat)
#### S dp = - rho g dz S = - g dm --> dm = - S dp / g
meanmass2D = np.empty_like(meanv2D)
nl = len(pniv)
dp = pniv[0:nl-2]-pniv[1:nl-1]
dlat,dlon = latit[1,0]-latit[0,0],longit[0,1]-longit[0,0]
surf = dlat*planet.deglength() * dlon*planet.deglength(lat)
for i in range(len(lat)):
  meanmass2D[0:nl-2,i]=surf[i]*(pniv[0:nl-2]-pniv[1:nl-1])/planet.g
  meanmass2D[nl-1,i]=surf[i]*pniv[nl-1]/planet.g
massmetric = acoslat2D*meanmass2D
####
##########################################

print "compute transport"
# [(vq)bar]   : transport total 
# [qbar][vbar]: transport mmc
# [(q'v')bar] : transport trs = [(vq)bar] - [qbar vbar]
tot_trans = ppcompute.mean(ppcompute.mean(v4D*u4D,axis=3),axis=0)*massmetric
mmc_trans = meanu2D*meanv2D*massmetric
eddy_trans = ppcompute.mean(ppcompute.mean(primv4D*primu4D,axis=3),axis=0)*massmetric

print "make plot"
fig = ppplot.figuref(x=20,y=6)
subv,subh = ppplot.definesubplot(3, fig)

pl = ppplot.plot2d()
pl.invert = True
pl.logy = True
pl.fmt = "%.1e"
pl.vmin = -2
pl.vmax = +2
pl.colorbar = "RdBu_r"
lat=latit[:,0]
pl.ylabel="Pressure (Pa)"
pl.xlabel="Latitude"

#fig.add_subplot(subv,subh,1)
#pl.f = tot_trans / factor
#pl.x = lat
#pl.y = pniv
#pl.title = "Total horizontal transport ["+str(factor)+" kg m$^3$ s$^{-2}$]"
#pl.make()
#
#fig.add_subplot(subv,subh,2)
#pl.f = mmc_trans / factor
#pl.x = lat
#pl.y = pniv
#pl.title = "MMC horizontal transport"
#pl.make()
#
#fig.add_subplot(subv,subh,3)
pl.f = eddy_trans / factor
pl.c = meanu2D
pl.x = lat
pl.y = pniv
pl.title = "Eddy horizontal transport"
pl.xmax = 60.
pl.xmin = 0.
pl.make()

ppplot.show()

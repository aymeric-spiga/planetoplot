#!/bin/env python

from ppclass import pp
import ppplot
import ppcompute
import planets


###############
planet=planets.Venus
fileAP="/donnees/sllmd/VENUS/Ete14/HR0/histmth.200_P.nc"
verb=True
verb=False
###############

print "get 4D fields"
temp4D,longit,latit,pniv,time=pp(file=fileAP,var="temp",verbose=verb).getfd()
v4D=pp(file=fileAP,var="vitv",verbose=verb).getf()

print "get 4D fields (zonal anomaly"
anotemp4D=pp(file=fileAP,var="temp",verbose=verb,compute="pert_x").getf()
anov4D=pp(file=fileAP,var="vitv",verbose=verb,compute="pert_x").getf()

print "get 4D fields (temporal mean)"
meantemp4D=pp(file=fileAP,var="temp",verbose=verb,t="0,1e15",x="-180,180").getf()
meanv4D=pp(file=fileAP,var="vitv",verbose=verb,t="0,1e15",x="-180,180").getf()

print "compute transport"
# [(vq)bar]   : transport total 
# [qbar][vbar]: transport mmc
# [(q'v')bar] : transport trs = [(vq)bar] - [qbar vbar]
tot_trans = ppcompute.mean(ppcompute.mean(v4D*temp4D,axis=3),axis=0)
mmc_trans = meantemp4D*meanv4D
eddy_trans = ppcompute.mean(ppcompute.mean(anov4D*anotemp4D,axis=3),axis=0)

print "make plot"
fig = ppplot.figuref(x=20,y=6)
subv,subh = ppplot.definesubplot(3, fig)

pl = ppplot.plot2d()
pl.invert = True
pl.logy = True
pl.fmt = "%.1f"
pl.vmin = -2
pl.vmax = +2
pl.colorbar = "RdBu_r"
lat=latit[:,0]

fig.add_subplot(subv,subh,1)
pl.f = tot_trans / 1000.
pl.x = lat
pl.y = pniv
pl.title = "Total transport /1e3"
pl.make()

fig.add_subplot(subv,subh,2)
pl.f = mmc_trans / 1000.
pl.x = lat
pl.y = pniv
pl.title = "MMC transport  /1e3"
pl.make()

fig.add_subplot(subv,subh,3)
pl.f = eddy_trans / 1000.
pl.x = lat
pl.y = pniv
pl.title = "Eddy transport /1e3"
pl.make()

ppplot.show()

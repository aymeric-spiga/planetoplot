#! /usr/bin/env python
from ppclass import pp
from ppplot import plot2d
from ppcompute import mean
import numpy as np

############################################
folder = "/home/aymeric/Remote/"
folder = "./"
suffix = ".nc"
zonalpert = True # zonal perturbations
zonalpert = False # temporal perturbations
############################################
ffftab = ["diagfired2",\
          "saturn_128x96x64_guided_hyperdiff_20000j",\
          "saturn_128x96x64_guided_hyperdiff_40000j",\
          "saturn_128x96x64_guided_hyperdiff_60000j"]
tabt = [350.,350.,350.,350.]
############################################
ffftab = ["saturn_0.5deg_dynamico"]
tabt = [1e20]
############################################
ffftab = ["DRAG90days_DISSIP50000_lat10_913286052-951338051_red"]
tabt = [300]

n = 0

for ff in ffftab:

   tt = tabt[n]
   n=n+1
   fi = folder + ff + suffix

   ###############################################################
   if zonalpert:
     # perturbation fields
     up = pp(var="u",file=fi,t=tt,compute="pert_x",changetime="correctls_noadd").getf()
     vp = pp(var="v",file=fi,t=tt,compute="pert_x",changetime="correctls_noadd").getf()
     tp = pp(var="temp",file=fi,t=tt,compute="pert_x",changetime="correctls_noadd").getf()
   else:
     dummy,dimx,dimy,dimz,dimt = pp(var="u",file=fi).getfd()
     #itemindex = np.where(dimt==tt)
     itemindex = np.argmin( np.abs( dimt - tt ) )
     # perturbation fields wrt to time
     # -- contrary to eddyflux.py this is 4D
     # -- it is not possible to reduce time dimension here
     #    because averaging / perturbing is over time
     up = pp(var="u",file=fi,compute="pert_t").getf()
     vp = pp(var="v",file=fi,compute="pert_t").getf()
     tp = pp(var="temp",file=fi,compute="pert_t").getf()
     # OK. now: once perturbation is taken over time axis
     # reduce dimension over this axis by choosing index it
     # reminder: axis are reversed: t / z / y / x
     up = up[itemindex,:,:,:]
     vp = vp[itemindex,:,:,:]
     tp = tp[itemindex,:,:,:]
   ###############################################################


   # coordinates
   try:
     press = pp(var="p",file=fi,t=tt,x=0,y=0,changetime="correctls_noadd").getf()
   except:
     press = pp(var="presnivs",file=fi,x=0,y=0,changetime="correctls_noadd").getf()
   lat = np.linspace(-90.,90.,up.shape[1])

   print press

   # compute <u'v'> (zonal mean)
   vptpm = mean(vp*tp,axis=2)
   upvpm = mean(up*vp,axis=2)

   # plot
   p = plot2d()
   p.f = upvpm
   p.c = vptpm
   p.x = lat
   p.y = press/100.
   p.title = r'$\overline{u^\prime v^\prime}$ (shaded)    $\overline{v^\prime t^\prime}$ (contours)'
   p.xlabel = "Latitude"
   p.ylabel = "Pressure (mb)"
   p.units = r'$m^{2}s^{-2}$'
   p.invert = True
   p.logy = True
   p.colorbar = "RdBu_r"
   p.fmt = "%.1e"
   p.makesave(mode="png",filename="eddyflux_"+ff)

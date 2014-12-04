#! /usr/bin/env python
from ppclass import pp
from ppplot import plot2d
from ppcompute import mean
import numpy as np

############################################
folder = "/home/aymeric/Remote/"
suffix = ".nc"
############################################
ffftab = ["diagfired2",\
          "saturn_128x96x64_guided_hyperdiff_20000j",\
          "saturn_128x96x64_guided_hyperdiff_40000j",\
          "saturn_128x96x64_guided_hyperdiff_60000j"]
tabt = [350.,350.,350.,350.]
############################################

ffftab = ["saturn_0.5deg_dynamico"]
tabt = [1e20]


n = 0

for ff in ffftab:

   tt = tabt[n]
   n=n+1
   fi = folder + ff + suffix

   # perturbation fields
   up = pp(var="u",file=fi,t=tt,compute="pert_x",changetime="correctls_noadd").getf()
   vp = pp(var="v",file=fi,t=tt,compute="pert_x",changetime="correctls_noadd").getf()
   tp = pp(var="temp",file=fi,t=tt,compute="pert_x",changetime="correctls_noadd").getf()

   # coordinates
   press = pp(var="p",file=fi,t=tt,x=0,y=0,changetime="correctls_noadd").getf()
   lat = np.linspace(-90.,90.,up.shape[1])

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

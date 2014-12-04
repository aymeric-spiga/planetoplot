#! /usr/bin/env python
from ppclass import pp
from ppplot import plot2d
from ppcompute import mean
import numpy as np

#############################
ff = "diagfi.nc"
it = 44
#############################

# perturbation fields wrt to time
# -- contrary to eddyflux.py this is 4D
# -- it is not possible to reduce time dimension here
#    because averaging / perturbing is over time
up = pp(var="u",file=ff,t=it,compute="pert_t").getf()
vp = pp(var="v",file=ff,compute="pert_t").getf()
tp = pp(var="temp",file=ff,compute="pert_t").getf()

# OK. now: once perturbation is taken over time axis
# reduce dimension over this axis by choosing index it
# reminder: axis are reversed: t / z / y / x
up = up[it,:,:,:]
vp = vp[it,:,:,:]
tp = tp[it,:,:,:]

# coordinates
press = pp(var="p",file=ff,t=tt,x=0,y=0).getf()
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
#p.title = r'$\overline{u^\prime v^\prime}$'
p.title = "" 
p.xlabel = "Latitude"
p.ylabel = "Pressure (mb)"
p.units = r'$m^{2}s^{-2}$'
p.invert = True
p.logy = True
p.colorbar = "RdBu_r"
p.fmt = "%.1f"
p.vmin = -7
p.vmax = 7
p.div = 28
p.makesave(mode="png",filename="eddyflux2")

#! /usr/bin/env python
from ppclass import pp

#############################
xx = 180.
xx = 0. 
xm = "-175.,180."
ff = "diagfi.nc"
tt = 240.
#############################

# field
u = pp(var="u",file=ff,t=tt,x=xx,quiet=True).get()
v = pp(var="v",file=ff,t=tt,x=xx,quiet=True).get()
t = pp(var="temp",file=ff,t=tt,x=xx,quiet=True).get()

# zonal mean
um = pp(var="u",file=ff,t=tt,x=xm,quiet=True).get()
vm = pp(var="v",file=ff,t=tt,x=xm,quiet=True).get()
tm = pp(var="temp",file=ff,t=tt,x=xm,quiet=True).get()

# generic plot settings (will be passed through operations)
u.logy = True
u.nxticks = 5
u.ycoeff = 0.01
u.ylabel = "pressure (mb)"
u.out = "png"
# (what follows changes key attributes but harmless in our case, while useful for plot settings)
v << u ; t << u

# perturbation
vp = v-vm
up = u-um
tp = t-tm

# coupling terms
vptp = vp*tp
upvp = up*vp
uptp = up*tp

# specific plot settings
vp.title = r"$v^\prime$"
up.title = r"$u^\prime$"
tp.title = r"$T^\prime$"
vptp.title = r"$v^\prime T^\prime$"
vptp.units = "ms$^{-1}$K"
vptp.colorbar = "RdBu_r"
upvp.title = r"$u^\prime v^\prime$"
upvp.units = "m$^2$s$^{-2}$"
upvp.colorbar = "RdBu_r"
uptp.title = r"$u^\prime T^\prime$"
uptp.units = "ms$^{-1}$K"
uptp.colorbar = "RdBu_r"

# composite plots
tp.filename = "saturn_perturbations"
tp.plot(extraplot=2)
up.plotin = tp ; up.plot()
vp.plotin = tp ; vp.plot()
vptp.filename = "saturn_coupling"
vptp.plot(extraplot=2)
uptp.plotin = vptp ; uptp.plot()
upvp.plotin = vptp ; upvp.plot()

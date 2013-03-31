#! /usr/bin/env python
from ppclass import pp

u = pp()
u.file = "/home/aymeric/Big_Data/DATAPLOT/diagfired.nc"
u.var = "u"
u.t = "0.5,0.8"
u.z = "10,20"
u.getdefineplot(extraplot=2) # prepare 2 extraplots (do not show)
u.p[0].proj = "ortho"
u.p[0].title = "$u$ (m s$^{-1}$)"
u.makeplot()

v = pp()
v << u  # NB: initialize v object with u object's attributes
v.var = "v"
v.get()
v.plotin = u # plotin must be defined before .defineplot()
v.defineplot()
v.p[1].proj = "ortho"
v.p[1].title = "$v$ (m s$^{-1}$)"
v.makeplot() # plot within the previous one (do not show)

wind = u**2 + v**2
wind = wind**0.5
wind.plotin = v
wind.filename = "windspeed"
wind.defineplot()
wind.p[2].title = "$\sqrt{u^2+v^2}$ (m s$^{-1}$)"
wind.p[2].proj = "ortho"
wind.makeplot() # plot within the previous one (show because complete)


## in this case it is not possible to make another plot afterwards...


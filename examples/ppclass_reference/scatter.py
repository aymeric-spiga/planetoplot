#! /usr/bin/env python
from ppclass import pp

tsurf = pp()
tsurf.file = "/home/aymeric/Big_Data/DATAPLOT/diagfired.nc"
tsurf.var = "tsurf"
tsurf.x = None
tsurf.y = 10.
tsurf.t = 2.
tsurf.getdefineplot()

ps = pp()
ps << tsurf
ps.var = "ps"
ps.getdefineplot()

S = ps.func(tsurf)
S.p[0].linestyle=""
S.p[0].marker="h"
S.p[0].color="g"
S.makeplot()

icetot = pp()
icetot << tsurf
icetot.var = "icetot"
icetot.getdefineplot()

S2 = icetot.func(tsurf)
S2.p[0].linestyle=""
S2.p[0].marker="D"
S2.p[0].color="r"
S2.makeplot()

u = pp()
u << tsurf
u.var = "u"
u.z = 1.
u.get()

v = pp()
v << u
v.var = "v"
v.get()

wind = u**2 + v**2
wind = wind**0.5
S3 = wind.func(ps)
S3.p[0].linestyle=""
S3.p[0].marker="o"
S3.p[0].color="k"
S3.p[0].ylabel="wind speed $\sqrt{u^{2}+v^{2}}$ (m s$^{-1}$)"
S3.filename="scatter"
S3.makeplot()

## multidim scatter also possible
## the important thing is forcedimplot

tsurf = pp()
tsurf.file = "/home/aymeric/Big_Data/DATAPLOT/diagfired.nc"
tsurf.var = "tsurf"
tsurf.x = None
tsurf.y = None
tsurf.t = 2.
tsurf.forcedimplot = 1
tsurf.getdefineplot()

ps = pp()
ps << tsurf
ps.var = "ps"
ps.getdefineplot()

S = ps.func(tsurf)
S.p[0].linestyle=""
S.p[0].marker="h"
S.p[0].color="g"
S.makeplot()

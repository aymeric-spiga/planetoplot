#! /usr/bin/env python
from ppclass import pp

gw = pp()
gw.file = "/home/aymeric/Big_Data/GALE/wrfout_d03_2024-06-09_00:00:00_z"
gw.var = ["W","tk"]
gw.vargoal = ["main","contour"]

gw.x = None
gw.y = -5.
gw.z = None
gw.t = 5.
gw.verbose = True

#gw.getdefineplot()
#gw.p[0].colorb = "RdBu_r"
#gw.p[0].vmin = -1.
#gw.p[0].vmax = 1.
#gw.p[0].div = 30
#gw.p[0].xlabel = "East longitude (deg)"
#gw.p[0].ycoeff = 1./1000.
#gw.p[0].ylabel = "Altitude above MOLA reference (km)"
#gw.makeplot()

gw.vargoal = ["main","main"]
gw.out = "png"
gw.getdefineplot()
for plot in gw.p:
    plot.xlabel = "East longitude (deg)"
    plot.ycoeff = 1./1000.
    plot.ylabel = "Altitude above MOLA reference (km)"
    plot.div = 30
    plot.nxticks = 5
gw.p[0].colorbar = "RdBu_r"
gw.p[0].vmin = -1.
gw.p[0].vmax = 1.
gw.makeplot()


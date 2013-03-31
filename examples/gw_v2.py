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
gw.filename = "gw"

# we define global settings for plot
# (must be done before defineplot)
gw.colorb = "RdBu_r"
gw.vmin = -1. ; gw.vmax = 1.
gw.div = 30
gw.xlabel = "East longitude (deg)"
gw.ycoeff = 1./1000.
gw.ylabel = "Altitude above MOLA reference (km)"

# we define, get field and define, get plot
gw.getplot()

# we keep the same object for another plot of a different kind
gw.vargoal = ["main","main"]
gw.out = "png"
gw.getdefineplot()

# NB: global plot settings are defined in gw !!
# ... so we redefine plot settings for the plot on the right
gw.p[1].define_from_var() # we define colors, etc... from variable
gw.p[1].vmin = None # we reinitialize to default for min bounds
gw.p[1].vmax = None # we reinitialize to default for max bounds
gw.makeplot()

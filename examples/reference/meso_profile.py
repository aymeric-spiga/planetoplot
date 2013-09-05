#! /usr/bin/env python
from ppclass import pp

# define object, file, var
m = pp()
m.file = "/home/aymeric/Big_Data/GALE/wrfout_d03_2024-06-09_00:00:00"
m.var = "W"

# define dimensions
m.x = "136.,139." # computing over x interval
m.y = -5. # setting a fixed y value
m.z = None # leaving z as a free dimension
m.t = [6.,9.,12.,15.,18.,21.,24.] # setting 4 fixed t values

# define settings
m.superpose = True # superpose 1D plots
m.verbose = True # making the programe verbose
#m.out = "pdf" # output format
m.colorb = "spectral" # color cycle according to a color map

# get data and make plot with default settings
m.getplot()

# get potential temperature at same point. 
# don't plot it. do an operation on it.
tpot = pp()
tpot << m
tpot.var = "T"
tpot.get()
tpot = tpot + 220.

# get geopotential at same point.
# don't plot it. do an operation on it (to get height).
geop = pp()
geop << m
geop.var = "PHTOT"
geop.get()
z = geop/3.72/1000.

# define potential temperature as a function of height
S = tpot.func(z)

# change a few plot settings
for curve in S.p: 
    curve.lstyle = "--"
    curve.marker = ""
S.p[0].swaplab = False
S.p[0].ylabel="Geopotential height (km)"
S.p[0].xlabel="Potential temperature (K)"
S.filename = "meso_profile"
S.colorb = None # come back to default color cycle

# make the plot
S.makeplot()

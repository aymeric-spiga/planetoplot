#! /usr/bin/env python
from ppclass import pp

prof = pp()
prof.verbose = True
#prof.file = "/home/aymeric/Big_Data/case_night.nc"
prof.file = "/home/aymeric/Big_Data/GALE/wrfout_d03_2024-06-09_00:00:00_z"
prof.var = "tk"
prof.t = [10.,20.]
prof.x = -5. #30
prof.y = 138. #30
prof.filename = "vertpress"

prof.get()
prof.superpose = True
prof.defineplot()
prof.p[0].logy = True
prof.makeplot()

## a small workaround if one wants to change axis names
## (it is difficult because of swapping axis)
prof.p[0].swaplab = False
prof.p[0].xlabel = "Temperature (K)"
prof.p[0].ylabel = "Pressure (Pa)"
prof.makeplot()





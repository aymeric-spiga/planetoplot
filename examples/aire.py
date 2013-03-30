#! /usr/bin/env python
from ppclass import pp

tpond = pp()
tpond.file = "/home/aymeric/Big_Data/DATAPLOT/diagfired.nc"
tpond.var = "temp"
tpond.t = "0.5"
tpond.x = "-180,180"
tpond.y = "-90,90"
tpond.verbose = True
tpond.compute = "meanarea"
tpond.superpose = True
tpond.getplot(extraplot=1)

tnormal = pp()
tnormal << tpond
tnormal.compute = "mean"
tnormal.plotin = tpond
tnormal.getplot()


#! /usr/bin/env python
from ppclass import pp
from ppplot import plot2d
import numpy as np

## load streamfunction from another file
## ** USE kind3d="tzy" --> because 3D file is different
pfield = pp(file="histmth.1401_PSI.nc",x=0,t=10,var="psi",verbose=True,kind3d="tzy").getf()

## get field and define a plot
t = pp()
t.file = "histmth.1401.nc"
t.x = "-180,180"
t.t = 10
t.var = "temp"
t.logy = True
t.get()
t.defineplot()

## add the streamfunction as contour
t.p[0].c = pfield

## make plot
t.makeplot()

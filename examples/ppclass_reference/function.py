#! /usr/bin/env python
from ppclass import pp

var1 = pp(\
file="/home/aymeric/Big_Data/DATAPLOT/diagfired.nc",\
var="ps",\
x=None,\
y=10.,\
t=2.)
var1.get()

var2 = pp(\
file="/home/aymeric/Big_Data/DATAPLOT/diagfired.nc",\
var="phisinit",\
x=None,\
y=10.,\
t=2.)
var2.get()
var2 = var2 / 3.72

S = var2.func(var1)

S.p[0].marker = 'o'
S.p[0].linestyle = ''
S.p[0].ylabel = "Surface geopotential height (km)"
S.p[0].ycoeff = 1./1000.
S.p[0].fmt = '%.0f'
S.filename = "function"
S.makeplot()


#! /usr/bin/env python
from ppclass import pp

## AVERAGED PROFILE
temp = pp()
temp.file = "/home/aymeric/Big_Data/DATAPLOT/diagfired.nc"
temp.var = "temp"
temp.x = "-180,175"
temp.y = "-90,90"
temp.t = "0,1"
temp.get()
temp.defineplot()
temp.title = "This is an averaged temperature profile"
temp.makeplot()

## ZONAL MEAN
u = pp()
u.file = "/home/aymeric/Big_Data/DATAPLOT/diagfired.nc"
u.var = "u"
u.x = "-180,175"
u.y = None
u.t = "0.5"
u.filename = "zonalmean"
u.get()
u.defineplot()
u.div = 30.
u.colorbar = "RdBu_r"
u.title = "This is a zonal mean"
u.makeplot()

## ZONAL MINIMUM
u.compute = "min"
u.filename = "zonalmin"
u.get()
u.defineplot()
u.div = 30.
u.colorbar = "cool"
u.title = "This is minimum over zonal axis"
u.makeplot()

## ZONAL MAXIMUM
u.compute = "max"
u.filename = "myplot"
u.get()
u.defineplot()
u.div = 30.
u.colorbar = "hot"
u.title = "This is maximum over zonal axis"
u.makeplot()

## MAP OF MAXIMUM OVER TIME
u.compute = "max"
u.x = None
u.y = None
u.t = "0,1"
u.z = 20000.
u.get()
u.defineplot()
u.div = 30.
u.title = "This is maximum over time"

u.out = "png"

u.makeplot()

#! /usr/bin/env python
from ppclass import pp

## A SIMPLE EXAMPLE
simple = pp()
simple.file = "/home/aymeric/Big_Data/DATAPLOT/diagfired.nc"
simple.var = "tsurf"
simple.t = "0.5"
simple.getplot()

## GET A TIME AVERAGE WITHOUT A PLOT
mean = pp()
mean << simple
mean.t = "0,1"
mean.get()

## COMPUTE AND PLOT RELATIVE ANOMALY in %
anomaly = ((simple-mean)/mean)*100.
anomaly.filename = "anomaly"
anomaly.defineplot()
anomaly.p[0].title = "surface temperature anomaly in %"
anomaly.makeplot()

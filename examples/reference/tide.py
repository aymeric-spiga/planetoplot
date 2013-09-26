#! /usr/bin/env python
from ppclass import pp

fi = "/home/aymeric/Big_Data/DATAPLOT/diagfired.nc"
v = ["ps","phisinit"]
vg = ["main","contour"]

meanps = pp(file=fi,var=v,vargoal=vg,t="0,1").get()
waveref = pp(file=fi,var=v,vargoal=vg,t=0.5).get() - meanps
waveref.smooth(5)
waveref.title = "$p_s$ diurnal anomaly"
waveref.proj = "moll"
waveref.vmin = -10
waveref.vmax = 10
waveref.div = 10
waveref.filename = "tide"
waveref.plot(extraplot=3)

eqmeanps = pp(file=fi,var="ps",t="0,1",y=0.).get() 
addplot = pp(file=fi,var="ps",t=0.5,y=0.).get() - eqmeanps
addplot.plotin = waveref
addplot.title = "Tide signal at the equator"
addplot.ylabel = "$p_s$ diurnal anomaly (Pa)"
addplot.linestyle = ""
addplot.marker = "."
addplot.plot()

for choice_t in [0.6,0.7]:
    wave = pp(file=fi,var=v,vargoal=vg,t=choice_t).get() - meanps
    wave << waveref
    wave.smooth(5)
    wave.plotin = waveref
    wave.plot()








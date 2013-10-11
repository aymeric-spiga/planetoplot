#! /usr/bin/env python
import ppplot
import numpy as np

# load data
# http://docs.scipy.org/doc/numpy/reference/generated/numpy.loadtxt.html
p,z,t,dew,hum,mix,wdir,wknot,thta,thte,thtv \
    = np.loadtxt("sounding.txt",skiprows=8,unpack=True)

# make a simple plot
sdg = ppplot.plot1d()
sdg.f = z
sdg.x = t
sdg.makeshow()

# Tweaking a little bit the plot
sdg.linestyle = '-'
sdg.marker = '.'
sdg.color = 'r'
sdg.ycoeff = 1.e-3
sdg.title = "A random terrestrial sounding"
sdg.legend = "Fort Smith 00Z 26 Sep 2013"
sdg.xlabel = "Temperature ($^{\circ}$C)"
sdg.ylabel = "Altitude (km)"
sdg.makeshow()

sdg.make()
ppplot.save(mode="png",filename="plot",res=50,includedate=False)



#! /usr/bin/env python

from ppclass import pp
import planets
import ppplot

fi = "DRAG90days_DISSIP10000_year7_912791376_512_z5"
u,lon,lat,p,t = pp(file=fi+".nc",var="u",t=0,z=0,x="-180,180").getfd()

sindex = planets.Saturn.superrot(u=u,lat=lat)

pl = ppplot.plot1d()
pl.f = sindex*1000.
pl.x = lat
pl.ylabel = r'Local superrotation index $s = \mathcal{M}_u / \Omega a^2$ ($\times 10^{-3}$)'
pl.xmin = -10.
pl.xmax = 10.
pl.fmt = '%.0f'
pl.ymin = -30.
pl.ymax = 10.
pl.marker = ''
pl.xlabel = "Latitude (deg)"
pl.make()
ppplot.save(mode="png",filename=fi+"_superrot")





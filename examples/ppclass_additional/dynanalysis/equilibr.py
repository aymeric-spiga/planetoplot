#! /usr/bin/env python

#pp.py precast_170.nc -t 0,1e10 -v acceddh -v accrmch -z 8e4 -x 999 -S -K None


from ppclass import pp
import ppcompute
import ppplot
zz = 8e4
zz = 3e3
ww = 30 #20
fifi = "precast_170.nc"
#fifi = "precast_200.nc"
fifi = "precast_171.nc"

edd,x,y,z,t = pp(file=fifi,t = "0,1e10",x = 999,z = zz,var = "acceddh").getfd()
edds = ppcompute.smooth1d(edd,window=ww)
ys = ppcompute.smooth1d(y,window=ww)

rmc = pp(file=fifi,t = "0,1e10",x = 999,z = zz,var = "accrmch").getf()
rmcs = ppcompute.smooth1d(rmc,window=ww)

#acc = pp(file=fifi,t = "0,1e10",x = 999,z = zz,var = "dudt").getf()
#accs = ppcompute.smooth1d(acc,window=ww)
#
ustart = pp(file=fifi,var="u"    ,t=0   ,x = 999, z=zz).getf()
uend   = pp(file=fifi,var="u"    ,t=1e10,x = 999, z=zz).getf()
dudt = (uend - ustart) / (1000.*38052.)
dudts = ppcompute.smooth1d(dudt,window=ww)


pl = ppplot.plot1d()
pl.f = edds
pl.x = ys
pl.marker = None
#pl.xmin = -80.
pl.xmin = 0.
pl.xmax = 80.
pl.fmt = "%.1e"
pl.legend = "eddies"
pl.make()
pl.f = rmcs
pl.legend = "residual circulation"
#pl.make()
#pl.f = dudts
#pl.make()
#pl.f = accs
pl.makeshow()

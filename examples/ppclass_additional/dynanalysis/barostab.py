#! /usr/bin/env python
from ppclass import pp
import ppplot
import ppcompute
import numpy as np
import planets

####################################################
myplanet = planets.Saturn
####################################################
fff = "/home/aymeric/Remote/diagfi.nc"
fff = "/home/aymeric/Remote/LATEST/saturn_128x96x64_guided_hyperdiff.nc"
fff = "/home/aymeric/Remote/saturn_128x96x64_guided_hyperdiff_20000j_last.nc"
fff = "/home/aymeric/Remote/saturn_128x96x64_guided_beg.nc"
fff = "/home/aymeric/Remote/red_60000j.nc"
#fff = "/home/aymeric/Remote/red_20000j.nc"
#fff = "/home/aymeric/Remote/yorgl.nc"
#fff = "/home/aymeric/Remote/saturn_0.5deg_dynamico.nc"
ttt = 1000
#####################################################

# useful terms
a = myplanet.a
conv = a*np.pi/180.
omega = 2.*np.pi/myplanet.day

# get coordinates and zonal mean temperature + zonal mean zonal wind
temp,lon,lat,p,t = pp(file=fff,var="temp",t=ttt,x="-180,180").getfd()
u = pp(file=fff,var="u",t=ttt,x="-180,180").getf()
press = pp(file=fff,var="p",t=ttt,x="-180,180").getf()

# convert deg to rad
phi = lat*np.pi/180.
sinphi = np.sin(phi)
cosphi = np.cos(phi)

# convert pressure to pseudo-altitude
p0 = 1.e5
z = myplanet.H()*np.log(p0/p)

# vertical derivatives
dTdphi,dTdz = ppcompute.deriv2d(temp,phi,z)
dudphi,dudz = ppcompute.deriv2d(u,phi,z)

# convert horizontal derivatives
dTdy = dTdphi/a
dudy = dudphi/a

# compute Brunt-Vaisala frequency
bv0 = myplanet.N2() 
bv = myplanet.N2(T0=temp,dTdz=dTdz)

# double horizontal derivative for u
dudydy,dummy = ppcompute.deriv2d(dudy,phi,z,fac=a)

# beta term
betaterm = 2.*omega*a*cosphi
# stability term 
dummy,sterm = ppcompute.deriv2d(dudz*press/bv,phi,z)
sterm = -sterm*4.*omega*omega*a*a*(sinphi**2)/press
# barotropic term
cosphi[np.where(cosphi < 1.e-15)] = 1. # to avoid problem with divide by zero at poles
term = cosphi*u
dterm,dummy = ppcompute.deriv2d(term,phi,z)
d2term,dummy = ppcompute.deriv2d(dterm/cosphi,phi,z)
barterm = -d2term

# meridional gradient of PV
tot = betaterm + barterm + sterm
#tot = barterm + betaterm
#tot = (2.*omega/a) - dudydy

## PLOT
pl=ppplot.plot2d()
pl.x = lat
pl.y = p0*np.exp(-z/myplanet.H())
pl.colorbar = "RdBu_r"
pl.logy = True
pl.invert = True

#pl.f = temp ; pl.c = tot


pl.f = tot
mm = np.min([np.abs(pl.f.min()),np.abs(pl.f.max())])
pl.vmin = -mm/2. ; pl.vmax = mm/2.
#pl.vmin = -mm/10. ; pl.vmax = mm/10.

pl.xlabel = "Latitude"
pl.ylabel = "Pressure (Pa)"

pl.makeshow()


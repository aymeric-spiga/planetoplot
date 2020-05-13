#! /usr/bin/env python
from ppclass import pp
import ppplot
import numpy as np
from ppcompute import smooth2diter

## options
from optparse import OptionParser ### TBR by argparse
parser = OptionParser()
parser.usage = "les.py netCDF_file(s)" #les.py becomes a commandline applicable to a netcdf file
(opt,files) = parser.parse_args()

## constants
t0 = 220.
p0 = 610.
r_cp = 1.0/4.4
grav = 3.72
R = 192.

## dimensions
foo1,foo2,foo3,z,t = pp(file=files[0],var="T",x=0,y=0).getfd()
nz,nt = z.size,t.size
nf = len(files)
print(nf)
ntt = nt*nf

## arrays for time series
tprimemean,tmean = np.zeros([ntt,nz]),np.zeros([ntt,nz])
wprimemean,wmean = np.zeros([ntt,nz]),np.zeros([ntt,nz])
wmax = np.zeros([ntt,nz])
vehfmean = np.zeros([ntt,nz])
pblh,pblh1,pblh2,pblh3 = np.zeros(ntt),np.zeros(ntt),np.zeros(ntt),np.zeros(ntt)

## loop
indt = -1

for ff in files:

 print("### PROCESSING FILE:",ff)
 print("### --> T")
 t = pp(file=ff,var="T").getf() + t0
 print("### --> W")
 w = pp(file=ff,var="W").getf()
 print("### --> PHTOT")
 geop = pp(file=ff,var="PHTOT").getf()
 geop = geop - np.mean(geop[:,0,:,:])
 
 for tt in range(nt):

  indt = indt + 1
 
  for zz in range(nz):
    #print tt, indt
    tprime = t[tt,zz,:,:]
    tmean[indt,zz] = np.mean(tprime)
    tprime = tprime - tmean[indt,zz]
 
    wprime = w[tt,zz,:,:]
    wmean[indt,zz] = np.mean(wprime)
    wmax[indt,zz] = np.max(np.abs(wprime))
    wprime = wprime - wmean[indt,zz]
 
    vehfmean[indt,zz] = np.mean(tprime*wprime)
    tprimemean[indt,zz] = np.mean(tprime)
    wprimemean[indt,zz] = np.mean(wprime)

  ###########
  ## method 1: potential temperature profile
  diff = np.abs(tmean[indt,2:]-tmean[indt,1])
  wheremin = np.argmin(diff) + 2
  print(indt, tt)
  pblh1[indt] = np.mean(geop[tt,wheremin,:,:]) / grav #Height of the geopotential of the boundary layer
 
  ###########
  ## method 2: minimum of vertical eddy heat flux
  wheremin = np.argmin(vehfmean[indt,:]) 
  pblh2[indt] = 1.1*np.mean(geop[tt,wheremin,:,:]) / grav
  ## sometimes spurious values caused by GW
  diff = 100.*np.abs(pblh2[indt]-pblh1[indt])/pblh1[indt]
  if diff > 33.: pblh2[indt] = pblh1[indt]

  ###########
  ## method 3: convective motions
  diff = np.abs(wmax[indt,8:]-np.max(wmax[indt,:]/3.))
  wheremin = np.argmin(diff) + 8
  pblh3[indt] = np.mean(geop[tt,wheremin,:,:]) / grav
  ## sometimes spurious values caused by GW
  diff = 100.*np.abs(pblh3[indt]-pblh1[indt])/pblh1[indt]
  if diff > 33.: pblh3[indt] = pblh1[indt]

  pblh[indt] = np.mean([pblh1[indt],pblh2[indt],pblh3[indt]])

## remove small or negative values 
pblh = pblh[pblh > 100.]

## compute mean height
altitude = np.mean(np.mean(np.mean(geop,axis=3),axis=2),axis=0)/grav
altitude = altitude[0:nz]

file1=open('pbl.txt','w')
for val in pblh:
  file1.write("%8.3f\n"%(val))
file1.close()

pl = ppplot.plot1d() # plot of the boundary layer height time evolution
pl.f = pblh
pl.makeshow()

pl = ppplot.plot2d() # shade of the vertical eddy heat flux time evolution
pl.f = vehfmean
pl.y = altitude
pl.makeshow()



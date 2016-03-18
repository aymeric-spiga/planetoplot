#! /usr/bin/env python

import numpy as np
from ppclass import pp
import netCDF4 as nc

####################################################
fileAP="DRAG90days_DISSIP10000_year1-10_512_every200d_zonmean_stride4lat.nc"
####################################################
p_upper,p_lower,nlev = 4.0e2,2.5e5,40 # whole atm
#p_upper,p_lower,nlev = 3e3, 2e5, 40 # tropo
targetp1d = np.logspace(np.log10(p_lower),np.log10(p_upper),nlev)
####################################################
targetp1d = np.array([5.,10.,20.,50.,100.,200.,500.,1000.,2000.])
targetp1d = targetp1d*100.
####################################################

#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------

####################################################
def interpolate(targetp1d,sourcep3d,fieldsource3d):
  nt,nz,nlat = fieldsource3d.shape
  coordsource3d = -np.log(sourcep3d) # interpolate in logp
  coordtarget1d = -np.log(targetp1d) # interpolate in logp
  nzt = coordtarget1d.size
  fieldtarget3d = np.zeros((nt,nzt,nlat))
  for nnn in range(nlat):
   for ttt in range(nt):
    fieldtarget3d[ttt,:,nnn] = np.interp(coordtarget1d,coordsource3d[ttt,:,nnn],fieldsource3d[ttt,:,nnn],left=np.nan,right=np.nan)
  return fieldtarget3d

####################################################
def fix_time_axis(tdim,period):
  ntt = tdim.size
  tdim = tdim % period 
  nperiod = 0
  corrected_tdim = np.empty(ntt)
  for iii in range(1,ntt):
    if tdim[iii] - tdim[iii-1] < 0: 
      nperiod = nperiod + 1
    corrected_tdim[iii] = float(nperiod) + float(tdim[iii]/period)
  return corrected_tdim

####################################################
print "... getting fields from file !"
press,xdim,ydim,zdim,tdim=pp(file=fileAP,var="p",x=999).getfd()
tdim = fix_time_axis(tdim,24430.)
u=pp(file=fileAP,var="u",x=999).getf()
v=pp(file=fileAP,var="v",x=999).getf()
temp=pp(file=fileAP,var="temp",x=999).getf()
vpup=pp(file=fileAP,var="vpup",x=999).getf()
vptp=pp(file=fileAP,var="vptp",x=999).getf()
uptp=pp(file=fileAP,var="uptp",x=999).getf()
upup=pp(file=fileAP,var="upup",x=999).getf() 
vpvp=pp(file=fileAP,var="vpup",x=999).getf() 
eke = 0.5 * (vpvp + upup)

####################################################
print "... interpolating !"
ui = interpolate(targetp1d,press,u)
vi = interpolate(targetp1d,press,v)
tempi = interpolate(targetp1d,press,temp)
vpupi = interpolate(targetp1d,press,vpup)
vptpi = interpolate(targetp1d,press,vptp)
uptpi = interpolate(targetp1d,press,uptp)
ekei = interpolate(targetp1d,press,eke)

####################################################
print "... creating the target file !"

f = nc.Dataset("precast.nc",'w',format='NETCDF3_CLASSIC')

nt,nz,nlat = ui.shape
nlon = 1
xdim = [999.]
dimx = 'longitude'
dimy = 'latitude'
dimz = 'p'
dimt = 'time_counter'

## this is OK with ncview-like tools
dimz = 'pseudoalt'
scaleheight = 50.
psalt = -scaleheight*np.log(targetp1d/targetp1d[0])

nam4 = (dimt,dimz,dimy,dimx)
shp4 = (nt,nz,nlat,nlon)
fie4 = (tdim,psalt,ydim,xdim)

for iii in range(len(shp4)):
  f.createDimension(nam4[iii],shp4[iii])
  var = f.createVariable(nam4[iii], 'd', nam4[iii])
  var[:] = fie4[iii]

## this is OK with planetoplot
var = f.createVariable('p', 'd', nam4[1])
var[:] = targetp1d
var = f.createVariable('u',  'd', nam4)
var[:,:,:,:]          = ui
var = f.createVariable('v',  'd', nam4)
var[:,:,:,:]          = vi
var = f.createVariable('temp',  'd', nam4)
var[:,:,:,:]          = tempi
var = f.createVariable('vpup',  'd', nam4)
var[:,:,:,:]          = vpupi
var = f.createVariable('vptp',  'd', nam4)
var[:,:,:,:]          = vptpi
var = f.createVariable('uptp',  'd', nam4)
var[:,:,:,:]          = uptpi
var = f.createVariable('eke',  'd', nam4)
var[:,:,:,:]          = ekei

####################################################
print "... adding 2D variables !"
ISR=pp(file=fileAP,var="ISR",x=999.).getf()
OLR=pp(file=fileAP,var="OLR",x=999.).getf()
var = f.createVariable('ISR',  'd', (nam4[0],nam4[2],nam4[3]))
var[:,:,:] = ISR
var = f.createVariable('OLR',  'd', (nam4[0],nam4[2],nam4[3]))
var[:,:,:] = OLR

f.close()

#! /usr/bin/env python

import numpy as np
from ppclass import pp
import ppcompute
import netCDF4 as nc
import planets
import time

fileAP="Xhistins_1900.nc"
#fileAP="Xhistins_450.nc"
fileAP="Xhistins_225.nc"
p_upper,p_lower,nlev = 4.0e2,2.5e5,40
p_upper,p_lower,nlev = 1e3,1e5,30
p_upper,p_lower,nlev = 1e2,3.5e5,40 # whole atm
targetp1d = np.logspace(np.log10(p_lower),np.log10(p_upper),nlev)
myp = planets.Saturn
charx = "0,360"
vartemp = "temperature"
outfile = "precast.nc"
use_spline = False

#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------

####################################################
def etape(charvar,time0):
  ttt = round(time.time()-time0,2)
  print "TIME=",ttt,"... ... done: "+charvar

####################################################
def interpolate(targetp1d,sourcep3d,fieldsource3d,spline=False):
  if spline:
    from scipy import interpolate
  nz,nlat = fieldsource3d.shape
  coordsource3d = -np.log(sourcep3d) # interpolate in logp
  coordtarget1d = -np.log(targetp1d) # interpolate in logp
  nzt = coordtarget1d.size
  fieldtarget3d = np.zeros((nzt,nlat))
  for nnn in range(nlat):
    xs = coordsource3d[:,nnn]
    ys = fieldsource3d[:,nnn]
    if not spline:
      fieldtarget3d[:,nnn] = np.interp(coordtarget1d,xs,ys,left=np.nan,right=np.nan)
    else:
      tck = interpolate.splrep(xs, ys, s=0)
      fieldtarget3d[:,nnn] = interpolate.splev(coordtarget1d, tck, der=0)
  return fieldtarget3d

####################################################
def addvar(filename,dimname,varname,varfield,time0=None):
  f = nc.Dataset(filename,'a',format='NETCDF3_CLASSIC')
  var = f.createVariable(varname, 'd', dimname) 
  if   len(dimname) == 4: var[:,:,:,:] = varfield
  elif len(dimname) == 3: var[:,:,:] = varfield
  elif len(dimname) == 2: var[:,:] = varfield
  varfield = None ; var = None
  f.close()
  if time0 is not None:
    etape(varname,time0)
  return

####################################################
####################################################
####################################################
####################################################

#lstab = pp(file=fileAP,var="ls",x=0,y=0,z=0).getf()
#lstab = lstab*180./np.pi
#lstab = fix_time_axis(lstab,360.) # in years

####################################################
time0 = time.time()

psbc = pp(file=fileAP,var="ps",t="0,1e10",x=charx).getf()
ny = psbc.shape[0]
ap,bp = np.loadtxt("apbp.txt",unpack=True)


nz = len(ap)
aps = 0.5*(ap[0:nz-1]+ap[1:nz])
bps = 0.5*(bp[0:nz-1]+bp[1:nz])
nz = len(aps)

press = np.zeros((nz,ny))
for kk in range(nz):
 for yy in range(ny):
  press[kk,yy] = aps[kk]+bps[kk]*psbc[yy]

####################################################
print "... getting fields from file !"
ubc,xdim,ydim,zdim,tdim = pp(file=fileAP,var="u"    ,t="0,1e10",x=charx).getfd() ; etape("ubc",time0)
vbc                     = pp(file=fileAP,var="v"    ,t="0,1e10",x=charx).getf()  ; etape("vbc",time0)
#tbc                     = pp(file=fileAP,var=vartemp,t="0,1e10",x=charx).getf()  ; etape("tbc",time0)

####################################################
etape("creating the target file",time0)
f = nc.Dataset(outfile,'w',format='NETCDF3_CLASSIC')
nam4 = ('press','latitude')
fie4 = (targetp1d,ydim)
for iii in range(len(nam4)):
  f.createDimension(nam4[iii],fie4[iii].shape[0])
  var = f.createVariable(nam4[iii], 'd', nam4[iii])
  var[:] = fie4[iii]
  var = None
f.close()
####################################################

ustart = pp(file=fileAP,var="u"    ,t=0   ,x=charx).getf()
uend   = pp(file=fileAP,var="u"    ,t=1e10,x=charx).getf()
dudt = (uend - ustart) / (1000.*38052.)



##
up  = pp(file=fileAP,var="u"    ,compute="pert_t").getf() ; etape("up",time0)
print up.shape
vp  = pp(file=fileAP,var="v"    ,compute="pert_t").getf() ; etape("vp",time0)
upvpb  = ppcompute.mean(up*vp ,axis=0) ; etape("upvpb" ,time0) ; del up ; del vp
print upvpb.shape
upvpbc = ppcompute.mean(upvpb ,axis=2) ; etape("upvpbc",time0) ; del upvpb

### stationary: negligible
#us  = pp(file=fileAP,var="u"    ,compute="pert_x").getf() ; etape("us",time0)
#vs  = pp(file=fileAP,var="v"    ,compute="pert_x").getf() ; etape("vs",time0)
#us_b  = ppcompute.mean(us    ,axis=0) ; etape("usb"  ,time0) ; del us
#vs_b  = ppcompute.mean(vs    ,axis=0) ; etape("vsb"  ,time0) ; del vs
#usvsbc = ppcompute.mean(us_b,axis=2)*ppcompute.mean(vs_b,axis=2) ; etape("usvsbc",time0)
#usvsbc = interpolate(targetp1d,press,usvsbc) ; etape("usvsbc",time0)


####################################################
print "... interpolating !"
ubc    = interpolate(targetp1d,press,ubc)    ; etape("ubc"   ,time0) ; addvar(outfile,nam4,'ubc',ubc)
vbc    = interpolate(targetp1d,press,vbc)    ; etape("vbc"   ,time0) ; addvar(outfile,nam4,'vbc',vbc)
#tbc    = interpolate(targetp1d,press,tbc)    ; etape("tbc"   ,time0)
upvpbc = interpolate(targetp1d,press,upvpbc) ; etape("upvpbc",time0)
dudt = interpolate(targetp1d,press,dudt)    ; etape("dudt"   ,time0) ; addvar(outfile,nam4,'dudt',dudt)

####################################################
nz,nlat = ubc.shape
nlon = 1
lat2d = np.tile(ydim,(nz,1))
acosphi2d = myp.acosphi(lat=lat2d)
cosphi2d = acosphi2d / myp.a
latrad,lat2drad = ydim*np.pi/180.,lat2d*np.pi/180.
beta = myp.beta(lat=lat2d)
f = myp.fcoriolis(lat=lat2d)
print f
tanphia = myp.tanphia(lat=lat2d)
etape("coordinates",time0)

####################################################
dubc_dy,dummy = ppcompute.deriv2d(ubc*cosphi2d,latrad,targetp1d) / acosphi2d
dummy,dubc_dp = ppcompute.deriv2d(ubc,latrad,targetp1d)
dupvpbc_dy,dummy = ppcompute.deriv2d(upvpbc*cosphi2d,latrad,targetp1d) / acosphi2d
#dusvsbc_dy,dummy = ppcompute.deriv2d(usvsbc*cosphi2d,latrad,targetp1d) / acosphi2d

#dubc_dy,dubc_dp = ppcompute.deriv2d(ubc,myp.a*latrad,targetp1d) 
#dupvpbc_dy,dummy = ppcompute.deriv2d(upvpbc,myp.a*latrad,targetp1d) 
##dusvsbc_dy,dummy = ppcompute.deriv2d(usvsbc,myp.a*latrad,targetp1d) 


vdubc_dy = -vbc*dubc_dy
fvbc = f*vbc
trans = - dupvpbc_dy
#stat = - dusvsbc_dy

summ = trans + vdubc_dy + fvbc #+ stat




# *** MASS STREAMFUNCTION ***
# *** AND VERTICAL VELOCITY ***
# NB: slide 7 in https://atmos.washington.edu/~dennis/501/501_Gen_Circ_Atmos.pdf
# term = 2 pi a cosphi / g
# --> PSIM = term * int_0^p v dp
import scipy
import scipy.integrate
psim = np.zeros((nz,nlat)) # mass streamfunction
alph = 2.*np.pi*acosphi2d/myp.g
w = np.isnan(vbc) # save NaN locations 
vbc[w] = 0. # necessary otherwise integrations fail
# integration loop
x = targetp1d[:] # coordinate
#x = np.insert(x,0,0) # JL: integrate from p=0 towards p
x = np.append(x,0) # JL: integrate from p=0 towards p
for yyy in range(nlat):
  y = vbc[:,yyy] # integrand
  #y = np.insert(y,0,y[0]) # JL: integrate from p=0 towards p
  y = np.append(y,y[-1]) # JL: integrate from p=0 towards p
  for zzz in range(0,nz):
#     the minus sign below comes from the fact that x is ordered by decreasing values of p
#          whereas the integral should be performed from 0 to p. 
    psim[zzz,yyy] = -scipy.integrate.simps(y[zzz:],x[zzz:])*alph[0,yyy]
    #psim[zzz,yyy] = scipy.integrate.simps(y[0:zzz+1],x[0:zzz+1])*alph[0,yyy]
etape("streamfunction",time0)
# reset to NaN after integration
vbc[w] = np.nan ; psim[w] = np.nan
# derivatives of streamfunction --> velocity (notably omega)
dpsim_dphi,dpsim_dp = ppcompute.deriv2d(psim,latrad,targetp1d)/alph
# meridional: v = 1/term dPSIM/dp
vphi = dpsim_dp
# vertical: omega = (-1/a) 1/term dPSIM/dphi
omegabc = -dpsim_dphi/myp.a



omegabcdubc_dp = - omegabc*dubc_dp
summ = summ + omegabcdubc_dp




#####################################################
etape("3D variables",time0)
addvar(outfile,nam4,'fvbc',fvbc)
addvar(outfile,nam4,'vdubc_dy',vdubc_dy)
addvar(outfile,nam4,'trans',trans)
#addvar(outfile,nam4,'stat',stat)
addvar(outfile,nam4,'summ',summ)
addvar(outfile,nam4,'omegabc',omegabc)
addvar(outfile,nam4,'vphi',vphi)
addvar(outfile,nam4,'omegabcdubc_dp',omegabcdubc_dp)
addvar(outfile,nam4,'psim',psim)



etape("",time0)

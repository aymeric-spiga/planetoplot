#! /usr/bin/env python

from ppclass import pp
import ppplot
import ppcompute
import planets
import numpy as np

###############
planet=planets.Saturn
fileAP="DRAG90days_DISSIP50000_lat10_913286052-951338051_red.nc"
fileAP="DRAG90days_DISSIP15000_2130950052_512.nc"
#planet=planets.Mars
#fileAP="diagfired.nc"
#fileAP="/media/WORK/SIMULATIONS/DRAG90days_DISSIP15000_2130950052-2169002051_512_z5.nc"
fileAP="/media/WORK/SIMULATIONS/DRAG90days_DISSIP15000_2130950052_256.nc"
#fileAP="DRAG90days_DISSIP15000_2130950052_512_z5.nc"
#fileAP="DRAG90days_DISSIP15000_14000-to50000-each1000.nc"
fileAP="DRAG90days_DISSIP15000_1902638052_z5.nc"
fileAP="DRAG90days_DISSIP15000_1902638052.nc"
#fileAP="UHD_DRAG18days_DISSIP5000.nc"
#fileAP="DRAG90days_DISSIP50000_lat10_z5.nc"
###############
# -- 1: meridional eddy momentum flux --> mass:N metric:N
# -- 2: meridional eddy angular momentum flux per unit mass --> mass:N metric:Y
# -- 3: meridional eddy angular momentum flux --> mass:Y metric:Y
typevar=3
###############

print "get full fields --> q"
verb = False
u4D,longit,latit,pniv,time=pp(file=fileAP,var="u",verbose=verb).getfd()
v4D=pp(file=fileAP,var="v",verbose=verb).getf()
print "   -- got shape",u4D.shape
lat=latit[:,0] # for further use

print "get axis characteristics"
timeaxis = (time.size > 1)
vertaxis = (pniv.size > 1)

print "get zonal anomaly fields --> q*=q-[q]"
staru4D=pp(file=fileAP,var="u",verbose=verb,compute="pert_x").getf()
starv4D=pp(file=fileAP,var="v",verbose=verb,compute="pert_x").getf()
print "   -- got shape",staru4D.shape

print "get temporal anomaly fields --> q'=q-qbar"
if timeaxis:
  primu4D=pp(file=fileAP,var="u",verbose=verb,compute="pert_t").getf()
  primv4D=pp(file=fileAP,var="v",verbose=verb,compute="pert_t").getf()
  print "   -- got shape",primu4D.shape
else:
  primu4D=staru4D
  primv4D=starv4D
  print "   -- no time axis: ASSUMING q' ~ q*"

print "get temporal+zonal mean fields --> [qbar]"
meanu2D=pp(file=fileAP,var="u",verbose=verb,t="0,1e15",x="-180,180").getf()
meanv2D=pp(file=fileAP,var="v",verbose=verb,t="0,1e15",x="-180,180").getf()
meant2D=pp(file=fileAP,var="temp",verbose=verb,t="0,1e15",x="-180,180").getf()
#meanmass2D=pp(file=fileAPM,var="dmass",verbose=verb,t="0,1e15",x="-180,180").getf()
print "   -- got shape",meanu2D.shape

##########################################
massmetric = 1.
# -- 1: meridional eddy momentum flux --> mass:N metric:N
# -- 2: meridional eddy angular momentum flux per unit mass --> mass:N metric:Y
# -- 3: meridional eddy angular momentum flux --> mass:Y metric:Y
##########################################
if typevar > 1:
  # compute metric terms (distance to axis)
  if vertaxis:
    acoslat2D=np.empty_like(meanv2D)
    for i in range(len(pniv)):
      # NB: thin shell approximation
      acoslat2D[i,:]=planet.acosphi(lat)
  else:
    acoslat2D=planet.acosphi(lat)
  massmetric = massmetric*acoslat2D
##########################################
if typevar > 2:
  if vertaxis:
    ## compute mass for each considered grid mesh
    ## S dp = - rho g dz S = - g dm --> dm = - S dp / g
    meanmass2D = np.empty_like(meanv2D)
    nl = len(pniv)
    dp = pniv[0:nl-2]-pniv[1:nl-1]
    dlat,dlon = latit[1,0]-latit[0,0],longit[0,1]-longit[0,0]
    surf = dlat*planet.deglength() * dlon*planet.deglength(lat)
    for i in range(len(lat)):
      meanmass2D[0:nl-2,i]=surf[i]*(pniv[0:nl-2]-pniv[1:nl-1])/planet.g
      meanmass2D[nl-1,i]=surf[i]*pniv[nl-1]/planet.g
    massmetric = massmetric*meanmass2D
  else:
    print "I need a vertical dimension to compute mass! Change file or choose typevar<3"
    exit()
##########################################

print "compute transport"
# ---- e.g. Lebonnois et al. 2010 equation 19
# [(vq)bar] : total meridional transport
# [qbar][vbar]: by mean meridional circulation
# [(qstar)bar][(vstar)bar] : by stationary waves
# [(q'v')bar] : by transients (eddies) = [(vq)bar] - [qbar vbar]
# -----
# 0. the easy bit (MMC)
mmc_trans = meanu2D*meanv2D
# 1. compute bar (temporal mean)
tot_trans = v4D*u4D
eddy_trans = primv4D*primu4D
if timeaxis:
  tot_trans = ppcompute.mean(tot_trans,axis=0)
  eddy_trans = ppcompute.mean(eddy_trans,axis=0)
  statu = ppcompute.mean(staru4D,axis=0)
  statv = ppcompute.mean(starv4D,axis=0)
# 2. compute [] (zonal mean)
if vertaxis: zeaxis=2
else: zeaxis=1
tot_trans = ppcompute.mean(tot_trans,axis=zeaxis)
eddy_trans = ppcompute.mean(eddy_trans,axis=zeaxis)
if timeaxis:
  stat_trans = ppcompute.mean(statu,axis=zeaxis)*ppcompute.mean(statv,axis=zeaxis)
else:
  stat_trans = None
# 3. add mass/metric term (computed above)
tot_trans = tot_trans*massmetric
mmc_trans = mmc_trans*massmetric
eddy_trans = eddy_trans*massmetric
if stat_trans is not None:
  stat_trans = stat_trans*massmetric
# 4. compute meridional divergence of eddy momentum flux
if not vertaxis:
  dETdy = ppcompute.deriv1d(eddy_trans,coord=lat*np.pi/180.)/planet.a  
else:
  dETdy,dETdp = ppcompute.deriv2d(eddy_trans,lat*np.pi/180.,pniv,fac=planet.a)



### PLOTS
### -- \langle q \rangle = \left[ \overline{q} \right]
print "plots"

### 1D plots
if not vertaxis:

 #####################################################################
 #####################################################################
 ## PLOT : check d <u>/dt = -d <u'v'> /dy
 pl = ppplot.plot1d()
 # get highest value and highest power
 absmax = np.abs(meanu2D*massmetric).mean()
 exponent=int(round(np.log10(absmax)))
 norm = 10.**exponent
 ## curve 1 : zonal flow
 pl.f = meanu2D*massmetric/norm
 pl.x = lat
 pl.marker = ''
 pl.xmin = -90.
 pl.xmax = 90.
 pl.xlabel = "Latitude"
 pl.linestyle = '--'
 pl.fmt = "%.1f"
 if typevar == 1:
   pl.legend = r'$\langle u \rangle$'
 elif typevar == 2:
   pl.legend = r'$\langle m_u \rangle$'
 pl.make()
 ## curve 2 : convergence of eddy momentum flux
 smoothwindow = 1 # no smooth
 smoothwindow = 30
 duration = planet.day*1000.
 pl.f = ppcompute.smooth1d(-dETdy*duration,window=smoothwindow)/norm
 pl.x = ppcompute.smooth1d(lat,window=smoothwindow)
 pl.linestyle = '-'
 if typevar == 1:
   pl.legend = r'$- \Delta t \, \frac{\partial}{\partial y} \langle v^\prime u^\prime \rangle$'
   pl.ylabel = r'$10^{%i}$ m s$^{-1}$' % (exponent)
 elif typevar == 2:
   pl.legend = r'$- \Delta t \, \frac{\partial}{\partial y} \langle v^\prime m_u^\prime \rangle$'
   pl.ylabel = r'$10^{%i}$ m$^2$ s$^{-1}$' % (exponent)
 pl.makeshow()
 #####################################################################
 #####################################################################

 #####################################################################
 #####################################################################
 ## PLOT : show the different contributions to meridional transport
 pl = ppplot.plot1d()
 # get highest value and highest power
 absmax = np.abs(tot_trans).mean()
 exponent=int(round(np.log10(absmax)))
 norm = 10.**exponent
 ##
 pl.f = tot_trans/norm
 pl.x = lat
 pl.xmin = -90.
 pl.xmax = 90.
 pl.marker = ''
 pl.xlabel = "Latitude"
 pl.fmt = "%.1f"
 if typevar == 1:
   pl.legend = r'total $\langle v \, u \rangle$'
   pl.ylabel = r'$10^{%i}$ m$^2$ s$^{-2}$' % (exponent)
 elif typevar == 2:
   pl.legend = r'total $\langle v \, m_u \rangle$'
   pl.ylabel = r'$10^{%i}$ m$^3$ s$^{-2}$' % (exponent)
 pl.make()
 ##
 pl.f = mmc_trans/norm
 if typevar == 1:
   pl.legend = r'MMC $\langle v \rangle \langle u \rangle$'
 elif typevar == 2:
   pl.legend = r'MMC $\langle v \rangle \langle m_u \rangle$'
 pl.make()
 ##
 #if stat_trans is not None:
 # pl.f = mmc_trans
 # pl.legend = r'stationary $\langle u^{\star} \rangle \langle v^{\star} \rangle$'
 # pl.make()
 ##
 pl.f = eddy_trans/norm
 if typevar == 1:
   pl.legend = r'eddies $\langle v^{\prime} u^{\prime} \rangle$'
 elif typevar == 2:
   pl.legend = r'eddies $\langle v^{\prime} m_u^{\prime} \rangle$'
 pl.makeshow()
 #####################################################################
 #####################################################################

### 2D plots
else:

 print "make plot"

 what_I_plot = dETdy
 what_I_plot = eddy_trans

 # get highest value and highest power
 absmax = np.abs(what_I_plot).mean()
 exponent=int(round(np.log10(absmax)))
 norm = 10.**exponent

 fig = ppplot.figuref(x=20,y=6)
 subv,subh = ppplot.definesubplot(2, fig)

 pl = ppplot.plot2d()
 pl.invert = True
 pl.logy = True
 pl.fmt = "%.1e"
 pl.vmin = -2.
 pl.vmax = +2.
 pl.colorbar = "seismic" #"RdBu_r"
 pl.ylabel="Pressure (Pa)"
 pl.xlabel="Latitude"
 pl.fmt = "%.1f"

 if 1 == 1:
  fig.add_subplot(subv,subh,1)
  pl.f = mmc_trans / norm
  pl.x = lat
  pl.y = pniv
  pl.units = '$10^{'+str(exponent)+'}$ s$^{-1}$'
  pl.title = "MMC transport"
  pl.make()
    
  fig.add_subplot(subv,subh,2)

 pl.f = what_I_plot / norm
 pl.c = meanu2D
 pl.x = lat
 pl.y = pniv
 pl.title = "Eddy transport"
 #pl.xmax = 60.
 #pl.xmin = 0.
 pl.make()
 ppplot.show()
 #pl.makesave(mode="png",filename="transport",includedate=False)

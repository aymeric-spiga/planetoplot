#! /usr/bin/env python
from ppclass import pp


### Active Clouds
### -------------
act = pp()
act.file = '/planeto/tnalmd/runs/newref/MY26.ifort.modelprecise/monthlybox/concat_LT14_P.nc'
act.var = 'temp'
act.x = "-180,180"
act.y = "-20,20"
act.z = 50.
###### Change time axis from Sols to Ls:
act.changetime="mars_sol2ls"
act.label="GCM active clouds"
act.title="Equatorial daytime temperature at 50 Pa"
act.marker='-'
act.color='g'
act.out = "png"
act.get()


### Inactive Clouds
### ---------------
inac = pp()
inac << act
inac.file = '/planeto/tnalmd/runs/newref/MY26.inac/monthlybox/concat_LT14_P.nc'
inac.label="GCM inactive clouds"
inac.marker=''
inac.color='r'
inac.get()


### TES data
### --------
tes = pp()
tes << act
tes.file = '/d5/emlmd/TES/TES.MappedClimatology.nadir.MY26.nc'
tes.var='T_nadir_day'
#### Do nothing with time axis (is usually the default option):
tes.changetime=None
tes.label="TES"
tes.marker='.'
tes.color='b'
tes.get()



### PLOT
### ----
act.superpose  = True
inac.superpose = True
tes.superpose  = True
### add two plots to act:
act.plot(extraplot=2)
### plot inac:
inac.plotin = act
inac.plot()
### plot tes:
tes.plotin = act
tes.plot()

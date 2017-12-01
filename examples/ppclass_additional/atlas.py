#! /usr/bin/env python

from ppclass import pp

######
gcm = pp()
gcm.file = "amipYEAR/IPSL-CM5A-MR.nc"
gcm.var = "rlut"
gcm.t = 1
#
gcm.out = "pdf"
gcm.filename = "atlas"
gcm.back = "coast"
gcm.proj = "robin"
gcm.fmt = "%.0f"
gcm.div = 20
gcm.vmin = 100
gcm.vmax = 300
gcm.colorbar = "hot"
gcm.title = "OLR (GCM)"
gcm.units = r'W m$^{-2}$'
#
gcm.getplot(extraplot=3)

######
obs = pp()
obs << gcm
obs.file = "amipYEAR/OBS.nc"
obs.plotin = gcm
obs.title = "OLR (OBS)"
obs.getplot()

######
diff = gcm - obs
diff.plotin = gcm
diff.vmin = -35
diff.vmax = 35
diff.colorbar = 'RdBu_r'
diff.title = "OLR (GCM-OBS)"
diff.plot()

######
gcmoy = pp()
gcmoy << gcm
gcmoy.x = "-180,180"
gcmoy.get()
#
obmoy = pp()
obmoy << obs
obmoy.x = "-180,180"
obmoy.get()
#
diff = gcmoy - obmoy
diff.plotin = gcm
diff.ylabel = r'OLR (W m$^{-2}$)'
diff.title = "GCM-OBS difference of zonal-mean OLR"
diff.plot()




exit()


####################################################
####################################################


#use "/d1/hourdin/ATLAS/MERGE/amipYEAR/IPSL-CM5A-MR.nc"
#use "/d1/hourdin/ATLAS/MERGE/amipYEAR/OBS.nc"
#let asr=rsdt-rsut
#let swcre=rsutcs-rsut
#let lwcre=rlutcs-rlut
#set v ul
#fill/l=1/lev=(-Inf)(0,400,10)(Inf)/x=-180:180/d=1/title="ASR W/mm2,  IPSL-CM5A-MR" asr ; go land
#set v ur
#fill/l=1/lev=(-Inf)(0,400,10)(Inf)/x=-180:180/d=2/title="ASR W/mm2, OBS" asr ; go land
#set v ll
#fill/l=1/lev=(-Inf)(-40,40,5)(Inf)/x=-180:180/d=2/title="DIFF W/m2 "/pal=blue_darkred asr[d=1]-asr[d=2] ; go land
#set v lr
#plot/l=1/vlim=50.:400/d=2 asr[d=2,i=@ave]
#plot/l=1/o asr[d=1,i=@ave]
#quit



from ppclass import pp


gcm1 = pp()
gcm1.file = "amipYEAR/IPSL-CM5A-MR.nc"
gcm1.var = "rsdt"
gcm1.x = "-180,180"
gcm1.z = 1
gcm1.t = 1
gcm1.get()

gcm2 = pp()
gcm2 << gcm1
gcm2.var = "rsut"
gcm2.get()

gcm = gcm2 - gcm1
gcm.plot()


obs1 = pp()
obs1 << gcm1
obs1.file = "amipYEAR/OBS.nc"
obs1.get()

obs2 = pp()
obs2 << gcm2
obs2.file = "amipYEAR/OBS.nc"
obs2.get()

obs = obs2 - obs1
obs.plot()




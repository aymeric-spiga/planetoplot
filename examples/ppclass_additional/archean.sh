#! /usr/bin/env python
import ppclass
import numpy as np

asr = ppclass.pp()
asr.file = "diagfi115.nc"
asr.var = "ASR"
asr.t = 360
asr.get()

isr = ppclass.pp()
isr << asr
isr.var = "ISR"
isr.get()

tau = ppclass.pp()
tau << asr
tau.var = "tau_col"
tau.get()

alb = ppclass.pp()
alb << asr
alb.var = "ALB"
alb.get()

tsurf = ppclass.pp()
tsurf << asr
tsurf.var = "tsurf"
add = tsurf.getf()



##alb.smooth(20)
#tau.smooth(2)


#albedo = asr/isr
##albedo = 1.-albedo
#albedo = albedo + taucol

albedo = alb + tau/100.

albedo.colorb = "Blues_r"
albedo.div = 20 #10
#albedo.vmin = 0.3
#albedo.vmax = 0.7
albedo.proj = "ortho"


albedo.defineplot()
albedo.p[0].blat = 15.
albedo.p[0].showcb = False
albedo.p[0].title = ""

albedo.p[0].addcontour = add - 273.15

#dafield = albedo.p[0].field
#dafield[np.where(dafield < 0.4)] = 0.0
#dafield[np.where(dafield > 0.6)] = 1.0
albedo.p[0].vmin = 0.1
albedo.p[0].vmax = 0.6

#albedo.p[0].back = "coast"

albedo.makeplot()

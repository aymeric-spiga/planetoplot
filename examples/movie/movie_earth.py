#! /usr/bin/env python
from ppclass import pp

thres = 0.1

ww = pp()
ww.xp = 24
ww.yp = 12
ww.file = "tata.nc"
ww.var = "E366"
ww.useindex = "1000"
ww.z = 0
ww.proj = "cyl"
ww.out = "jpg" #"png"
ww.res = 100 #150 #300 #75
ww.fmt = "%.0f"
ww.colorbar = "gist_gray"
ww.quiet = True
ww.showcb = False
ww.div = 50 #10 #30
ww.vmin = -thres
ww.vmax = 1
#ww.trans = 0.5
ww.back = "blue_local"
ww.nopickle = True

# loop
nn=360
nn=25
tabcase = ["full","np","sp","trop"]
tabcase = ["full"]
for case in tabcase:
 ##
 if case == "full": ww.proj = "cyl"
 elif case == "np": ww.proj = "npstere" ; ww.blat = 50
 elif case == "sp": ww.proj = "spstere" ; ww.blat = -50
 elif case == "trop": ww.proj = "cyl" ; ww.ymin = -30. ; ww.ymax = 30.
 ##
 count = 1
 for ttt in range(nn):
  print(count)
  ww.t = ttt
  ############################
  ff = pp(file=ww.file,var="TIME_COUNTER",useindex="1111",x=0,y=0,z=0,t=ttt).getf()
  day = int(ff/86400.)
  hour = int((ff - day*86400.)/3600.)
  day = day + 1
  ww.title = "day %02i @ %02i00 UTC" % (day,hour)
  ww.title = "July %i" % (day)
  print(ww.title)
  ############################
  ww.filename = "mov_"+case+"_%03d"%(count)
  ww.get()

  import numpy as np
  w = np.where(ww.f < thres)
  ww.f[w] = np.nan  

  ww.plot()
  count = count+1


## animate
## -- mencoder -mc 0 -noskip -skiplimit 0 -ovc lavc -lavcopts vcodec=mpeg4:vhq:trell:mbd=2:vmax_b_frames=1:v4mv:vb_strategy=0:cmp=6:subcmp=6:precmp=6:predia=3:dia=3:vme=4:vqscale=1 "mf://mov_full*.png" -mf type=png:fps=18 -o mov_full.avi

## vlelim=0:vcelim=0:

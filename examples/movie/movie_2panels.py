#! /usr/bin/env python
from ppclass import pp

ww = pp()
ww.file = "series30m_1km.nc"
ww.var = ["HGT","Um","Vm","TSURF","Um","Vm"]
ww.vargoal = ["main","vector","vector","main","vector","vector"]
ww.z = 0.
#ww.vmin = [-6500.,150.] #does not work
#ww.vmax = [-1500.,300.] #does not work
ww.proj = "npstere"
ww.blat = 67.
#ww.title = ["topo+winds","surftemp+winds"] #does not work
ww.out = "png"
ww.includedate = False
ww.res = 75
ww.wscale = 10.
ww.svx = 3
ww.svy = 3
ww.quiet = True

# loop
count = 0
for ttt in range(474):
  print count
  ww.filename = "mov"+"%03d"%(count)
  ww.t = ttt
  ww.getdefineplot()
  ww.p[0].vmin = -6500. 
  ww.p[0].vmax = -1500.
  ww.p[1].vmin = 140.
  ww.p[1].vmax = 290.
  ww.p[0].title = "Topography"
  ww.p[1].title = "Surface temperature"  
  ww.makeplot()
  count = count+1


## animate
## -- mencoder -mc 0 -noskip -skiplimit 0 -ovc lavc -lavcopts vcodec=mpeg4:vhq:trell:mbd=2:vmax_b_frames=1:v4mv:vb_strategy=0:vlelim=0:vcelim=0:cmp=6:subcmp=6:precmp=6:predia=3:dia=3:vme=4:vqscale=1 "mf://*.png" -mf type=png:fps=18 -o output.avi

#! /usr/bin/env python
from ppclass import pp

############################################
folder = "/home/aymeric/Remote/"
suffix = ".nc"
############################################
ffftab = ["diagfired2",\
          "saturn_128x96x64_guided_hyperdiff_20000j",\
          "saturn_128x96x64_guided_hyperdiff_40000j",\
          "saturn_128x96x64_guided_hyperdiff_60000j"]
tabtabt = [[180.,270.,360.],\
           [90.,135.,180.],\
           [90.,135.,180.],\
           [90.,135.,180.]]
############################################

motif = r"$L_s = %.0f ^{\circ}$"

n = 0

for fff in ffftab:

   tabt = tabtabt[n]
   n=n+1

   print fff
   fp = pp()
   #fp.quiet = True
   fp.verbose = True
   fp.file = folder + fff + suffix
   fp.var = "temp"
   fp.changetime = "correctls_noadd"
   fp.x = 0.
   fp.logy = True
   fp.invert = True
   fp.colorbar = "gist_rainbow_r"
   fp.xmin = -45.
   fp.xmax = 45.
   fp.ycoeff = 0.01
   fp.ymax = 50.
   fp.ymin = 2.e-3
   fp.vmin = 85.
   fp.vmax = 160.
   fp.div = 30
   fp.fmt = "%.0f"
   fp.ylabel = "Pressure (mbar)"
   fp.out = "png"
   fp.filename = "fouchetplot_"+fff
   fp.xp = 22
   fp.t = tabt[0]
   fp.title = motif % (fp.t)
   fp.getplot(extraplot=2)
   
   fp2 = pp()
   fp2 << fp
   fp2.plotin = fp
   fp2.t = tabt[1]
   fp2.title = motif % (fp2.t)
   fp2.ylabel = ""
   fp2.getplot()
   
   fp3 = pp()
   fp3 << fp
   fp3.plotin = fp
   fp3.t = tabt[2]
   fp3.title = motif % (fp3.t)
   fp3.ylabel = ""
   fp3.getplot()

#! /usr/bin/env python
from ppclass import pp

fi = "diagfi.nc"
fi = "/home/aymeric/Remote/saturn_128x96x64_guided_hyperdiff_60000j_last.nc"
fi = "/home/aymeric/Remote/saturn_128x96x64_guided_last.nc"
tt = 90
tt = 1000

mm = pp()
mm.quiet = True
mm.file = fi
mm.var = "temp"
mm.t = tt
mm.z = 1.e4
mm.div = 20
mm.proj = "robin" #"ortho"
mm.blat = 15.
mm.vmin = 70.
mm.vmax = 100.
mm.out = "png"
mm.filename = "t_100mb"
mm.fmt = "%.0f"
mm.title = ""
mm.getplot()

mm.z = 1.e3
mm.vmin = 105.
mm.vmax = 125.
mm.div = 20
mm.filename = "t_10mb"
mm.getplot()

mm.var = "u"
mm.z = 1.e4
mm.vmin = -100.
mm.vmax = 500.
mm.div = 30
mm.colorbar = "gist_ncar"
mm.filename = "u_100mb"
mm.getplot()

mm.z = 1.e3
mm.filename = "u_10mb"
mm.getplot()

pe = pp()
pe.quiet = True
pe.file = fi
pe.var = "temp"
pe.t = tt
pe.z = 1.e4
pe.div = 20
pe.compute = "pert_x"
pe.proj = "robin" #"ortho"
pe.blat = 15.
pe.colorbar = "RdBu_r"
pe.vmin = -1.
pe.vmax = 1.
pe.out = "png"
pe.filename = "t_100mb_pert"
pe.title = ""
pe.getplot()

pe.z = 1.e3
pe.vmin = -0.5
pe.vmax = 0.5
pe.filename = "t_10mb_pert"
pe.getplot()

pe.z = 1.e4
pe.var = "u"
pe.vmin = -5.
pe.vmax = 5.
pe.filename = "u_100mb_pert"
pe.getplot()

pe.z = 1.e3
pe.vmin = -2.
pe.vmax = 2.
pe.filename = "u_10mb_pert"
pe.getplot()

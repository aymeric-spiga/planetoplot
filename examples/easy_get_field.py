#! /usr/bin/env python
from ppclass import pp
import numpy as np
import matplotlib.pyplot as mpl

test = pp()
test.file = [ \
  "/planeto/jllmd/Earth/Nmix/F1525/output/diagfi4.nc", \
  "/planeto/jllmd/Earth/Nmix/F1550/output/diagfi4.nc", \
  "/planeto/jllmd/Earth/Nmix/F1500/output/diagfi4.nc" ]
test.var = ["tsurf","ASR"]
test.x = "-180.,175."
test.y = "-90.,90."
test.t = "0,1000"
#test.verbose = True
test.compute = "meanarea"
test.get()

# -----------------------------
# self.allfield contains data !
# -----------------------------
# index ordering: 
# - file
# - var
# - request t
# - request z
# - request y
# - request x
# - field dimensions
# -----------------------------
# ALL DIM == 1 HAVE BEEN SQUEEZED
# -----------------------------

# for instance here, shape of allfield is (3,2)
# because the user requested 3 files, 2 vars, 0D values

tsurf = test.allfield[:,0]
asr = test.allfield[:,1]

mpl.plot(tsurf,asr)
mpl.show()

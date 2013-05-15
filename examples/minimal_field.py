#! /usr/bin/env python
from ppclass import pp
import matplotlib.pyplot as mpl
icetot = pp(file="/home/aymeric/Big_Data/DATAPLOT/diagfired.nc",var="icetot",t=0.5).getf()
mpl.pcolor(icetot)
mpl.show()

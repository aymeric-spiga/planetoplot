#! /usr/bin/env python
from ppclass import pp
import matplotlib.pyplot as mpl
import numpy as np

icetot = pp(file="/home/aymeric/Big_Data/DATAPLOT/diagfired.nc",var="icetot",t="0.5,2.5",compute="nothing").getf()
mpl.hist(np.ravel(icetot))
mpl.show()

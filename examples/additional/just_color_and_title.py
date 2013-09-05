#! /usr/bin/env python
from ppclass import pp

anewversion = pp()
anewversion.defineplot(loadfile="../../demo_data/gw.ppobj")
for plot in anewversion.p:
    plot.colorb = "Paired"
    plot.title = "Psychedelic plot!"
anewversion.makeplot()

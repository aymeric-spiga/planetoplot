#! /usr/bin/env python
from ppclass import pp

# we define a "planetoplot" object
u = pp()

# we define its attributes
u.file = "/home/aymeric/Big_Data/DATAPLOT/diagfired.nc"
u.var = "u"
u.t = "0.5" #NB: also works without quotes
u.z = "10" #NB: also works without quotes

# we get data
u.get()

## we smooth the field a little bit
#u.smooth(15)

# we define the plot, then set a few personal stuff
u.defineplot() 
u.p[0].title = "This is what we name $u$ (m s$^{-1}$)"
u.p[0].proj = "robin"
u.filename = "simple"

# we plot
u.makeplot()

# we simply change the colorbar
# ... no need to reload data
u.p[0].colorb = "RdBu"
u.filename = "myplot"
u.makeplot()

# we remove map projection
# ... idem, no need to reload data
# ... but have to redefine plot
u.noproj = True
u.defineplot()
u.makeplot()

# we multiply the field by two
# ... and redefine+remake the plot
u = u * 2.
u.defineplot()
u.makeplot()

#! /usr/bin/env python
from ppclass import pp
import numpy as np

## settings
#fi = ["nest1_rad.nc","nest2_rad.nc","nest3_rad.nc"]
#fi = ["PHOENIX_water/my29_start03_top10km/wrfout_d01_9999-01-01_02:03:21"]
#fi = ["PHOENIX_water/old_to_compare_with/wrfout_d01_9999-01-01_07:11:42"]
#fi = ["wrfout_d01_9999-01-01_06:10:00"]
fi = ["wrfout_d01_9999-01-01_07:11:42_z"]
#limvmr=5.0

# get radius
rad=pp()
rad.file=fi
rad.var=["RICE","QH2O_ICE"]
rad.vargoal=["main","contour"]
#rad.x=-120.
#rad.t=4
#rad.x = 0
#rad.t = 30
rad.x = 28
rad.y = 46
#rad.changetime = "mars_meso_lt"
#rad.verbose = True
rad.filename = "radius"
rad.out = "png"
#rad.includedate = False
rad.get()

# radius in microns
rad = rad*1.e6
rad.units = "$\mu$m"
rad.fmt = "%.0f"

# define plot
rad.defineplot()

# loop on all plots
for plotobj in rad.p:

    # compute limit
    fac = 100.
    limvmr = np.mean(plotobj.c)/fac
    print("mixing ratio limit is ", limvmr)
    # mask plotted field with contoured field
    plotobj.f[np.where(plotobj.c<limvmr)]=np.nan
    # a few plot settings 
    # (some of them could have been actually set globally in rad object)
    plotobj.vmin = 0.
    plotobj.vmax = 40. 
    plotobj.div = 20
#    plotobj.xmin = 12.
#    plotobj.xmax = 16.
    plotobj.colorbar = "gist_stern"
    plotobj.ycoeff = 1.e-3
    plotobj.ylabel = 'Altitude above MOLA$_0$ (km)'
#    plotobj.xlabel = 'Latitude ($^{\circ}$N)'
#    plotobj.xcoeff = 0.2
    plotobj.xlabel = 'Time (seconds)'
    plotobj.xcoeff = 100.
#    plotobj.nxticks = 4
#    plotobj.logy = True

## set different titles
#rad.p[0].title = "Ice radius. Nest 1."
#rad.p[1].title = "Ice radius. Nest 2."
#rad.p[2].title = "Ice radius. Nest 3."


# make final plot
rad.makeplot()

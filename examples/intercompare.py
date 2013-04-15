#! /usr/bin/env python
from ppclass import pp

### LMD simulation results
### ----------------------
### set object
lmdu = pp()
### set file var coord
lmdu.file = 'LMD/wrfout_d01_2024-09-08_01:00:00_z'
lmdu.var = 'Um'
lmdu.x = -6.
lmdu.y = -2.
lmdu.t = 0.
### get data
lmdu.get()

### MRAMS simulation results
### ------------------------
### set object
mramsu = pp()
### get attributes from lmd object
mramsu << lmdu
### OK just change what changes in MRAMS
mramsu.file = 'MRAMS/mramsout_d01_2024-09-08_01:00:00_z'
mramsu.var = 'u_areo'
### get data
mramsu.get()

### THE SAME BUT FOR V
### ------------------
lmdv = pp()
lmdv << lmdu
lmdv.var = 'Vm'
lmdv.get()
mramsv = pp()
mramsv << mramsu
mramsv.var = 'v_areo'
mramsv.get()

### COMPUTE WIND SPEED
### ------------------
lmd = (lmdu**2+lmdv**2)**0.5
mrams = (mramsu**2+mramsv**2)**0.5

### NOW PLOT
### --------
### define plot for LMD
### ... and prepare it so that MRAMS will be superimposed
lmd.superpose = True
lmd.defineplot(extraplot=1)
### ... add a few personal settings
lmd.p[0].title = "LMD (blue) vs. MRAMS (red)"
lmd.p[0].swaplab = False
lmd.p[0].xlabel = "Wind speed (m s$^{-1})$"
lmd.p[0].ylabel = "Altitude above MOLA zero datum (km)"
lmd.p[0].ycoeff = 1./1000.
### ... and make plot (no output, wait for extraplot)
lmd.makeplot()
### --------
### define plot for MRAMS
### ... say we will plot this in lmd figure
mrams.plotin = lmd
mrams.superpose = True
### ... and make plot (now there is an output)
mrams.out = 'png'
mrams.filename = 'wind_intercomp'
mrams.plot()

### COMPUTE REL DIFF in % AND PLOT IT INDEPENDENTLY
### -----------------------------------------------
diff = ((lmd-mrams)/(lmd*0.5+mrams*0.5))*100.
diff.filename = 'wind_intercomp_diff'
diff.superpose = False
diff.out = 'png'
diff.defineplot()
diff.p[0].title = ""
diff.p[0].swaplab = False
diff.p[0].xlabel = "Relative difference in wind LMD vs. MRAMS (%)"
diff.p[0].ylabel = "Altitude above MOLA zero datum (km)"
diff.p[0].ycoeff = 1./1000.
diff.makeplot()

### SAME FOR TEMPERATURE
### --------------------
lmdt = pp()
lmdt << lmdu
lmdt.var = 'tk'
lmdt.get()
mramst = pp()
mramst << mramsu
mramst.var = 'tk'
mramst.get()
difft = ((lmdt-mramst)/(lmdt*0.5+mramst*0.5))*100.
difft.defineplot()
difft.filename = 'temp_intercomp_diff'
difft.superpose = False
difft.out = 'png'
difft.defineplot()
difft.p[0].title = ""
difft.p[0].swaplab = False
difft.p[0].xlabel = "Relative difference in temperature LMD vs. MRAMS (%)"
difft.p[0].ylabel = "Altitude above MOLA zero datum (km)"
difft.p[0].ycoeff = 1./1000.
difft.makeplot()

### TEMPERATURE
###############

int = pp()
int.file = ['LMD/wrfout_d01_2024-09-08_01:00:00_z','MRAMS/mramsout_d01_2024-09-08_01:00:00_z']
int.var = 'tk'
int.x = -6.
int.y = -2.
int.t = 0.
int.out = 'png'
int.filename = 'temp_intercomp'
int.superpose = True
int.getdefineplot()

int.p[0].title = "LMD (blue) vs. MRAMS (red)"
int.p[0].swaplab = False
int.p[0].ylabel = "Altitude above MOLA zero datum (km)"
int.p[0].ycoeff = 1./1000.
int.p[0].title = lmd.p[0].title
int.p[0].xlabel = "Temperature (K)"

int.makeplot()


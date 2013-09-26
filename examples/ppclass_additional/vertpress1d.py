#! /usr/bin/env python
from ppclass import pp

##################################################################
##################################################################

files = ["diagfi_start150.nc","diagfi_start175.nc","diagfi_start200.nc"]
labs = ['$T_{ini}=150$ K','$T_{ini}=175$ K','$T_{ini}=200$ K']
name = "profile_tini"

#files = ["diagfi_start200.nc","diagfi_start200_Rayleigh.nc","diagfi_start200_A034.nc"]
#labs = ['ref ($T_{ini}=200$ K)','ref + Rayleigh scatt.','ref + $A_{surf} = 0.34$']
#name = "profile_tests"

##################################################################
##################################################################

# OK retrieve pressure and temperatures
press = pp(file=files,var="p",x=0,y=0,t=1.e15).get()
temp = pp(file=files,var="temp",x=0,y=0,t=1.e15).get()

# this must be done before the next step
# ... because the next step is doing defineplot
temp.superpose = True
#temp.out = "png"
temp.filename = name
temp.includedate = False
temp.colorbar = "hot"

# this is doing defineplot amongst other things
temp.func(press)

# customize plots
count=0
for yeah in temp.p:
    yeah.logy = True
    yeah.invert = True
    yeah.swaplab = False
    yeah.ylabel = "Pressure (Pa)"
    yeah.xlabel = "Temperature (K)"
    yeah.ymin = 1.
    yeah.xmax = 230.
    yeah.xmin = 70.
    yeah.div = 32
    yeah.legend = labs[count]
    yeah.title = "Equilibrium profile with a 1D globally-averaged model after 50 Saturn years"
    count = count+1

# finally make the plot
temp.makeplot()

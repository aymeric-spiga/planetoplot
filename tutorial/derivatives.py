
# coding: utf-8

# ### Use PLANETOPLOT to compute and plot derivatives
# 
# *Author: [Aymeric SPIGA](http://www.lmd.jussieu.fr/~aslmd)*
# 
# The goal of this tutorial is to use PLANETOPLOT to compute (and plot) derivatives of a given field (say, temperature $T$).
# 
# -- Note that we make use here of the convenient `planets` library to get planetary constants. You can download and learn about [here](https://github.com/aymeric-spiga/planets)
# 
# First and foremost, import all the necessary librairies. Note that you need the `ppcompute` utility from PLANETOPLOT: derivatives functions are in this library.

# In[88]:

from ppclass import pp
import ppplot
import ppcompute
import planets
import numpy as np


# In[89]:

# This line configures matplotlib to show figures embedded in the notebook.
get_ipython().magic(u'matplotlib inline')


# Set up the planet and file to consider. Here we use the file used in the main PLANETOPLOT tutorial.

# In[90]:

myplanet = planets.Mars
fff = "diagfired.nc"


# Get the temperature field at a given time and zonally averaged

# In[91]:

temp,lon,lat,z,t = pp(file=fff,var="temp",t=0.8,x="-180,180").getfd()


# Convert degree to radian

# In[92]:

phi = lat*np.pi/180.
sinphi = np.sin(phi)
cosphi = np.cos(phi)


# Convert altitudes from km to m

# In[93]:

z = z*1000.


# Compute latitudinal & vertical derivatives using `ppcompute`

# In[94]:

dTdphi,dTdz = ppcompute.deriv2d(temp,phi,z)


# Plot with ppplot the meridional gradient of temperature as shaded field and temperature as contoured field 

# In[95]:

pl = ppplot.plot2d()
pl.f = dTdphi
pl.c = temp 
pl.x = lat
pl.y = z
pl.ylabel = 'altitude (m)'
pl.xlabel = 'latitude ($^{\circ}$N)'
pl.title = '$dT/d\phi$'
pl.makeshow()


# Compute the Brunt-Väisälä frequency
# $$ N^2 = \frac{g}{T_0} \left( \frac{-g}{c_p} + \frac{\text{d}T}{\text{d}z} \right) $$
# We use here the function defined in the `planets` library

# In[96]:

bv = myplanet.N2(T0=temp,dTdz=dTdz)


# Plot with `ppplot` the BV frequency as shaded field and temperature as contoured field. We can actually use the `plot2d()` object defined above since a bunch of settings are in common.

# In[97]:

pl.f = bv 
pl.title = '$N^2$'
pl.vmax = 2.e-4
pl.makeshow()


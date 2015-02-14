
# coding: utf-8

# ### Operations on PLANETOPLOT requests: differences, ratios, etc...
# 
# *Author: [Aymeric SPIGA](http://www.lmd.jussieu.fr/~aslmd)*

# In[45]:

# This line configures matplotlib to show figures embedded in the notebook.
get_ipython().magic(u'matplotlib inline')


# Operations on PLANETOPLOT objects `pp()` make it easy to show the difference or ratio between two PLANETOPLOT requests: this can be difference between two simulations, ratio between two variables, etc... The five operations -- plus minus multiply divide power -- are coded in the `pp()` class. Below we give several examples
# 
# This is done by setting up a `python` script or using `ipython` for interactive use. First import `pp()` class.

# In[46]:

from ppclass import pp


# Now perform two distinct requests and apply the `get()` method to load data from netCDF file(s). NB: we use the same data file as in the main tutorial.

# In[47]:

req1 = pp(file="diagfired.nc",var="tsurf",t=0.7).get()
req2 = pp(file="diagfired.nc",var="tsurf",t=0.9).get()


# Now create a new `pp()` object containing the difference between the two requests.

# In[48]:

diff = req2-req1


# It is then easy to plot the difference between the two requested fields! Simply call the `plot()` method for the `diff` object.

# In[49]:

diff.plot()


# Operations with actual numbers are also supported. For instance, show surface temperature in degrees Celsius instead of Kelvin.
# 
# Note that the computed object (here, `cels`) gets all his attributes from the original `req2` object. Before the plotting command, you might want to provide settings more suitable for `ratio`. For instance here we change units.

# In[54]:

cels = -273.15 + req1
cels.units = "$^{\circ}C$"
cels.plot()


# Computing a ratio is not the least difficult than computing additions and differences. In the example below, we changed the formatting of values to get 2 decimals since the `fmt` attribute imported from the `req2` object is not suitable (float with no decimal).

# In[44]:

ratio = req2/req1
ratio.fmt = "%.2f"
ratio.plot()


# One last example, suppose you need to plot horizontal wind modulus at a given altitude
# $$ U = \sqrt{u^2+v^2} $$

# In[30]:

u = pp(file="diagfired.nc",var="u",t=0.7,z=120.).get()
v = pp(file="diagfired.nc",var="v",t=0.7,z=120.).get()
wind = (u**2 + v**2)**0.5
wind.proj = "moll"
wind.plot()


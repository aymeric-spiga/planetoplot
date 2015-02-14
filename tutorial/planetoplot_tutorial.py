
# coding: utf-8

# # PLANETOPLOT tutorial
# 
# *Author: [Aymeric SPIGA](http://www.lmd.jussieu.fr/~aslmd)*
# 
# PLANETOPLOT is a plotting/mapping tool based on popular Python librairies. Here the great work done by the [matplotlib](http://matplotlib.org/) and [basemap](http://matplotlib.org/basemap) teams shall be acknowledged.
# 
# The tool is [available on Github](https://github.com/aymeric-spiga/planetoplot) under a GNU GPL licence. Please refer to [this page for download instructions](https://github.com/aymeric-spiga/planetoplot/wiki).
# 
# Initially I developed it to learn more about Python and to build an unified tool I could use for my research (and within our research team). I hope it will be useful to you as well. Note that I am a research scientist, not a computer programmer: this tool comes as is, without guarantee. 
# 
# Below is a quick tutorial (and a sample gallery) to discover the tool. Download [here](http://www.lmd.jussieu.fr/~aslmd/planetoplot/diagfired.nc) the file used in this tutorial (a 32M file with predictions from our Mars Global Climate Model!). 
# 
# This tutorial cover the three kinds of use of PLANETOPLOT:
# * As a quick and convenient [command line](#commandline) tool
# * As a versatile Python [library](#library) to use in your own scripts (or in `ipython`)
# * As a source of interesting libraries for your work -- in a [modular](#modular) approach
# 
# Enjoy!

# <small>NB: If you have `ipython notebook` installed on your computer, download and use this tutorial interactively with the command `ipython notebook planetoplot_tutorial.ipynb`</small>

# In[1]:

# This line configures matplotlib to show figures embedded in the notebook, 
# instead of opening a new window for each figure.  
# If you are using an old version of IPython, try using '%pylab inline' instead.
get_ipython().magic(u'matplotlib inline')


# <a id='commandline'></a>
# ## 1. Quick and convenient: Use PLANETOPLOT as a command line tool

# You can use PLANETOPLOT as a simple command line interface.
# 
# This allows you to easily and quickly plot a field saved in a netCDF file.
# 
# ### Important general remarks for a start
# 
# * First thing to do is to discover all the available options is
# 
#     `pp.py -h`
#     
#     Note that options related to reading fields are lower case while those related to plotting fields are upper case.
#     
#     
# * Assume you have a netCDF file named `diagfired.nc`, the command 
# 
#     `pp.py diagfired.nc`
# 
#     will give you some information on the available variables and on the xyzt dimensions that the program is able to recognize in the netCDF file 
#  
#  
# * The general use of the command is
# 
#     `pp.py [options] file(s)` or, equivalently, `pp.py file(s) [options]`
#     
#     To obtain the same plot side-by-side for two different files `file1.nc` and `file2.nc`
#     
#     `pp.py [options] file1.nc file2.nc`
#     
#     This works for multiples files. And, of course, regular expressions can be used, as well as automatic completion
#     
#     `pp.py [options] file*.nc`
#     `pp.py [options] file?.nc`
#     `pp.py [options] file[1-6].nc`
#     
#     
# * In any example that follows, adding the option `--verbose` make the program describe what it is doing. This can be useful to understand how PLANETOPLOT works.
# 
# 
# * In any example that follows, the command can be saved to a txt file named `your_choice.sh` by using the option `-o your_choice`. This can be useful to store long commands for further reference.
# 
# 
# * In any example that follows, the figure is output through the great `matplotlib` GUI. To save the figure in one of the usual compressed or vector format, use the corresponding option: `-O png` or `-O jpg` or `-O eps` or `-O ps` or `-O svg` or `-O pdf`
# 
# 

# ### Tutorial examples
# 
# A guided tutorial is worth a thousands words. 
# 
# Try the examples below with the file you just downloaded (if not, check out [here](http://www.lmd.jussieu.fr/~aslmd/planetoplot/diagfired.nc)) 
# 
# *NB: in what follows, type each command in a terminal without `%run` (except if you use `ipython`)*

# **Example 1: Map a time-varying 2D variable tsurf by slicing at time value 0.9**

# In[3]:

get_ipython().magic(u'run pp.py diagfired.nc -v tsurf -t 0.9')


# ** Example 1bis: Same as Example 1 except use another [map projection](http://matplotlib.org/basemap/api/basemap_api.html) **

# *Robinson projection*

# In[3]:

get_ipython().magic(u'run pp.py diagfired.nc -v tsurf -t 0.9 -P robin')


# *Orthographic projection with point of view centered in longitude 90W and latitude 30N*

# In[4]:

get_ipython().magic(u'run pp.py diagfired.nc -v tsurf -t 0.9 -P ortho --blon -90 --blat 30')


# ** Example 1ter: Same as Example 1 except use another [colormap](http://matplotlib.org/examples/color/colormaps_reference.html) **

# In[5]:

get_ipython().magic(u'run pp.py diagfired.nc -v tsurf -t 0.9 -C hot')


# **Example 2: Map a time-varying 2D variable tsurf by slicing at several time values**

# In[6]:

get_ipython().magic(u'run pp.py diagfired.nc -v tsurf -t 0.5 -t 0.6 -t 0.7 -t 0.8')


# **Example 3: Reduce dimension further by slicing at longitude value 45 in addition to time value 0.9**

# In[7]:

get_ipython().magic(u'run pp.py diagfired.nc -v tsurf -t 0.9 -x 45')


# **Example 4: Superimpose two 1D slices at different time values**

# In[8]:

get_ipython().magic(u'run pp.py diagfired.nc -v tsurf -t 0.5 -t 0.7 -t 0.9 -x 45 -S')


# The color and legend for each curve can be set by using respectively the option `-Q` and `-E`

# In[9]:

get_ipython().magic(u'run pp.py diagfired.nc -v tsurf -t 0.5 -t 0.7 -t 0.9 -x 45 -S -Q red -E afternoon -Q purple -E evening -Q black -E night')


# ** Example 5: Obtain a yz section at longitude 0 and time value 0.9 of a time-varying 3D variable u **

# In[10]:

get_ipython().magic(u'run pp.py diagfired.nc -v u -t 0.9 -x 0')


# ** Example 6: Same as Example 5 except request a zonal average between -175° and 180° longitude **

# In[11]:

get_ipython().magic(u'run pp.py diagfired.nc -v u -t 0.9 -x -175,180')


# ** Example 7: Same as example 6 except request minimum values between -175° and 180° longitude **

# In[12]:

get_ipython().magic(u'run pp.py diagfired.nc -v u -t 0.9 -x -175,180 -u min')


# ** Example 8: Same as example 6 except change min and max displayed values, number of color levels, title, x-axis title, and y-title **

# In[13]:

get_ipython().magic(u'run pp.py diagfired.nc -v u -t 0.9 -x -175,180       -N -100 -M 100 -D 50       -T "Zonal mean for zonal wind $\\langle u \\rangle$" -X \'Latitude\' -Y \'Altitude (km)\'')


# ** Example 9: Same as example 8 except request two different 3D variables to be displayed with different settings **

# In[14]:

get_ipython().magic(u'run pp.py diagfired.nc -t 0.9 -x -175,180 -X \'Latitude\' -Y \'Altitude (km)\' -D 30 -F \'%.0f\'       -v u -N -100 -M 100 -T "Zonal wind $\\langle u \\rangle$"       -v temp -N 100 -M 250 -T "Temperature $\\langle T \\rangle$"')


# Note that we also used the option `-F` which impose the format of values in the colorbar. Here we set floats without decimals, but any `python` format can be used.

# ** Example 10: Map temperature at altitude 120 km and time value 0.9 with wind vectors superimposed (two components: u and v) **

# In[15]:

get_ipython().magic(u'run pp.py diagfired.nc -v temp -z 120 -t 0.9 -i u -j v')


# <a id='library'></a>
# ## 2. Versatile and complete: Make PLANETOPLOT python scripts

# Line-by-line description
# * Script header
# * Import `pp` object from `ppclass`
# * Define a `pp` object named `ex`
#   (and make it quiet)
# * Define attributes
#   - netCDF file to read
#   - variable to read
#   - value on y axis
#   - value on z axis
# * Get data from netCDF file + Make plot

# In[16]:

#! /usr/bin/env python
from ppclass import pp
ex = pp()
ex.quiet = True
ex.file = 'diagfired.nc'
ex.var = 'temp'
ex.y = 15.
ex.z = 75.
ex.getplot()


# Here is an example on how to customize this plot

# In[17]:

ex.title = r'Hovmoller plot for $T_{75km}$'
ex.xcoeff = 24.
ex.xlabel = 'Time (hours)'
ex.ylabel = 'Altitude (km)'
ex.vmin = 130.
ex.vmax = 170.
ex.div = 50
ex.fmt = '%.0f'
ex.colorbar = "spectral"
ex.plot() # remake plot


# The first example in the [command line](#commandline) section
# 
# `pp.py diagfired.nc -v tsurf -t 0.9`
# 
# is equivalent to the following (relatively minimalist) script

# In[18]:

#! /usr/bin/env python
from ppclass import pp
mini = pp()
mini.quiet = True
mini.file = "diagfired.nc"
mini.var = "tsurf"
mini.t = 0.9
mini.xp = 4 # this line is just added to make figure below smaller
mini.yp = 3 # this line is just added to make figure below smaller
mini.getplot()


# and this script is equivalent to the following one using single-line syntax

# In[19]:

#! /usr/bin/env python
from ppclass import pp
mini = pp(quiet=True,file="diagfired.nc",var="tsurf",t=0.9,xp=4,yp=3).getplot()


# <a id='modular'></a>
# ## 3. Modular and powerful: Use individual PLANETOPLOT components in your scripts

# ### Reading netCDF only
# 
# In `ppclass`, the `pp()` class features methods, on the one hand, to retrieve fields from netCDF files and, on the other hand, to plot those fields. However, it is easy to use only the former capability and not the latter. This is often the case when your own plotting recipes are preferred, or simply that complex operations or data analysis/filtering are needed prior to displaying the fields. This is easily done with the method `getf()`.
# 
# For instance, reading the field `icetot` in the file `diagfired.nc` at time value 0.5 takes only one line

# In[20]:

#! /usr/bin/env python
from ppclass import pp
icetot = pp(file="diagfired.nc",var="icetot",t=0.5).getf()
print icetot[10,10]


# A more complex and versatile use of the getf() method -- and getfl() for labelling -- is summarized in this [example script](https://github.com/aymeric-spiga/planetoplot/blob/master/examples/ppclass_additional/easy_get_field.py)

# ### Plotting only: the `ppplot` class

# The `ppplot` class can be used without the `ppclass` class. This is for instance the case if one wants to plot fields that are not stored in netCDF format, or simply to plot anything within a `python` script with a quick and convenient solution.
# 
# In the following example script, using a set of librairies including `ppplot` 

# In[21]:

#! /usr/bin/env python
import urllib
import numpy as np
import ppplot


# a radiosounding stored somewhere on the web is imported and saved in a local file `sounding.txt`

# In[22]:

url="https://raw.githubusercontent.com/aymeric-spiga/planetoplot/master/examples/ppplot/sounding.txt"
sounding = urllib.urlopen(url).read()
soundingfile = open('sounding.txt','w').write(sounding)


# and loaded using `numpy.loadtxt()`

# In[23]:

press,z_alt,temp,dew,hum,mix,wdir,wknot,thta,thte,thtv     = np.loadtxt("sounding.txt",skiprows=8,unpack=True)


# before the variables can be easily displayed using `ppplot` through a one-line instruction

# In[24]:

ppplot.plot1d(f=temp).makeshow()


# A more elaborate plot can be obtained as well by setting more attributes for the `plot1d` class object

# In[25]:

sd = ppplot.plot1d()
sd.f = z_alt
sd.x = temp
sd.linestyle = '-'
sd.marker = '.'
sd.color = 'r'
sd.ycoeff = 1.e-3
sd.title = "A random terrestrial sounding"
sd.legend = "Fort Smith 00Z 26 Sep 2013"
sd.xlabel = "Temperature ($^{\circ}$C)"
sd.ylabel = "Altitude (km)"
sd.makeshow()


# The command line script `asciiplot.py` was created on this basis to plot something from an ASCII file containing columns. This script shares common options with the PLANETOPLOT `pp.py` command line script.

# In[26]:

get_ipython().magic(u'run asciiplot.py -h')


# Here is a command line example doing something similar to the above script

# In[27]:

get_ipython().magic(u'run asciiplot.py sounding.txt -s 8 -x 3 -y 2')


# For the exact same result as what we got with the above script, a few options can be added

# In[28]:

get_ipython().magic(u'run asciiplot.py sounding.txt -s 8 -x 3 -y 2         -L \'-\' -K \'.\' -Q \'r\' --ycoeff 1.e-3         -T "A random terrestrial sounding" -E "Fort Smith 00Z 26 Sep 2013"         -X "Temperature ($^{\\circ}$C)" -Y "Altitude (km)"')


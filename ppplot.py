###############################################
## PLANETOPLOT                               ##
## --> PPPLOT                                ##
###############################################
## Author: Aymeric Spiga. 02-03/2013         ##
###############################################
# python built-in librairies
import time as timelib
# added librairies
import numpy as np
import matplotlib.pyplot as mpl
from matplotlib.cm import get_cmap
from mpl_toolkits.basemap import Basemap
from matplotlib.ticker import FormatStrFormatter,MaxNLocator
# personal librairies
import ppcompute
###############################################

#################################
# global variables and settings #
#################################

# matplotlib settings
# http://matplotlib.org/users/customizing.html
# -------------------------------
mpl.rcParams['font.family'] = "serif"
mpl.rcParams['axes.color_cycle'] = "b,r,g,k"
mpl.rcParams['contour.negative_linestyle'] = "dashed" # ou "solid"
mpl.rcParams['verbose.level'] = "silent"
mpl.rcParams['lines.linewidth'] = 1.5
mpl.rcParams['lines.markersize'] = 10
mpl.rcParams['xtick.major.pad'] = 10
mpl.rcParams['ytick.major.pad'] = 10

# global variables
# -------------------------------
# - where settings files are located
#   (None means planetoplot_v2 in PYTHONPATH)
whereset = None
# - some good default settings.
# (bounds)
how_many_sigma = 3.0 
# (contours)
ccol = 'black'
cline = 0.55
# (vectors)
widthvec = 0.002
reducevec = 30.
# (colorbar)
zeorientation="vertical"
zefrac = 0.05
# (save figures)
pad_inches_value=0.25
# (negative mode)
def_negative = False
###############################################

### settings for 'negative-like' mode
if def_negative:
    mpl.rc('figure', facecolor='k', edgecolor='k')
    mpl.rcParams['text.color'] = 'w'
    mpl.rc('axes',labelcolor='w',facecolor='k',edgecolor='w')
    mpl.rcParams['xtick.color'] = 'w'
    mpl.rcParams['ytick.color'] = 'w'
    mpl.rcParams['grid.color'] = 'w'
    mpl.rc('savefig',facecolor='k',edgecolor='k')

##########################
# executed when imported #
##########################
###########################################
# we load user-defined automatic settings #
###########################################
# initialize the warning variable about file not present...
files_not_present = ""
# ... and the whereset variable
whereset = ppcompute.findset(whereset)

# - variable settings
# -------------------------------
zefile = "set_var.txt"
vf = {} ; vc = {} ; vl = {} ; vu = {}
try: 
    f = open(whereset+zefile, 'r')
    for line in f:
        if "#" in line: pass
        else:
            var, format, colorb, label, units = line.strip().split(';')
            ind = var.strip() 
            vf[ind] = format.strip()
            vc[ind] = colorb.strip()
            vl[ind] = label.strip()
            vu[ind] = units.strip()
    f.close()
except IOError: 
    files_not_present = files_not_present + zefile + " "

## - file settings
## -------------------------------
#zefile = "set_file.txt"
#prefix_t = {} ; suffix_t = {} ; name_t = {} ; proj_t = {} ; vecx_t = {} ; vecy_t = {}
#try:
#    f = open(whereset+zefile, 'r')
#    for line in f:
#        if "#" in line: pass
#        else:
#            prefix, suffix, name, proj, vecx, vecy = line.strip().split(';')
#            ind = name.strip() 
#            prefix_t[ind] = prefix.strip()
#            suffix_t[ind] = suffix.strip()
#            #name_t[ind] = name.strip()
#            proj_t[ind] = proj.strip()
#            vecx_t[ind] = vecx.strip()
#            vecy_t[ind] = vecy.strip()
#    f.close()
#except IOError:
#    files_not_present = files_not_present + zefile + " "

# - multiplot settings
# -------------------------------
zefile = "set_multiplot.txt"
subv_t = {} ; subh_t = {} ; wspace_t = {} ; hspace_t = {} ; font_t = {}
try:
    f = open(whereset+zefile, 'r')
    for line in f:
        if "#" in line: pass
        else:
            num, subv, subh, wspace, hspace, font = line.strip().split(';')
            ind = int(num.strip())
            subv_t[ind] = int(subv.strip())
            subh_t[ind] = int(subh.strip())
            wspace_t[ind] = float(wspace.strip())
            hspace_t[ind] = float(hspace.strip())
            font_t[ind] = float(font.strip())
    f.close()
except IOError:
    files_not_present = files_not_present + zefile + " "

# - background settings
# -------------------------------
zefile = "set_back.txt"
back = {}
try:
    f = open(whereset+zefile, 'r')
    for line in f:
        if "#" in line: pass
        else:
            name, link = line.strip().split(';')
            ind = name.strip() 
            back[ind] = link.strip()
    f.close()
except IOError:
    files_not_present = files_not_present + zefile + " "

# - area settings
# -------------------------------
zefile = "set_area.txt"
area = {}
try:
    f = open(whereset+zefile, 'r')
    for line in f:
        if "#" in line: pass
        else:
            name, wlat1, wlat2, wlon1, wlon2 = line.strip().split(';')
            area[name.strip()] = [[float(wlon1.strip()),float(wlon2.strip())],\
                                  [float(wlat1.strip()),float(wlat2.strip())]]
    f.close()
except IOError:
    files_not_present = files_not_present + zefile + " "
# A note to the user about missing files
if files_not_present != "":
    print "warning: files "+files_not_present+" not in "+whereset+" ; those presets will be missing"

# TBD: should change vector color with colormaps
#"gist_heat":    "white",\
#"hot":          "white",\
#"gray":         "red",\

####################
# useful functions #
####################

# a function to change color lines according to a color map (idea by A. Pottier)
# ... two optional arguments: color maps + a number telling how much intervals are needed
def rainbow(cb="jet",num=8):
    ax = mpl.gca()
    pal = mpl.cm.get_cmap(name=cb)
    ax._get_lines.set_color_cycle([pal(i) for i in np.linspace(0,0.9,num)])

# a function to define subplot
# ... user can change settings in set_multiplot.txt read above
# -------------------------------
def definesubplot(numplot, fig):
    try: 
        mpl.rcParams['font.size'] = font_t[numplot]
    except: 
        mpl.rcParams['font.size'] = 18
    try: 
        fig.subplots_adjust(wspace = wspace_t[numplot], hspace = hspace_t[numplot])
        subv, subh = subv_t[numplot],subh_t[numplot]
    except: 
        print "!! WARNING !! no preset found from set_multiplot.txt, or this setting file was not found."
        subv = 1 ; subh = numplot
    return subv,subh

# a function to calculate automatically bounds (or simply prescribe those)
# -------------------------------
def calculate_bounds(field,vmin=None,vmax=None):
    # prescribed cases first
    zevmin = vmin
    zevmax = vmax
    # computed cases
    if zevmin is None or zevmax is None:
       # select values
       ind = np.where(field < 9e+35)
       fieldcalc = field[ ind ] # field must be a numpy array
       # calculate stdev and mean
       dev = np.std(fieldcalc)*how_many_sigma
       damean = ppcompute.mean(fieldcalc)
       # fill min/max if needed
       if vmin is None: zevmin = damean - dev
       if vmax is None: zevmax = damean + dev
       # special case: negative values with stddev while field is positive
       if zevmin < 0. and ppcompute.min(fieldcalc) >= 0.: zevmin = 0.
    # check that bounds are not too tight given the field
    amin = ppcompute.min(field)
    amax = ppcompute.max(field)
    if np.abs(amin) < 1.e-15: 
        cmin = 0.
    else:
        cmin = 100.*np.abs((amin - zevmin)/amin)
    cmax = 100.*np.abs((amax - zevmax)/amax)
    if cmin > 150. or cmax > 150.:
        print "!! WARNING !! Bounds are a bit too tight. Might need to reconsider those."
        print "!! WARNING !! --> actual",amin,amax,"adopted",zevmin,zevmax
    return zevmin, zevmax    
    #### treat vmin = vmax for continuity 
    #if vmin == vmax:  zevmin = damean - dev ; zevmax = damean + dev

# a function to solve the problem with blank bounds !
# -------------------------------
def bounds(what_I_plot,zevmin,zevmax,miss=9e+35):
    small_enough = 1.e-7
    if zevmin < 0: what_I_plot[ what_I_plot < zevmin*(1.-small_enough) ] = zevmin*(1.-small_enough)
    else:          what_I_plot[ what_I_plot < zevmin*(1.+small_enough) ] = zevmin*(1.+small_enough)
    what_I_plot[ what_I_plot > miss  ] = -miss
    what_I_plot[ what_I_plot > zevmax ] = zevmax*(1.-small_enough)
    return what_I_plot

# a generic function to show (GUI) or save a plot (PNG,EPS,PDF,...)
# -------------------------------
def save(mode="gui",filename="plot",folder="./",includedate=True,res=150,custom=False):
    if mode != "nothing":
      # a few settings
      possiblesave = ['eps','ps','svg','png','jpg','pdf'] 
      # now the main instructions
      if mode == "gui": 
          mpl.show()
      elif mode in possiblesave:
          ## name of plot
          name = folder+'/'+filename
          if includedate:
              for ttt in timelib.gmtime():
                  name = name + "_" + str(ttt)
          name = name +"."+mode
          ## save file
          print "**** Saving file in "+mode+" format... Please wait."
          if not custom:
              # ... regular plots
              mpl.savefig(name,dpi=res,pad_inches=pad_inches_value,bbox_inches='tight')
          else:
              # ... mapping mode, adapted space for labels etc...
              mpl.savefig(name,dpi=res)
      else:
          print "!! ERROR !! File format not supported. Supported: ",possiblesave

##################################
# a generic class to make a plot #
##################################
class plot():

    # print out a help string when help is invoked on the object
    # -------------------------------
    def __repr__(self):
        whatprint = 'plot object. \"help(plot)\" for more information\n'
        return whatprint

    # default settings
    # -- user can define settings by two methods. 
    # -- 1. yeah = plot2d(title="foo")
    # -- 2. yeah = pp() ; yeah.title = "foo"
    # -------------------------------
    def __init__(self,\
                 var=None,\
                 field=None,\
                 absc=None,\
                 xlabel="",\
                 ylabel="",\
                 div=20,\
                 logx=False,\
                 logy=False,\
                 swap=False,\
                 swaplab=True,\
                 invert=False,\
                 xcoeff=1.,\
                 ycoeff=1.,\
                 fmt=None,\
                 colorb="jet",\
                 units="",\
                 title=""):
        ## what could be defined by the user
        self.var = var
        self.field = field
        self.absc = absc
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.title = title
        self.div = div
        self.logx = logx
        self.logy = logy
        self.swap = swap
        self.swaplab = swaplab # NB: swaplab only used if swap=True
        self.invert = invert
        self.xcoeff = xcoeff
        self.ycoeff = ycoeff
        self.fmt = fmt
        self.units = units
        self.colorb = colorb
        ## other useful arguments
        ## ... not used here in ppplot but need to be attached to plot object
        self.axisbg = "white"
        self.superpose = False

    # check
    # -------------------------------
    def check(self):
        if self.field is None: print "!! ERROR !! Please define a field to be plotted" ; exit()

    # define_from_var
    # ... this uses settings in set_var.txt
    # -------------------------------
    def define_from_var(self):
        if self.var is not None:
         if self.var.upper() in vl.keys():
          self.title = vl[self.var.upper()]
          self.units = vu[self.var.upper()]

    # make
    # this is generic to all plots
    # -------------------------------
    def make(self):
        self.check()
        # labels, title, etc...
        mpl.xlabel(self.xlabel)
        mpl.ylabel(self.ylabel)
        if self.swap:
         if self.swaplab:
           mpl.xlabel(self.ylabel)
           mpl.ylabel(self.xlabel)
        mpl.title(self.title)
        # if masked array, set masked values to filled values (e.g. np.nan) for plotting purposes
        if type(self.field).__name__ in 'MaskedArray':
            self.field[self.field.mask] = self.field.fill_value

################################
# a subclass to make a 1D plot #
################################
class plot1d(plot):

    # print out a help string when help is invoked on the object
    # -------------------------------
    def __repr__(self):
        whatprint = 'plot1d object. \"help(plot1d)\" for more information\n'
        return whatprint

    # default settings
    # -- user can define settings by two methods. 
    # -- 1. yeah = plot1d(title="foo")
    # -- 2. yeah = pp() ; yeah.title = "foo"
    # -------------------------------
    def __init__(self,\
                 lstyle=None,\
                 color=None,\
                 marker='x',\
                 label=None):
        ## get initialization from parent class
        plot.__init__(self)
        ## what could be defined by the user
        self.lstyle = lstyle
        self.color = color
        self.marker = marker
        self.label = label
        self.xmin = None
        self.ymin = None
        self.xmax = None
        self.ymax = None

    # define_from_var
    # ... this uses settings in set_var.txt
    # -------------------------------
    def define_from_var(self):
        # get what is done in the parent class
        plot.define_from_var(self)
        # add specific stuff
        if self.var is not None:
         if self.var.upper() in vl.keys():
          self.ylabel = vl[self.var.upper()] + " (" + vu[self.var.upper()] + ")"
          self.title = ""
          self.fmt = vf[self.var.upper()]

    # make
    # -------------------------------
    def make(self):
        # get what is done in the parent class
        plot.make(self)
        if self.fmt is None: self.fmt = '%.0f'
        # add specific stuff
        mpl.grid(color='grey')
        if self.lstyle == "": self.lstyle = " " # to allow for no line at all with ""
        # set dummy x axis if not defined
        if self.absc is None: 
            self.absc = np.array(range(self.field.shape[0]))
            print "!! WARNING !! dummy coordinates on x axis"
        # swapping if requested
        if self.swap:  x = self.field ; y = self.absc
        else:          x = self.absc ; y = self.field
        # coefficients on axis
        x=x*self.xcoeff ; y=y*self.ycoeff
        # check axis
        if x.size != y.size:
            print "!! ERROR !! x and y sizes don't match on 1D plot.", x.size, y.size
            exit()
        # make the 1D plot
        # either request linestyle or let matplotlib decide
        if self.lstyle is not None and self.color is not None:
            mpl.plot(x,y,self.color+self.lstyle,marker=self.marker,label=self.label)
        elif self.color is not None:
            mpl.plot(x,y,color=self.color,marker=self.marker,label=self.label)
        elif self.lstyle is not None:
            mpl.plot(x,y,linestyle=self.lstyle,marker=self.marker,label=self.label)
        else:
            mpl.plot(x,y,marker=self.marker,label=self.label)
        # make log axes and/or invert ordinate
        # ... this must be after plot so that axis bounds are well-defined
        # ... also inverting must be after making the thing logarithmic
        if self.logx: mpl.xscale("log") # not mpl.semilogx() because excludes log on y
        if self.logy: mpl.yscale("log") # not mpl.semilogy() because excludes log on x
        if self.invert: ax = mpl.gca() ; ax.set_ylim(ax.get_ylim()[::-1])
        # add a label for line(s)
        if self.label is not None:
            if self.label != "":
                mpl.legend(loc="best",fancybox=True)
        # AXES
        ax = mpl.gca()
        # format labels
        if self.swap: ax.xaxis.set_major_formatter(FormatStrFormatter(self.fmt))
        else: ax.yaxis.set_major_formatter(FormatStrFormatter(self.fmt))
        # plot limits: ensure that no curve would be outside the window
        # x-axis
        x1, x2 = ax.get_xbound()
        xmin = ppcompute.min(x)
        xmax = ppcompute.max(x)
        if xmin < x1: x1 = xmin
        if xmax > x2: x2 = xmax
        if self.xmin is not None: x1 = self.xmin
        if self.xmax is not None: x2 = self.xmax
        ax.set_xbound(lower=x1,upper=x2)
        # y-axis
        y1, y2 = ax.get_ybound()
        ymin = ppcompute.min(y)
        ymax = ppcompute.max(y)
        if ymin < y1: y1 = ymin
        if ymax > y2: y2 = ymax
        if self.ymin is not None: y1 = self.ymin
        if self.ymax is not None: y2 = self.ymax
        ax.set_ybound(lower=y1,upper=y2)
        ## set with .div the number of ticks. (is it better than automatic?)
        if not self.logx:
            ax.xaxis.set_major_locator(MaxNLocator(self.div/2))
        else:
            print "!! WARNING. in logx mode, ticks are set automatically."

################################
# a subclass to make a 2D plot #
################################
class plot2d(plot):

    # print out a help string when help is invoked on the object
    # -------------------------------
    def __repr__(self):
        whatprint = 'plot2d object. \"help(plot2d)\" for more information\n'
        return whatprint

    # default settings
    # -- user can define settings by two methods. 
    # -- 1. yeah = plot2d(title="foo")
    # -- 2. yeah = pp() ; yeah.title = "foo"
    # -------------------------------
    def __init__(self,\
                 ordi=None,\
                 mapmode=False,\
                 proj="cyl",\
                 back=None,\
                 trans=1.0,\
                 addvecx=None,\
                 addvecy=None,\
                 addcontour=None,\
                 blon=None,\
                 blat=None,\
                 area=None,\
                 vmin=None,\
                 vmax=None,\
                 showcb=True,\
                 wscale=None,\
                 stridevecx=1,\
                 stridevecy=1,\
                 colorvec="black"):
        ## get initialization from parent class
        plot.__init__(self)
        ## what could be defined by the user
        self.ordi = ordi
        self.mapmode = mapmode
        self.proj = proj
        self.back = back
        self.trans = trans
        self.addvecx = addvecx
        self.addvecy = addvecy
        self.colorvec = colorvec
        self.addcontour = addcontour
        self.blon = blon ; self.blat = blat
        self.area = area
        self.vmin = vmin ; self.vmax = vmax
        self.showcb = showcb
        self.wscale = wscale
        self.stridevecx = stridevecx
        self.stridevecy = stridevecy

    # define_from_var
    # ... this uses settings in set_var.txt
    # -------------------------------
    def define_from_var(self):
        # get what is done in the parent class
        plot.define_from_var(self)
        # add specific stuff
        if self.var is not None:
         if self.var.upper() in vl.keys():
          self.colorb = vc[self.var.upper()]
          self.fmt = vf[self.var.upper()]

    # make
    # -------------------------------
    def make(self):
        # get what is done in the parent class...
        plot.make(self)
        if self.fmt is None: self.fmt = "%.1e"
        # ... then add specific stuff
        ############################################################################################
        ### PRE-SETTINGS
        ############################################################################################
        # set dummy xy axis if not defined
        if self.absc is None: 
            self.absc = np.array(range(self.field.shape[0]))
            self.mapmode = False
            print "!! WARNING !! dummy coordinates on x axis"
        if self.ordi is None: 
            self.ordi = np.array(range(self.field.shape[1]))
            self.mapmode = False
            print "!! WARNING !! dummy coordinates on y axis"
        # transposing if necessary
        shape = self.field.shape
        if shape[0] != shape[1]:
         if len(self.absc) == shape[0] and len(self.ordi) == shape[1]:
            print "!! WARNING !! Transposing axes"
            self.field = np.transpose(self.field)
        # bound field
        zevmin, zevmax = calculate_bounds(self.field,vmin=self.vmin,vmax=self.vmax)
        what_I_plot = bounds(self.field,zevmin,zevmax)
        # define contour field levels. define color palette
        ticks = self.div + 1
        zelevels = np.linspace(zevmin,zevmax,ticks)
        palette = get_cmap(name=self.colorb)
        # do the same thing for possible contourline entries
        if self.addcontour is not None:
            # if masked array, set masked values to filled values (e.g. np.nan) for plotting purposes
            if type(self.addcontour).__name__ in 'MaskedArray':
               self.addcontour[self.addcontour.mask] = self.addcontour.fill_value
            zevminc, zevmaxc = calculate_bounds(self.addcontour)
            what_I_contour = bounds(self.addcontour,zevminc,zevmaxc)
            ticks = self.div + 1
            zelevelsc = np.linspace(zevminc,zevmaxc,ticks)
        ############################################################################################
        ### MAIN PLOT
        ### NB: contour lines are done before contour shades otherwise colorar error
        ############################################################################################
        if not self.mapmode:
            ## A SIMPLE 2D PLOT
            ###################
            # swapping if requested
            if self.swap:  x = self.ordi ; y = self.absc
            else:          x = self.absc ; y = self.ordi
            # coefficients on axis
            x=x*self.xcoeff ; y=y*self.ycoeff
            # make shaded and line contours
            if self.addcontour is not None: 
                objC = mpl.contour(x, y, what_I_contour, \
                            zelevelsc, colors = ccol, linewidths = cline)
                #mpl.clabel(objC, inline=1, fontsize=10)
            mpl.contourf(x, y, \
                         self.field, \
                         zelevels, cmap=palette)
            # make log axes and/or invert ordinate
            if self.logx: mpl.semilogx()
            if self.logy: mpl.semilogy()
            if self.invert: ax = mpl.gca() ; ax.set_ylim(ax.get_ylim()[::-1])
        else:
            ## A 2D MAP USING PROJECTIONS (basemap)
            #######################################
            mpl.xlabel("") ; mpl.ylabel("")
            # additional security in case self.proj is None here
            # ... we set cylindrical projection (the simplest one)
            if self.proj is None: self.proj = "cyl"
            # get lon and lat in 2D version.
            # (but first ensure we do have 2D coordinates)
            if self.absc.ndim == 1: 	[self.absc,self.ordi] = np.meshgrid(self.absc,self.ordi)
            elif self.absc.ndim > 2:    print "!! ERROR !! lon and lat arrays must be 1D or 2D"
            # get lat lon intervals and associated settings
            wlon = [np.min(self.absc),np.max(self.absc)]
            wlat = [np.min(self.ordi),np.max(self.ordi)]
            # -- area presets are in set_area.txt
            if self.area is not None:
             if self.area in area.keys():
                wlon, wlat = area[self.area]
            # -- settings for meridians and parallels
            steplon = int(abs(wlon[1]-wlon[0])/6.)
            steplat = int(abs(wlat[1]-wlat[0])/3.)
            #mertab = np.r_[wlon[0]:wlon[1]:steplon] ; merlab = [0,0,0,1]
            #partab = np.r_[wlat[0]:wlat[1]:steplat] ; parlab = [1,0,0,0]
            mertab = np.r_[-180.:180.:steplon] ; merlab = [0,0,0,1] #-360:360.
            partab = np.r_[-90.:90.:steplat] ; parlab = [1,0,0,0]
            format = '%.1f'
            # -- center of domain and bounding lats
            lon_0 = 0.5*(wlon[0]+wlon[1])
            lat_0 = 0.5*(wlat[0]+wlat[1])
            # some tests, bug fixes, and good-looking settings
            # ... cyl is good for global and regional
            if self.proj == "cyl":
                format = '%.0f'
            # ... global projections
            elif self.proj in ["ortho","moll","robin"]:
                wlat[0] = None ; wlat[1] = None ; wlon[0] = None ; wlon[1] = None
                steplon = 30. ; steplat = 30.
                if self.proj in ["robin","moll"]: steplon = 60.
                mertab = np.r_[-360.:360.:steplon]
                partab = np.r_[-90.:90.:steplat]
                if self.proj == "ortho": 
                    merlab = [0,0,0,0] ; parlab = [0,0,0,0]
                    # in ortho projection, blon and blat can be used to set map center
                    if self.blon is not None: lon_0 = self.blon
                    if self.blat is not None: lat_0 = self.blat
                elif self.proj == "moll":
                    merlab = [0,0,0,0]
                format = '%.0f'
            # ... regional projections
            elif self.proj in ["lcc","laea","merc"]:
                if self.proj in ["lcc","laea"] and wlat[0] == -wlat[1]: 
                    print "!! ERROR !! with Lambert lat1 must be different than lat2" ; exit()
                if wlat[0] < -80. and wlat[1] > 80.:
                    print "!! ERROR !! set an area (not global)" ; exit()
                format = '%.0f'
            elif self.proj in ["npstere","spstere"]:
                # in polar projections, blat gives the bounding lat
                # if not set, set something reasonable
                if self.blat is None:   self.blat = 60.
                # help the user who forgets self.blat would better be negative in spstere
                # (this actually serves for the default setting just above)
                if self.proj == "spstere" and self.blat > 0: self.blat = -self.blat
                # labels
                mertab = np.r_[-360.:360.:15.]
                partab = np.r_[-90.:90.:5.]
            # ... unsupported projections
            else:
                print "!! ERROR !! unsupported projection. supported: "+\
                      "cyl, npstere, spstere, ortho, moll, robin, lcc, laea, merc"
            # finally define projection
            m = Basemap(projection=self.proj,\
                        lat_0=lat_0,lon_0=lon_0,\
                        boundinglat=self.blat,\
                        llcrnrlat=wlat[0],urcrnrlat=wlat[1],\
                        llcrnrlon=wlon[0],urcrnrlon=wlon[1])
            # draw meridians and parallels
            ft = int(mpl.rcParams['font.size']*3./4.)
            zelatmax = 85.
            m.drawmeridians(mertab,labels=merlab,color='grey',linewidth=0.75,fontsize=ft,fmt=format,latmax=zelatmax)
            m.drawparallels(partab,labels=parlab,color='grey',linewidth=0.75,fontsize=ft,fmt=format,latmax=zelatmax)
            # define background (see set_back.txt)
            if self.back is not None:
              if self.back in back.keys():
                 print "**** info: loading a background, please wait.",self.back
                 if self.back not in ["coast","sea"]:   m.warpimage(back[self.back],scale=0.75)
                 elif self.back == "coast":             m.drawcoastlines()
                 elif self.back == "sea":               m.drawlsmask(land_color='white',ocean_color='aqua')
              else:
                 print "!! ERROR !! requested background not defined. change name or fill in set_back.txt" ; exit()
            # define x and y given the projection
            x, y = m(self.absc, self.ordi)
            # contour field. first line contour then shaded contour.
            if self.addcontour is not None: 
                objC2 = m.contour(x, y, what_I_contour, \
                            zelevelsc, colors = ccol, linewidths = cline)
                #mpl.clabel(objC2, inline=1, fontsize=10)
            m.contourf(x, y, what_I_plot, zelevels, cmap = palette, alpha = self.trans)
        ############################################################################################
        ### COLORBAR
        ############################################################################################
        if self.trans > 0. and self.showcb:
            ## draw colorbar. settings are different with projections. or if not mapmode.
            if not self.mapmode: orientation=zeorientation ; frac = 0.075 ; pad = 0.03 ; lu = 0.5
            elif self.proj in ['moll']: orientation="horizontal" ; frac = 0.08 ; pad = 0.03 ; lu = 1.0
            elif self.proj in ['cyl']: orientation="vertical" ; frac = 0.023 ; pad = 0.03 ; lu = 0.5
            else: orientation = zeorientation ; frac = zefrac ; pad = 0.03 ; lu = 0.5
            zelevpal = np.linspace(zevmin,zevmax,num=min([ticks/2+1,21]))
            zecb = mpl.colorbar(fraction=frac,pad=pad,\
                                format=self.fmt,orientation=orientation,\
                                ticks=zelevpal,\
                                extend='neither',spacing='proportional')
            if zeorientation == "horizontal": zecb.ax.set_xlabel(self.title) ; self.title = ""
            # colorbar title --> units
            if self.units not in ["dimless",""]:
                zecb.ax.set_title("["+self.units+"]",fontsize=3.*mpl.rcParams['font.size']/4.,x=lu,y=1.025)

        ############################################################################################
        ### VECTORS. must be after the colorbar. we could also leave possibility for streamlines.
        ############################################################################################
        ### not expecting NaN in self.addvecx and self.addvecy. masked arrays is just enough.
        if self.addvecx is not None and self.addvecy is not None: 
                # vectors on map projection or simple 2D mapping
                if self.mapmode: [vecx,vecy] = m.rotate_vector(self.addvecx,self.addvecy,self.absc,self.ordi) # for metwinds only ?
                else: vecx,vecy = self.addvecx,self.addvecy ; x,y = np.meshgrid(x,y)
                # reference vector is scaled
                if self.wscale is None: self.wscale = ppcompute.mean(np.sqrt(self.addvecx*self.addvecx+self.addvecy*self.addvecy))
                # make vector field
                if self.mapmode: 
                    q = m.quiver( x[::self.stridevecy,::self.stridevecx],y[::self.stridevecy,::self.stridevecx],\
                                  vecx[::self.stridevecy,::self.stridevecx],vecy[::self.stridevecy,::self.stridevecx],\
                                  angles='xy',color=self.colorvec,pivot='tail',\
                                  scale=self.wscale*reducevec,width=widthvec )
                else:
                    q = mpl.quiver( x[::self.stridevecy,::self.stridevecx],y[::self.stridevecy,::self.stridevecx],\
                                    vecx[::self.stridevecy,::self.stridevecx],vecy[::self.stridevecy,::self.stridevecx],\
                                    angles='xy',color=self.colorvec,pivot='tail',\
                                    scale=self.wscale*reducevec,width=widthvec )
                # make vector key.
                #keyh = 1.025 ; keyv = 1.05 # upper right corner over colorbar
                keyh = 0.97 ; keyv = 1.05
                #keyh = -0.03 ; keyv = 1.08 # upper left corner
                p = mpl.quiverkey(q,keyh,keyv,\
                                  self.wscale,str(int(self.wscale)),\
                                  fontproperties={'size': 'small'},\
                                  color='black',labelpos='S',labelsep = 0.07)
        ############################################################################################
        ### TEXT. ANYWHERE. add_text.txt should be present with lines x ; y ; text ; color
        ############################################################################################
        try:
            f = open("add_text.txt", 'r')
            for line in f:
              if "#" in line: pass
              else:
                  userx, usery, usert, userc = line.strip().split(';')
                  userc = userc.strip()
                  usert = usert.strip()
                  userx = float(userx.strip())
                  usery = float(usery.strip())
                  if self.mapmode: userx,usery = m(userx,usery)
                  mpl.text(userx,usery,usert,\
                           color = userc,\
                           horizontalalignment='center',\
                           verticalalignment='center')
            f.close()
        except IOError:
            pass

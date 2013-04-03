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
mpl.rcParams['axes.color_cycle'] = "r,g,b,k"
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
vf = {} ; vc = {} ; vl = {}
try: 
    f = open(whereset+zefile, 'r')
    for line in f:
        if "#" in line: pass
        else:
            var, format, colorb, label = line.strip().split(';')
            ind = var.strip() 
            vf[ind] = format.strip()
            vc[ind] = colorb.strip()
            vl[ind] = label.strip()
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
       if zevmin < 0. and ppcompute.min(fieldcalc) > 0.: zevmin = 0.
    # check that bounds are not too tight given the field
    amin = ppcompute.min(field)
    amax = ppcompute.max(field)
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
                 lstyle='-',\
                 color='b',\
                 marker='x',\
                 label=None):
        ## get initialization from parent class
        plot.__init__(self)
        ## what could be defined by the user
        self.lstyle = lstyle
        self.color = color
        self.marker = marker
        self.label = label

    # define_from_var
    # ... this uses settings in set_var.txt
    # -------------------------------
    def define_from_var(self):
        # get what is done in the parent class
        plot.define_from_var(self)
        # add specific stuff
        if self.var is not None:
         if self.var.upper() in vl.keys():
          self.ylabel = vl[self.var.upper()]
          self.title = ""

    # make
    # -------------------------------
    def make(self):
        # get what is done in the parent class
        plot.make(self)
        # add specific stuff
        mpl.grid(color='grey')
        if self.lstyle == "": self.lstyle = " " # to allow for no line at all with ""
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
        else:
            mpl.plot(x,y,marker=self.marker,label=self.label)
        # make log axes and/or invert ordinate
        # ... this must be after plot so that axis bounds are well-defined
        # ... also inverting must be after making the thing logarithmic
        if self.logx: mpl.semilogx()
        if self.logy: mpl.semilogy()
        if self.invert: ax = mpl.gca() ; ax.set_ylim(ax.get_ylim()[::-1])
        # add a label for line(s)
        if self.label is not None:
            if self.label != "":
                mpl.legend(loc="best",fancybox=True)
        ## TBD: set with .div the number of ticks
        ## TBD: be able to control plot limits
        #ticks = self.div + 1
        #ax = mpl.gca()
        #ax.get_xaxis().set_ticks(np.linspace(ppcompute.min(x),ppcompute.max(x),ticks/2+1))

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
                 colorb="jet",\
                 trans=1.0,\
                 addvecx=None,\
                 addvecy=None,\
                 addcontour=None,\
                 fmt="%.2e",\
                 blon=None,\
                 blat=None,\
                 area=None,\
                 vmin=None,\
                 vmax=None,\
                 colorvec="black"):
        ## get initialization from parent class
        plot.__init__(self)
        ## what could be defined by the user
        self.ordi = ordi
        self.mapmode = mapmode
        self.proj = proj
        self.back = back
        self.colorb = colorb
        self.trans = trans
        self.addvecx = addvecx
        self.addvecy = addvecy
        self.colorvec = colorvec
        self.addcontour = addcontour
        self.fmt = fmt
        self.blon = blon ; self.blat = blat
        self.area = area
        self.vmin = vmin ; self.vmax = vmax

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
        # ... then add specific stuff
        ############################################################################################
        ### PRE-SETTINGS
        ############################################################################################
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
                mpl.contour(x, y, what_I_contour, \
                            zelevelsc, colors = ccol, linewidths = cline)
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
            steplon = abs(wlon[1]-wlon[0])/6.
            steplat = abs(wlat[1]-wlat[0])/3.
            mertab = np.r_[wlon[0]:wlon[1]:steplon] ; merlab = [0,0,0,1]
            partab = np.r_[wlat[0]:wlat[1]:steplat] ; parlab = [1,0,0,0]
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
            m.drawmeridians(mertab,labels=merlab,color='grey',linewidth=0.75,fontsize=ft,fmt=format)
            m.drawparallels(partab,labels=parlab,color='grey',linewidth=0.75,fontsize=ft,fmt=format)
            # define background (see set_back.txt)
            if self.back is not None:
              if self.back in back.keys():
                 print "**** info: loading a background, please wait.",self.back
                 if back not in ["coast","sea"]:   m.warpimage(back[self.back],scale=0.75)
                 elif back == "coast":             m.drawcoastlines()
                 elif back == "sea":               m.drawlsmask(land_color='white',ocean_color='aqua')
              else:
                 print "!! ERROR !! requested background not defined. change name or fill in set_back.txt" ; exit()
            # define x and y given the projection
            x, y = m(self.absc, self.ordi)
            # contour field. first line contour then shaded contour.
            if self.addcontour is not None: 
                m.contour(x, y, what_I_contour, \
                            zelevelsc, colors = ccol, linewidths = cline)
            m.contourf(x, y, what_I_plot, zelevels, cmap = palette, alpha = self.trans)
        ############################################################################################
        ### COLORBAR
        ############################################################################################
        if self.trans > 0.:
            ## draw colorbar. settings are different with projections. or if not mapmode.
            if not self.mapmode: orientation=zeorientation ; frac = 0.075 ; pad = 0.03
            elif self.proj in ['moll']: orientation="horizontal" ; frac = 0.075 ; pad = 0.03
            elif self.proj in ['cyl']: orientation="vertical" ; frac = 0.023 ; pad = 0.03
            else: orientation = zeorientation ; frac = zefrac ; pad = 0.03
            zelevpal = np.linspace(zevmin,zevmax,num=min([ticks/2+1,21]))
            zecb = mpl.colorbar(fraction=frac,pad=pad,\
                                format=self.fmt,orientation=orientation,\
                                ticks=zelevpal,\
                                extend='neither',spacing='proportional') 
            if zeorientation == "horizontal": zecb.ax.set_xlabel(self.title) ; self.title = ""
        ############################################################################################
        ### VECTORS. must be after the colorbar. we could also leave possibility for streamlines.
        ############################################################################################
        ### not expecting NaN in self.addvecx and self.addvecy. masked arrays is just enough.
        stride = 3
        if self.addvecx is not None and self.addvecy is not None and self.mapmode:
                ## for metwinds only ???
                [vecx,vecy] = m.rotate_vector(self.addvecx,self.addvecy,self.absc,self.ordi) 
                # reference vector is scaled
                zescale = ppcompute.mean(np.sqrt(self.addvecx*self.addvecx+self.addvecy*self.addvecy))
                # make vector field
                q = m.quiver( x[::stride,::stride],y[::stride,::stride],\
                              vecx[::stride,::stride],vecy[::stride,::stride],\
                              angles='xy',color=self.colorvec,pivot='middle',\
                              scale=zescale*reducevec,width=widthvec )
                # make vector key. default is on upper left corner.
                keyh = 1.025 ; keyv = 1.05
                #keyh = -0.03 ; keyv = 1.08
                p = mpl.quiverkey(q,keyh,keyv,\
                                  zescale,str(int(zescale)),\
                                  color='black',labelpos='S',labelsep = 0.07)

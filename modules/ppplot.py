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
#   (None means planetoplot in PYTHONPATH)
whereset = None
# - some good default settings.
# (contours)
cline = 0.55
#cline = 0.8
# (vectors)
widthvec = 0.002
reducevec = 30.
# (colorbar)
zeorientation="vertical"
zefrac = 0.05
# (save figures)
pad_inches_value=0.25
###############################################



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
            ind = var.strip().upper()
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

# a function to quickly obtain a plot
def quickplot(field):
    ddd = field.ndim
    if ddd == 1: plot1d(f=field).makeshow()
    elif ddd == 2: plot2d(f=field).makeshow()

# continuity with matplotlib
def close():
    mpl.close()

# a function for predefined figure sizes
def figuref(x=16,y=9):
    fig = mpl.figure(figsize=(x,y))
    return fig

# a function to change color lines according to a color map (idea by A. Pottier)
# ... two optional arguments: color maps + a number telling how much intervals are needed
def rainbow(cb="jet",num=8):
    ax = mpl.gca()
    pal = mpl.cm.get_cmap(name=cb)
    ax._get_lines.set_color_cycle([pal(i) for i in np.linspace(0,0.9,num)])

# a function to define subplot
# ... user can change settings in set_multiplot.txt read above
# -------------------------------
def definesubplot(numplot, fig, factor=1., sup=False):
    # default
    fontsize = 18 ; subv = 1 ; subh = 1
    # multiplot mode
    if not sup:
     try:
       fontsize = font_t[numplot]
       fig.subplots_adjust(wspace = wspace_t[numplot], hspace = hspace_t[numplot])
       subv, subh = subv_t[numplot],subh_t[numplot]
     except:
       print "!! WARNING !! no preset found from set_multiplot.txt, or this setting file was not found."
       subv = 1 ; subh = numplot
    # commensurate font
    if factor != 0.: fontsize = fontsize / factor
    # change parameter
    mpl.rcParams['font.size'] = fontsize
    return subv,subh

# a function to calculate automatically bounds (or simply prescribe those)
# -------------------------------
def calculate_bounds(field,vmin=None,vmax=None,sigma=None):
    # prescribed cases first
    zevmin = vmin
    zevmax = vmax
    # particular case: only nan in field
    if False not in np.isnan(field):
      zevmin = np.nan
      zevmax = np.nan
    # GENERAL cases to be computed
    else:
     if zevmin is None or zevmax is None:
       # calculate min and max
       ind = np.where(np.abs(field) < 9e+35) # select values
       fieldcalc = field[ ind ] # field must be a numpy array
       amin = ppcompute.min(field)
       amax = ppcompute.max(field)
       # default case: sigma is None, take min and max
       if sigma is None:
          zevmin = amin
          zevmax = amax
       # if sigma not None: calculate stdev and mean
       else:
          dev = np.std(fieldcalc)*sigma
          damean = ppcompute.mean(fieldcalc)
          # fill min/max if needed
          if vmin is None: zevmin = damean - dev
          if vmax is None: zevmax = damean + dev
          # special case: negative values with stddev while field is positive
          if zevmin < 0. and ppcompute.min(fieldcalc) >= 0.: zevmin = 0.
          # check that bounds are not too tight given the field
          if np.abs(amin) < 1.e-15: cmin = 0.
          else: cmin = 100.*np.abs((amin - zevmin)/amin)
          cmax = 100.*np.abs((amax - zevmax)/amax)
          if cmin > 150. or cmax > 150.:
            print "!! WARNING !! Bounds are a bit too tight. Might need to reconsider those."
            print "!! WARNING !! --> actual",amin,amax,"adopted",zevmin,zevmax
    return zevmin, zevmax    

# a function to solve the problem with blank bounds !
# -------------------------------
def bounds(what_I_plot,zevmin,zevmax,miss=9e+35):
    cond = (not np.isnan(zevmin)) or (not np.isnan(zevmax))
    if cond:
      small_enough = 1.e-7
      if zevmin < 0: what_I_plot[ what_I_plot < zevmin*(1.-small_enough) ] = zevmin*(1.-small_enough)
      else:          what_I_plot[ what_I_plot < zevmin*(1.+small_enough) ] = zevmin*(1.+small_enough)
      what_I_plot[ what_I_plot > miss  ] = -miss
      what_I_plot[ what_I_plot > zevmax ] = zevmax*(1.-small_enough)
    return what_I_plot

# a function to change labels with modulo
# ---------------------------------------
def labelmodulo(ax,mod):
    mpl.draw()
    strtab = []
    for tick in ax.get_xaxis().get_ticklabels():
        onetick = tick.get_text()
        if len(onetick) > 0:
          num = float(onetick)
          strtab.append(str(num % mod))
        elif len(onetick)==0:
          strtab.append('')
    ax.get_xaxis().set_ticklabels(strtab)
    return ax

# a function to output an ascii file
# ----------------------------------
def writeascii (field=None,absc=None,name=None):
    field = np.array(field)
    absc = np.array(absc)
    if name is None:
        name = "prof"
        for ttt in timelib.gmtime():
            name = name + "_" + str(ttt)
    if field is None or absc is None:
        print "!! WARNING !! Not printing the file, incorrect field or absc."
    else:
        if field.ndim == 1:
            myfile = open(name, 'w')
            for ix in range(len(absc)):
                myfile.write("%15.5e%15.5e\n" % ( absc[ix], field[ix] ))
            myfile.close()
        else:
            print "!! WARNING !! Not printing the file, 2D fields not supported yet."
    return

# a necessary addition for people used to matplotlib
def show():
    mpl.show()

# a generic function to show (GUI) or save a plot (PNG,EPS,PDF,...)
# -------------------------------
def save(mode=None,filename=None,folder="./",includedate=False,res=150,custom=False):
    # no filename or no mode set --> graphical user interface
    if filename is None or mode is None:
      show()
    # otherwise --> an image is saved
    else:
      # a few settings
      possiblesave = ['eps','ps','svg','png','jpg','pdf'] 
      # now the main instructions
      if mode in possiblesave:
          ## name of plot
          name = folder+'/'+filename
          if includedate:
              for ttt in timelib.gmtime():
                  name = name + "_" + "%03d" % (ttt)
          name = name +"."+mode
          ## save file
          print "**** Saving file in "+mode+" format... Please wait."
          if not custom:
              # ... regular plots
              mpl.savefig(name,dpi=res,pad_inches=pad_inches_value,bbox_inches='tight')
          else:
              # ... mapping mode, adapted space for labels etc...
              # TBD: not impacted by pad_inches_value. why?
              mpl.savefig(name,dpi=res) #,bbox_inches='tight')
      else:
          print "!! ERROR !! File format not supported. Supported: ",possiblesave

## make it negative (black background / white axes)
def negative(howblack="black"):
    mpl.rc('figure', facecolor=howblack, edgecolor=howblack)
    mpl.rcParams['text.color'] = 'w'
    mpl.rcParams['lines.color'] = 'w'
    mpl.rc('axes',labelcolor='w',facecolor=howblack,edgecolor='w')
    mpl.rcParams['xtick.color'] = 'w'
    mpl.rcParams['ytick.color'] = 'w'
    mpl.rcParams['grid.color'] = 'w'
    mpl.rc('savefig',facecolor=howblack,edgecolor=howblack)

### settings for xkcd mode (only with matplotlib 1.3)
### ... you must have Humori-Sans Font installed
def xkcd():
    mpl.xkcd()

## for command-line use: set up parser options for plot, plot1D and plot2D
## -- must take as an input a parser object
## -- return a modified parser object
def opt(parser):
  parser.add_option('-T','--title',action='append',dest='title',type="string",default=None,help="change 'title'")
  parser.add_option('-X','--xlabel',action='append',dest='xlabel',type="string",default=None,help="change 'xlabel'")
  parser.add_option('-Y','--ylabel',action='append',dest='ylabel',type="string",default=None,help="change 'ylabel'")
  parser.add_option('-D','--div',action='store',dest='div',type="int",default=20,help="integer for number of divisions")
  parser.add_option('-H','--trans',action='store',dest='trans',type="float",default=1.0,help="float for transparency (0=transp,1=opaque)")
  parser.add_option('-Z','--logy',action='store_true',dest='logy',default=False,help="set log for vertical axis")
  parser.add_option('--logx',action='store_true',dest='logx',default=False,help="set log for horizontal axis")
  parser.add_option('-O','--save',action='store',dest='out',type="string",default=None,help="save mode: 'gui' 'png' 'pdf' 'eps' 'svg' 'ps'")
  parser.add_option('-V','--void',action='store_true',dest='void',default=False,help="no colorbar, no title, no labels")
  parser.add_option('-U','--units',action='append',dest='units',type="string",default=None,help="units for the field")
  parser.add_option('-F','--fmt',action='append',dest='fmt',type="string",default=None,help="values formatting. ex: '%.0f' '%3.1e'")
  parser.add_option('--xcoeff',action='append',dest='xcoeff',type="float",default=None,help="multiply x axis [not for 2D map]")
  parser.add_option('--ycoeff',action='append',dest='ycoeff',type="float",default=None,help="multiply y axis [not for 2D map]")
  parser.add_option('--xmin',action='append',dest='xmin',type="float",default=None,help="min bound x axis [not for 2D map]")
  parser.add_option('--ymin',action='append',dest='ymin',type="float",default=None,help="min bound y axis [not for 2D map]")
  parser.add_option('--xmax',action='append',dest='xmax',type="float",default=None,help="max bound x axis [not for 2D map]")
  parser.add_option('--ymax',action='append',dest='ymax',type="float",default=None,help="max bound y axis [not for 2D map]")
  parser.add_option('--nxticks',action='append',dest='nxticks',type="float",default=None,help="ticks for x axis [not for 2D map]")
  parser.add_option('--nyticks',action='append',dest='nyticks',type="float",default=None,help="ticks for y axis [not for 2D map]")
  parser.add_option('--xp',action='store',dest='xp',type="int",default=None,help="x size of figure (integer)")
  parser.add_option('--yp',action='store',dest='yp',type="int",default=None,help="y size of figure (integer)")
  parser.add_option('--res',action='store',dest='res',type="float",default=150.,help="change resolution of figure")
  parser.add_option('--swap',action='store_true',dest='swap',default=False,help='swap axis')
  return parser
def opt1d(parser):
  parser.add_option('-L','--linestyle',action='append',dest='linestyle',type="string",default=None,help="[1D] linestyle: '-' '--' '.' '..'")
  parser.add_option('-Q','--color',action='append',dest='color',type="string",default=None,help="[1D] color: 'b' 'g' 'r' etc")
  parser.add_option('-K','--marker',action='append',dest='marker',type="string",default=None,help="[1D] marker: None '' 'x' 'o' etc")
  parser.add_option('-S','--superpose',action='store_true',dest='superpose',default=False,help="[1D] use same axis for all plots")
  parser.add_option('-E','--legend',action='append',dest='legend',type="string",default=None,help="[1D] legend for line")
  parser.add_option('--modx',action='append',dest='modx',type="float",default=None,help="[1D] change xticks with a modulo")
  return parser
def opt2d(parser):
  parser.add_option('-C','--colorbar',action='append',dest='colorbar',type="string",default=None,help="[2D] colormap: http://micropore.files.wordpress.com/2010/06/colormaps.png")
  parser.add_option('-P','--proj',action='store',dest='proj',type="string",default=None,help="[2D] map projection: 'cyl' 'npstere' 'spstere' 'ortho' 'moll' 'robin' 'lcc' 'laea' 'merc'")
  parser.add_option('-B','--back',action='append',dest='back',type="string",default=None,help='[2D] predefined map background (cf. set_back.txt)')
  parser.add_option('-A','--area',action='append',dest='area',type="string",default=None,help='[2D] predefined region of mapping (cf. set_area.txt)')
  parser.add_option('-I','--blon',action='append',dest='blon',type="float",default=None,help='[2D] float: bounding longitude for stere (or center longitude for ortho)')
  parser.add_option('-J','--blat',action='append',dest='blat',type="float",default=None,help='[2D] float: bounding latitude for stere (or center latitude for ortho) ')
  parser.add_option('-N','--vmin',action='append',dest='vmin',type="float",default=None,help='[2D] float: minimum value for displayed field')
  parser.add_option('-M','--vmax',action='append',dest='vmax',type="float",default=None,help='[2D] float: maximum value for displayed field')
  parser.add_option('-W','--wscale',action='append',dest='wscale',type="float",default=None,help='[2D] float: set size of reference wind vector')
  parser.add_option('--svx',action='store',dest='svx',type="int",default=None,help="Define an abscissa stride on vectors only -- not on field")
  parser.add_option('--svy',action='store',dest='svy',type="int",default=None,help="Define an ordinate stride on vectors only -- not on field")
  parser.add_option('--cbticks',action='append',dest='cbticks',type="float",default=None,help="ticks for colorbar")
  parser.add_option('--sigma',action='store',dest='sigma',type="float",default=None,help="[2D] enhance contrast between max/min (default: min/max limits)")
  parser.add_option('--nocb',action='store_false',dest='showcb',default=True,help="[2D] do not show colorbar")
  return parser

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
                 f=None,\
                 x=None,\
                 xlabel="",\
                 ylabel="",\
                 div=20,\
                 logx=False,\
                 logy=False,\
                 swap=False,\
                 swaplab=True,\
                 invert=False,\
                 xcoeff=None,\
                 ycoeff=None,\
                 fmt=None,\
                 colorbar="jet",\
                 units="",\
                 modx=None,\
                 xmin=None,\
                 ymin=None,\
                 xmax=None,\
                 ymax=None,\
                 nxticks=10,\
                 nyticks=10,\
                 cbticks=None,\
                 xdate=False,\
                 title=""):
        ## what could be defined by the user
        self.var = var
        self.f = f
        self.x = x
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
        self.colorbar = colorbar
        self.modx = modx
        self.xmin = xmin
        self.ymin = ymin 
        self.xmax = xmax
        self.ymax = ymax
        self.nxticks = nxticks
        self.nyticks = nyticks
        self.cbticks = cbticks
        self.xdate = xdate
        ## other useful arguments
        ## ... not used here in ppplot but need to be attached to plot object
        self.axisbg = "white"
        self.superpose = False

    def transopt(self,objfrom,num=None):
      # method to transfer not-None attributes
      # from a destination object to a plot object
      # e.g. useful in command line scripts with options 
      # ------------------------------------------------
      # for all common attributes ...
      commonattr = [(k,v) for k,v in vars(objfrom).items() for kt,vt in vars(self).items() if k == kt]
      for k,v in commonattr:
        # ... if value is not None ...
        if v is not None:
          # proceed
          # 1. detect if array attributes (obtained with append)
          try: l = len(v)
          except: l = 0
          # 2. intermediate transfer
          if l>0:
            if num is not None: 
              if num >= l: 
                vv = v[l-1]
              else: 
                vv = v[num]
            else: 
              vv = v
          else:
            vv = v
          # 3. set attibutes in destination
          setattr(self,k,vv)

    # check
    # -------------------------------
    def check(self):
        if self.f is None: print "!! ERROR !! Please define a field to be plotted" ; exit()
        else: self.f = np.array(self.f) # ensure this is a numpy array

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
        mpl.title(self.title,y=1.01) # raise a little bit for subscript
        # if masked array, set masked values to filled values (e.g. np.nan) for plotting purposes
        if type(self.f).__name__ in 'MaskedArray':
            self.f[self.f.mask] = self.f.fill_value

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
                 linestyle=None,\
                 color=None,\
                 marker='x',\
                 legend=None,\
                 *args, **kwargs):
        ## get initialization from parent class
        plot.__init__(self,*args, **kwargs)
        ## what could be defined by the user
        self.linestyle = linestyle
        self.color = color
        self.marker = marker
        self.legend = legend

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
        if self.linestyle == "": self.linestyle = " " # to allow for no line at all with ""
        # set dummy x axis if not defined
        if self.x is None: 
            self.x = np.array(range(self.f.shape[0]))
            print "!! WARNING !! dummy coordinates on x axis"
        else:
            self.x = np.array(self.x) # ensure this is a numpy array
        # swapping if requested
        if self.swap:  x = self.f ; y = self.x
        else:          x = self.x ; y = self.f
        # coefficients on axis
        if self.xcoeff is not None: x=x*self.xcoeff
        if self.ycoeff is not None: y=y*self.ycoeff
        # check axis
        if x.size != y.size:
            print "!! ERROR !! x and y sizes don't match on 1D plot.", x.size, y.size
            exit()
        # make the 1D plot
        # either request linestyle or let matplotlib decide
        if self.linestyle is not None and self.color is not None:
            mpl.plot(x,y,self.color+self.linestyle,marker=self.marker,label=self.legend)
        elif self.color is not None:
            mpl.plot(x,y,color=self.color,marker=self.marker,label=self.legend)
        elif self.linestyle is not None:
            mpl.plot(x,y,linestyle=self.linestyle,marker=self.marker,label=self.legend)
        else:
            mpl.plot(x,y,marker=self.marker,label=self.legend)
        # AXES
        ax = mpl.gca()
        # make log axes and/or invert ordinate
        # ... this must be after plot so that axis bounds are well-defined
        # ... also inverting must be after making the thing logarithmic
        if self.logx: mpl.xscale("log") # not mpl.semilogx() because excludes log on y
        if self.logy: mpl.yscale("log") # not mpl.semilogy() because excludes log on x
        if self.invert: ax.set_ylim(ax.get_ylim()[::-1])
        if self.xmin is not None and self.xmax is not None:
          if self.xmin > self.xmax:
            ax.set_xlim(ax.get_xlim()[::-1])
            self.xmin,self.xmax = self.xmax,self.xmin
        # add a label for line(s)
        if self.legend is not None:
            if self.legend != "":
                mpl.legend(loc="best",fancybox=True)
        # format labels
        if self.swap: 
          if not self.logx:
            ax.xaxis.set_major_formatter(FormatStrFormatter(self.fmt))
        else: 
          if not self.logy: 
            ax.yaxis.set_major_formatter(FormatStrFormatter(self.fmt))
        # plot limits: ensure that no curve would be outside the window
        # x-axis
        if self.xdate:
          import matplotlib.dates as mdates
          ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%y'))
          #ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%b-%d %H:%M:%S'))
          ax.xaxis.set_major_locator(mdates.DayLocator())
          mpl.setp(mpl.xticks()[1], rotation=30, ha='right') # rotate the x labels
        else:
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
        ## set with .div the number of ticks.
        if not self.logx:
            ax.xaxis.set_major_locator(MaxNLocator(self.nxticks))
        else:
            # in logx mode, ticks are set automatically
            pass
        if not self.logy:
            ax.yaxis.set_major_locator(MaxNLocator(self.nyticks))
        else:
            # in logy mode, ticks are set automatically
            pass
        ## specific modulo labels
        if self.modx is not None:
            ax = labelmodulo(ax,self.modx)

    # makeshow = make + show
    # -------------------------------
    def makeshow(self):
        self.make()
        mpl.show()
    # makesave = make + save
    # -------------------------------
    def makesave(self,mode="png",filename="plot",includedate=False,res=150):
        self.make()
        save(mode=mode,filename=filename,includedate=includedate,res=res)
        close()


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
                 y=None,\
                 mapmode=False,\
                 proj=None,\
                 back=None,\
                 trans=1.0,\
                 vx=None,\
                 vy=None,\
                 c=None,\
                 blon=None,\
                 blat=None,\
                 area=None,\
                 sigma=None,\
                 vmin=None,\
                 vmax=None,\
                 showcb=True,\
                 wscale=None,\
                 svx=3,\
                 svy=3,\
                 leftcorrect=False,\
                 clev=None,\
                 cfmt=None,\
                 ccol="black",\
                 colorvec="black",\
                 *args, **kwargs):
        ## get initialization from parent class
        plot.__init__(self,*args, **kwargs)
        ## what could be defined by the user
        self.y = y
        self.mapmode = mapmode
        self.proj = proj
        self.back = back
        self.trans = trans
        self.vx = vx
        self.vy = vy
        self.colorvec = colorvec
        self.c = c
        self.blon = blon ; self.blat = blat
        self.area = area
        self.vmin = vmin ; self.vmax = vmax
        self.sigma = sigma
        self.showcb = showcb
        self.wscale = wscale
        self.svx = svx
        self.svy = svy
        self.leftcorrect = leftcorrect
        self.clev = clev
        self.cfmt = cfmt
        self.ccol = ccol

    # define_from_var
    # ... this uses settings in set_var.txt
    # -------------------------------
    def define_from_var(self):
        # get what is done in the parent class
        plot.define_from_var(self)
        # add specific stuff
        if self.var is not None:
         if self.var.upper() in vl.keys():
          self.colorbar = vc[self.var.upper()]
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
        # if projection is set, set mapmode to True
        if self.proj is not None:
            self.mapmode = True
        # set dummy xy axis if not defined
        if self.x is None: 
            self.x = np.array(range(self.f.shape[0]))
            self.mapmode = False
            print "!! WARNING !! dummy coordinates on x axis"
        if self.y is None: 
            self.y = np.array(range(self.f.shape[1]))
            self.mapmode = False
            print "!! WARNING !! dummy coordinates on y axis"
        # check sizes
        if self.c is not None:
            if self.c.ndim != 2:
                print "!! WARNING !! Contour is not a 2D field. No contour.",self.c.ndim
                self.c = None
        if self.f.ndim != 2:
            print "!! ERROR !! Field is not two-dimensional" ; exit()
        # transposing if necessary
        shape = self.f.shape
        if shape[0] != shape[1]:
         if len(self.x) == shape[0] and len(self.y) == shape[1]:
            #print "!! WARNING !! Transposing axes"
            self.f = np.transpose(self.f)
            if self.c is not None: 
              self.c = np.transpose(self.c)
        # bound field
        zevmin, zevmax = calculate_bounds(self.f,vmin=self.vmin,vmax=self.vmax,sigma=self.sigma)
        what_I_plot = bounds(self.f,zevmin,zevmax)
        # define contour field levels. define color palette
        ticks = self.div + 1
        zelevels = np.linspace(zevmin,zevmax,ticks)
        palette = get_cmap(name=self.colorbar)
        # do the same thing for possible contourline entries
        if self.c is not None:
            # if masked array, set masked values to filled values (e.g. np.nan) for plotting purposes
            if type(self.c).__name__ in 'MaskedArray':
               self.c[self.c.mask] = self.c.fill_value
            # set levels for contour lines
            if self.clev is None:
              zevminc, zevmaxc = calculate_bounds(self.c)
              what_I_contour = bounds(self.c,zevminc,zevmaxc)
              ticks = self.div + 1
              self.clev = np.linspace(zevminc,zevmaxc,ticks)
            else:
              what_I_contour = self.c
            # formatting
            ft = int(mpl.rcParams['font.size']*0.55)
            if self.cfmt is None: self.cfmt = "%.2g"
              
        ############################################################################################
        ### MAIN PLOT
        ### NB: contour lines are done before contour shades otherwise colorar error
        ############################################################################################
        if not self.mapmode:
            ## A SIMPLE 2D PLOT
            ###################
            # swapping if requested
            if self.swap:  x = self.y ; y = self.x
            else:          x = self.x ; y = self.y
            # coefficients on axis
            if self.xcoeff is not None: x=x*self.xcoeff
            if self.ycoeff is not None: y=y*self.ycoeff
            # make shaded and line contours
            if self.c is not None: 
                objC = mpl.contour(x, y, what_I_contour, \
                            self.clev, colors = self.ccol, linewidths = cline)
                ft = int(mpl.rcParams['font.size']*0.55)
                mpl.clabel(objC, inline=1, fontsize=ft,\
                             inline_spacing=1,fmt=self.cfmt)
            mpl.contourf(x, y, \
                         self.f, \
                         zelevels, cmap=palette)
            #mpl.pcolor(x,y,\
            #             self.f, \
            #             cmap=palette)
            # make log axes and/or invert ordinate
            ax = mpl.gca()
            if self.xdate:
              import matplotlib.dates as mdates
              ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%y %Hh'))
              #ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%b-%d %H:%M:%S'))
              ax.xaxis.set_major_locator(mdates.DayLocator())
              mpl.setp(mpl.xticks()[1], rotation=30, ha='right') # rotate the x labels
            if self.logx: mpl.semilogx()
            if self.logy: mpl.semilogy()
            if self.invert: ax.set_ylim(ax.get_ylim()[::-1])
            if self.xmin is not None and self.xmax is not None:
              if self.xmin > self.xmax: 
                ax.set_xlim(ax.get_xlim()[::-1])
                self.xmin,self.xmax = self.xmax,self.xmin
            if self.xmin is not None: ax.set_xbound(lower=self.xmin)
            if self.xmax is not None: ax.set_xbound(upper=self.xmax)
            if self.ymin is not None: ax.set_ybound(lower=self.ymin)
            if self.ymax is not None: ax.set_ybound(upper=self.ymax)
            # use back attributes to set a background
            if self.back is not None: ax.set_axis_bgcolor(self.back)
            # set the number of ticks
            if not self.logx:
                ax.xaxis.set_major_locator(MaxNLocator(self.nxticks))
            else:
                pass
                #print "!! WARNING. in logx mode, ticks are set automatically."
            if not self.logy:
                ax.yaxis.set_major_locator(MaxNLocator(self.nyticks))
            else:
                pass
                #print "!! WARNING. in logy mode, ticks are set automatically."
            ## specific modulo labels
            if self.modx is not None:
                ax = labelmodulo(ax,self.modx)
        else:
            ## A 2D MAP USING PROJECTIONS (basemap)
            #######################################
            mpl.xlabel("") ; mpl.ylabel("")
            # additional security in case self.proj is None here
            # ... we set cylindrical projection (the simplest one)
            if self.proj is None: self.proj = "cyl"
            # get lon and lat in 2D version.
            # (but first ensure we do have 2D coordinates)
            if self.x.ndim == 1:  [self.x,self.y] = np.meshgrid(self.x,self.y)
            elif self.x.ndim > 2: print "!! ERROR !! lon and lat arrays must be 1D or 2D"
            # get lat lon intervals and associated settings
            wlon = [np.min(self.x),np.max(self.x)]
            wlat = [np.min(self.y),np.max(self.y)]
            # -- area presets are in set_area.txt
            if self.area is not None:
             if self.area in area.keys():
                wlon, wlat = area[self.area]
            # -- user-defined limits
            if self.xmin is not None: wlon[0] = self.xmin
            if self.xmax is not None: wlon[1] = self.xmax
            if self.ymin is not None: wlat[0] = self.ymin
            if self.ymax is not None: wlat[1] = self.ymax
            # -- settings for meridians and parallels
            steplon = int(abs(wlon[1]-wlon[0])/3.)
            steplat = int(abs(wlat[1]-wlat[0])/3.)
            #mertab = np.r_[wlon[0]:wlon[1]:steplon] ; merlab = [0,0,0,1]
            #partab = np.r_[wlat[0]:wlat[1]:steplat] ; parlab = [1,0,0,0]
            if steplon < 1: steplon = 1
            if steplat < 1: steplat = 1
            if np.abs(wlon[0]) < 180.1 and np.abs(wlon[1]) < 180.1:
                mertab = np.r_[-180.:180.:steplon]
            else:
                mertab = np.r_[0.:360.:steplon]
            merlab = [0,0,0,1]
            partab = np.r_[-90.:90.+steplat:steplat] ; parlab = [1,0,0,0]
            format = '%.1f'
            # -- center of domain and bounding lats
            lon_0 = 0.5*(wlon[0]+wlon[1])
            lat_0 = 0.5*(wlat[0]+wlat[1])
            # some tests, bug fixes, and good-looking settings
            # ... cyl is good for global and regional
            if self.proj == "cyl":
                format = '%.0f'
                partab = np.r_[-90.:90.+15.:15.]
            # ... global projections
            elif self.proj in ["ortho","moll","robin"]:
                wlat[0] = None ; wlat[1] = None ; wlon[0] = None ; wlon[1] = None
                lon_0 = np.ceil(lon_0) # reverse map if lon_0 is slightly below 180 with [0,360]
                steplon = 30. ; steplat = 30.
                if self.proj in ["moll"]: steplon = 60.
                if self.proj in ["robin"]: steplon = 90.
                mertab = np.r_[-360.:360.:steplon]
                #partab = np.r_[-90.:90.+steplat:steplat]
                partab = np.r_[-60.,-30.,0.,30.,60.]
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
	    try:
	      from mpl_toolkits.basemap import Basemap
	    except:
              print "!! ERROR !! basemap is not available."
	      print "... either install it or use another plot type."
	      exit()
            m = Basemap(projection=self.proj,\
                        lat_0=lat_0,lon_0=lon_0,\
                        boundinglat=self.blat,\
                        llcrnrlat=wlat[0],urcrnrlat=wlat[1],\
                        llcrnrlon=wlon[0],urcrnrlon=wlon[1])
            # in some case need to translated to the left for colorbar + labels
            # TBD: break stuff. a better solution should be found.
            if self.leftcorrect:
                ax = mpl.gca()
                pos = ax.get_position().bounds
                newpos = [0.,pos[1],pos[2],pos[3]]
                ax.set_position(newpos)
            # draw meridians and parallels
            ft = int(mpl.rcParams['font.size']*3./4.)
            zelatmax = 85.
            m.drawmeridians(mertab,labels=merlab,color='grey',linewidth=0.75,fontsize=ft,fmt=format,latmax=zelatmax)
            m.drawparallels(partab,labels=parlab,color='grey',linewidth=0.75,fontsize=ft,fmt=format,latmax=zelatmax)
            # define background (see set_back.txt)
            if self.back is not None:
              if self.back in back.keys():
                 print "**** info: loading a background, please wait.",self.back
                 if self.back not in ["coast","sea"]:
                    try: m.warpimage(back[self.back],scale=0.75)
                    except: print "!! ERROR !! no background image could be loaded. probably not connected to the internet?"
                 elif self.back == "coast":
                    m.drawcoastlines()
                 elif self.back == "sea":
                    m.drawlsmask(land_color='white',ocean_color='aqua')
              else:
                 print "!! ERROR !! requested background not defined. change name or fill in set_back.txt" ; exit()
            # define x and y given the projection
            x, y = m(self.x, self.y)
            # contour field. first line contour then shaded contour.
            if self.c is not None: 
                #zelevelsc = np.arange(900.,1100.,5.)
                objC2 = m.contour(x, y, what_I_contour, \
                            self.clev, colors = self.ccol, linewidths = cline)
                #mpl.clabel(objC2, inline=1, fontsize=10,manual=True,fmt='-%2.0f$^{\circ}$C',colors='r')
                #mpl.clabel(objC2, inline=0, fontsize=8, fmt='%.0f',colors='r', inline_spacing=0) 
            m.contourf(x, y, what_I_plot, zelevels, cmap = palette, alpha = self.trans, antialiased=True)
        ############################################################################################
        ### COLORBAR
        ############################################################################################
        if self.trans > 0. and self.showcb:
            ## draw colorbar. settings are different with projections. or if not mapmode.
            #if not self.mapmode: orientation=zeorientation ; frac = 0.075 ; pad = 0.03 ; lu = 0.5
            if not self.mapmode: orientation=zeorientation ; frac = 0.15 ; pad = 0.04 ; lu = 0.5
            elif self.proj in ['moll']: orientation="horizontal" ; frac = 0.08 ; pad = 0.03 ; lu = 1.0
            elif self.proj in ['robin']: orientation="horizontal" ; frac = 0.07 ; pad = 0.1 ; lu = 1.0
            elif self.proj in ['cyl']: orientation="vertical" ; frac = 0.023 ; pad = 0.03 ; lu = 0.5
            else: orientation = zeorientation ; frac = zefrac ; pad = 0.03 ; lu = 0.5
            if self.cbticks is None:
                self.cbticks = min([ticks/2+1,21])
            zelevpal = np.linspace(zevmin,zevmax,num=self.cbticks)
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
        ### not expecting NaN in self.vx and self.vy. masked arrays is just enough.
        if self.vx is not None and self.vy is not None: 
                # vectors on map projection or simple 2D mapping
                if self.mapmode: 
                   try:
                     #[vecx,vecy] = m.rotate_vector(self.vx,self.vy,self.x,self.y) # for metwinds only ?
                     vecx,vecy = self.vx,self.vy
                   except:
                     print "!! ERROR !! Problem with field shapes for vector?" 
                     print self.vx.shape,self.vy.shape,self.x.shape,self.y.shape
                     exit()
                else:
                   vecx,vecy = self.vx,self.vy 
                   if x.ndim < 2 and y.ndim < 2: x,y = np.meshgrid(x,y)
                # reference vector is scaled
                if self.wscale is None:
                    self.wscale = ppcompute.mean(np.sqrt(self.vx*self.vx+self.vy*self.vy))
                # make vector field
                if self.mapmode: 
                    q = m.quiver( x[::self.svy,::self.svx],y[::self.svy,::self.svx],\
                                  vecx[::self.svy,::self.svx],vecy[::self.svy,::self.svx],\
                                  angles='xy',color=self.colorvec,pivot='middle',\
                                  scale=self.wscale*reducevec,width=widthvec )
                else:
                    q = mpl.quiver( x[::self.svy,::self.svx],y[::self.svy,::self.svx],\
                                    vecx[::self.svy,::self.svx],vecy[::self.svy,::self.svx],\
                                    angles='xy',color=self.colorvec,pivot='middle',\
                                    scale=self.wscale*reducevec,width=widthvec )
                # make vector key.
                #keyh = 1.025 ; keyv = 1.05 # upper right corner over colorbar
                keyh = 0.97 ; keyv = 1.06
                keyh = 0.97 ; keyv = 1.11
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

    # makeshow = make + show
    # -------------------------------
    def makeshow(self,grid="off"):
        self.make()
        if grid != "off": mpl.grid(color=grid)
        mpl.show()

    # makesave = make + save
    # -------------------------------
    def makesave(self,mode="png",filename="plot",includedate=False,res=150):
        self.make()
        save(mode=mode,filename=filename,includedate=includedate,res=res)
        close()


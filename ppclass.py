###############################################
## PLANETOPLOT                               ##
## --> PPCLASS                               ##
## A generic and versatile Python module     ##
## ... to read netCDF files and plot         ##
###############################################
## Author: Aymeric Spiga. 02-03/2013         ##
###############################################
# python built-in librairies
import os
import time as timelib
import pickle
# added librairies
import numpy as np
import netCDF4
import matplotlib.pyplot as mpl
# personal librairies
import ppplot
import ppcompute
###############################################

###################################
#### HEADER                      ##
#### ... executed when imported  ##
###################################
print "************************************"
print "**** WELCOME TO PLANETOPLOT 2.0 ****"
print "************************************"
# where settings files are located...
whereset = None
whereset = ppcompute.findset(whereset)
# ... we load user-defined automatic settings from set_ppclass.txt
zefile = "set_ppclass.txt"
glob_listx = [] ; glob_listy = [] ; glob_listz = [] ; glob_listt = []
glob_listarea = []
try: 
    f = open(whereset+zefile, 'r') ; lines = f.readlines()
    for stuff in lines[5].strip().split(';'): glob_listx.append(stuff)
    for stuff in lines[8].strip().split(';'): glob_listy.append(stuff)
    for stuff in lines[11].strip().split(';'): glob_listz.append(stuff)
    for stuff in lines[14].strip().split(';'): glob_listt.append(stuff)
    for stuff in lines[17].strip().split(';'): glob_listarea.append(stuff)
except IOError: 
    print "warning: "+zefile+" not in "+whereset+" ; no presets."

##################################
#### USEFUL GENERIC FUNCTIONS ####
##################################

# inspect variables and dimensions in a netCDF file
def inspect(filename):
    print "**** INSPECT FILE",filename
    test = netCDF4.Dataset(filename)
    print "**** VARIABLES: ",test.variables.keys()
    for dim in test.dimensions.keys():
        output = "**** DIMENSION: "+str(dim)+" "+str(len(test.dimensions[dim]))
        try: output = output + " ----> "+str(test.variables[dim][0])+"  "+str(test.variables[dim][-1])
        except: pass
        print output ; output = ""

# check a tab and exit if wrong. if just one string make it a list.
# (if allownumber, convert this into a string).
def checktab(tab,mess="",allownone=False,allownumber=False):
    if tab is None: 
      if not allownone:  print "pp.define: no "+mess ; exit()
      else: pass
    else:
      if not isinstance(tab, list):
        if isinstance(tab, str): 
            tab = [tab]
        elif (isinstance(tab, int) or isinstance(tab, float)) and allownumber: 
            tab = [str(tab)] 
        else: 
            print "pp.define: "+mess+" should be either a string or a list of strings!" ; exit()
      elif isinstance(tab, list):
        if isinstance(tab[0],str): 
            pass
        elif (isinstance(tab[0], int) or isinstance(tab[0], float)) and allownumber:
            for iii in range(len(tab)): tab[iii] = str(tab[iii])
        else: 
            print "pp.define: "+mess+" should be either a string or a list of strings!" ; exit()
    return tab

# determine which method is to be applied to a given dimension
def findmethod(tab):
    if tab is None:              output = "free"
    elif tab[0,0] != tab[0,1]:   output = "comp"
    else:                        output = "fixed"
    return output

# read what is given by the user (version of T. Navarro)
def readslices(saxis):
    if saxis == None:
        zesaxis = None
    else:
        zesaxis = np.empty((len(saxis),2))
        for i in range(len(saxis)):
            a = separatenames(saxis[i])
            if len(a) == 1:
                zesaxis[i,:] = float(a[0])
            else:
                zesaxis[i,0] = float(a[0])
                zesaxis[i,1] = float(a[1])          
    return zesaxis
# (needed by readslices)
# look for comas in the input name to separate different names (files, variables,etc ..)
def separatenames (name):
    if name is None: names = None
    else:
      names = [] ; stop = 0 ; currentname = name
      while stop == 0:
        indexvir = currentname.find(',')
        if indexvir == -1: stop = 1 ; name1 = currentname
        else: name1 = currentname[0:indexvir]
        names = np.concatenate((names,[name1]))
        currentname = currentname[indexvir+1:len(currentname)]
    return names

#######################
### THE MAIN OBJECT ###
#######################
class pp():

    # print out a help string when help is invoked on the object
    def __repr__(self):
        whatprint = 'pp object. \"help(pp)\" for more information\n'
        return whatprint

    # default settings
    # -- user can define settings by two methods. 
    # -- 1. yeah = pp(file="file.nc")
    # -- 2. yeah = pp() ; yeah.file = "file.nc"
    def __init__(self,file=None,var="notset",\
                      filegoal=None,vargoal=None,\
                      x=None,y=None,z=None,t=None,\
                      stridex=1,stridey=1,\
                      stridez=1,stridet=1,\
                      compute="mean",\
                      verbose=False,noproj=False,\
                      superpose=False,\
                      plotin=None,\
                      forcedimplot=-1,\
                      out="gui",\
                      filename="myplot",\
                      folder="./"):
        self.request = None
        self.nfin = 0 ; self.nvin = 0
        self.nplotx = None ; self.nploty = None
        self.nplotz = None ; self.nplott = None
        self.status = "init"
        self.fig = None ; self.subv = None ; self.subh = None 
        self.n = 0 ; self.howmanyplots = 0
        self.nplot = 0
        self.p = None
        self.customplot = False
        ## what could be defined by the user
        self.file = file
        self.var = var
        self.filegoal = filegoal
        self.vargoal = vargoal
        self.x = x ; self.y = y   ## if None, free dimension
        self.z = z ; self.t = t   ## if None, free dimension
        self.stridex = stridex ; self.stridey = stridey
        self.stridez = stridez ; self.stridet = stridet
        self.compute = compute
        self.verbose = verbose
        self.noproj = noproj
        self.plotin = plotin
        self.superpose = superpose
        self.forcedimplot = forcedimplot
        self.out = out
        self.filename = filename
        self.folder = folder

    # print status
    def printstatus(self):
        if self.filename == "THIS_IS_A_CLONE":
            pass
        else:
            print "**** Done step: " + self.status

    #####################################################
    # EMULATE OPERATORS + - * / ** << FOR PP() OBJECTS  #
    #####################################################

    # check the compatibility of two objects for operations
    # --> if other is a pp class, test sizes and return isnum = False
    # --> if other is an int or a float, return isnum = True
    # --> otherwise, just print an error and exit
    def checktwo(self,other):
        if other.__class__.__name__ == "pp":
          isnum = False
          if self.status in ["init","defined"] or other.status in ["init","define"]: 
             print "!! ERROR !! Please use .retrieve to get fields for plots with one of your pp operands." ; exit()
          if self.nfin   != other.nfin   or \
             self.nvin   != other.nvin   or \
             self.nplott != other.nplott or \
             self.nplotz != other.nploty or \
             self.nploty != other.nploty or \
             self.nplotx != other.nplotx :
               print "!! ERROR !! The two operands do not have the same number of files, variables, t z y x requests."
               exit()
        elif isinstance(other,int) or isinstance(other,float):
          isnum = True
        else:
          print "!! ERROR !! The operand is neither a pp class nor an integer or a float." ; exit()
        return isnum

    # define a selective copy of a pp() object for operations
    # ... copy.copy() is not conservative (still acts like a pointer)
    # ... copy.deepcopy() does not work with netCDF objects
    # so what is done here is a copy of everything except
    # (to avoid sharing with self and therefore modifying self through operations)
    # - request attribute of pp() object
    # - field attribute of the onerequest() objects
    def selective_copy(self):
        if self.status in ["init","defined"]:
            print "!! ERROR !! Please use .retrieve to get fields for the object you want to copy from." ; exit()
        the_clone = pp()
        for k, v in vars(self).items():
           if k != "request":
               setattr(the_clone,k,v)
        the_clone.verbose = False
        the_clone.filename = "THIS_IS_A_CLONE" # trick to avoid additional outputs
        the_clone.define()
        for i in range(self.nfin):
         for j in range(self.nvin):
          for t in range(self.nplott):
           for z in range(self.nplotz):
            for y in range(self.nploty):
             for x in range(self.nplotx):
              obj_ref = self.request[i][j][t][z][y][x]
              obj = the_clone.request[i][j][t][z][y][x]
              for k, v in vars(obj_ref).items():
               if k != "field":
                setattr(obj,k,v)
        the_clone.status = "retrieved"
        return the_clone

    # define the operation + on two objects. or with an int/float.
    # ... with selective_copy the self object is not modified.
    def __add__(self,other):
        isnum = self.checktwo(other)
        the_clone = self.selective_copy()
        for i in range(self.nfin):
         for j in range(self.nvin):
          for t in range(self.nplott):
           for z in range(self.nplotz):
            for y in range(self.nploty):
             for x in range(self.nplotx):
              obj = the_clone.request[i][j][t][z][y][x]
              obj_ref = self.request[i][j][t][z][y][x]
              if not isnum:   
                  ope = other.request[i][j][t][z][y][x].field
                  if obj_ref.field.shape != ope.shape:
                    print "!! ERROR !! The two fields for operation do not have the same shape.",self.field.shape,other.field.shape
                    exit()
              else:           
                  ope = other
              if "vector" in self.vargoal[j] + self.filegoal[i]:
                  print "!! ERROR !! we do not operate on vectors yet."
                  exit()
              else:
                  obj.field = obj_ref.field + ope
        return the_clone

    # define the operation - on two objects. or with an int/float.
    # ... with selective_copy the self object is not modified.
    def __sub__(self,other):
        isnum = self.checktwo(other)
        the_clone = self.selective_copy()
        for i in range(self.nfin):
         for j in range(self.nvin):
          for t in range(self.nplott):
           for z in range(self.nplotz):
            for y in range(self.nploty):
             for x in range(self.nplotx):
              obj = the_clone.request[i][j][t][z][y][x]
              obj_ref = self.request[i][j][t][z][y][x]
              if not isnum:
                  ope = other.request[i][j][t][z][y][x].field
                  if obj_ref.field.shape != ope.shape:
                    print "!! ERROR !! The two fields for operation do not have the same shape.",self.field.shape,other.field.shape
                    exit()
              else:
                  ope = other
              if "vector" in self.vargoal[j] + self.filegoal[i]:
                  print "!! ERROR !! we do not operate on vectors yet."
                  exit()
              else:
                  obj.field = obj_ref.field - ope
        return the_clone

    # define the operation * on two objects. or with an int/float.
    # ... with selective_copy the self object is not modified.
    def __mul__(self,other):
        isnum = self.checktwo(other)
        the_clone = self.selective_copy()
        for i in range(self.nfin):
         for j in range(self.nvin):
          for t in range(self.nplott):
           for z in range(self.nplotz):
            for y in range(self.nploty):
             for x in range(self.nplotx):
              obj = the_clone.request[i][j][t][z][y][x]
              obj_ref = self.request[i][j][t][z][y][x]
              if not isnum:
                  ope = other.request[i][j][t][z][y][x].field
                  if obj_ref.field.shape != ope.shape:
                    print "!! ERROR !! The two fields for operation do not have the same shape.",self.field.shape,other.field.shape
                    exit()
              else:
                  ope = other
              if "vector" in self.vargoal[j] + self.filegoal[i]:
                  print "!! ERROR !! we do not operate on vectors yet."
                  exit()
              else:
                  obj.field = obj_ref.field * ope
        return the_clone

    # define the operation / on two objects. or with an int/float.
    # ... with selective_copy the self object is not modified.
    def __div__(self,other):
        isnum = self.checktwo(other)
        the_clone = self.selective_copy()
        for i in range(self.nfin):
         for j in range(self.nvin):
          for t in range(self.nplott):
           for z in range(self.nplotz):
            for y in range(self.nploty):
             for x in range(self.nplotx):
              obj = the_clone.request[i][j][t][z][y][x]
              obj_ref = self.request[i][j][t][z][y][x]
              if not isnum:
                  ope = other.request[i][j][t][z][y][x].field
                  if obj_ref.field.shape != ope.shape:
                    print "!! ERROR !! The two fields for operation do not have the same shape.",self.field.shape,other.field.shape
                    exit()
              else:
                  ope = other
              if "vector" in self.vargoal[j] + self.filegoal[i]:
                  print "!! ERROR !! we do not operate on vectors yet."
                  exit()
              else:
                  obj.field = obj_ref.field / ope
        return the_clone

    # define the reverse operation float/int + object
    def __radd__(self,other):
        isnum = self.checktwo(other)
        if not isnum: print "!! ERROR !! Operand should be a number" ; exit()
        return self.__add__(other)

    # define the reverse operation float/int - object
    def __rsub__(self,other):
        isnum = self.checktwo(other)
        if not isnum: print "!! ERROR !! Operand should be a number" ; exit()
        return self.__sub__(other)

    # define the reverse operation float/int * object
    def __rmul__(self,other):
        isnum = self.checktwo(other)
        if not isnum: print "!! ERROR !! Operand should be a number" ; exit()
        return self.__mul__(other)

    # define the reverse operation float/int / object
    def __rdiv__(self,other):
        isnum = self.checktwo(other)
        if not isnum: print "!! ERROR !! Operand should be a number" ; exit()
        return self.__div__(other)

    # define the operation ** on one object.
    # ... with selective_copy the self object is not modified.
    def __pow__(self,num):
        the_clone = self.selective_copy()
        if isinstance(num,int) or isinstance(num,float):
            for i in range(self.nfin):
             for j in range(self.nvin):
              for t in range(self.nplott):
               for z in range(self.nplotz):
                for y in range(self.nploty):
                 for x in range(self.nplotx):
                  obj  = the_clone.request[i][j][t][z][y][x]
                  obj_ref = self.request[i][j][t][z][y][x]
                  if "vector" in self.vargoal[j] + self.filegoal[i]:
                      print "!! ERROR !! we do not operate on vectors yet."
                      exit()
                  else:
                      obj.field = obj_ref.field ** num
        else:
            print "!! ERROR !! To define a power, either an int or a float is needed." ; exit()
        return the_clone

    # define the operation <<
    # ... e.g. obj2 << obj1
    # ... means: get init for pp object obj2 from another pp object obj1
    # ... (this should solve the affectation trap obj2 = obj1)
    def __lshift__(self,other):
        if other.__class__.__name__ == "pp":
            self.file = other.file
            self.var = other.var
            self.filegoal = other.filegoal
            self.vargoal = other.vargoal
            self.x = other.x ; self.y = other.y   ## if None, free dimension
            self.z = other.z ; self.t = other.t   ## if None, free dimension
            self.stridex = other.stridex ; self.stridey = other.stridey
            self.stridez = other.stridez ; self.stridet = other.stridet
            self.verbose = other.verbose
            self.noproj = other.noproj
            self.plotin = other.plotin
            self.superpose = other.superpose
            self.forcedimplot = other.forcedimplot
            self.out = other.out
            self.filename = other.filename
            self.folder = other.folder
        else:
            print "!! ERROR !! argument must be a pp object." ; exit()

    ##############################################################################################
    # define method
    # ---------
    # ... (file and var are either one string or a vector of strings)
    # ... the goal of define is to define a 2D array of onerequest() objects (see class below)
    #     given the number of file, var, x, y, z, t asked by the user
    # ... objectives for file or var are given through filegoal and vargoal
    #     --> possible values: main contour vector
    # ---------
    # ... then onerequest() objects are being defined more precisely
    #     by getting index_x index_y index_z index_t
    #     and setting method_x method_y method_z method_t to either
    #      - "free" for free dimensions (plot dimensions) 
    #      - "comp" for averages, max, min
    #      - "fixed" for fixed dimensions (possibly several i.e. multislice)
    ##############################################################################################
    def define(self):
        self.printstatus()
        # initial check and get dimensions
        self.file = checktab(self.file,mess="file")
        self.nfin = len(self.file)
        if self.verbose:
            for i in range(self.nfin): inspect(self.file[i])
        self.var = checktab(self.var,mess="var")
        self.nvin = len(self.var)
        # check goal tabs for files and variables
        # ... default is to plot everything
        if self.filegoal is None: self.filegoal = ["main"]*self.nfin
        if self.vargoal is None:  self.vargoal  = ["main"]*self.nvin
        self.filegoal = checktab(self.filegoal, mess="filegoal")
        self.vargoal  = checktab(self.vargoal,  mess="vargoal")
        if len(self.filegoal) != self.nfin:  print "!! ERROR !! filegoal must be the same size as file." ; exit()
        if len(self.vargoal)  != self.nvin:  print "!! ERROR !! vargoal must be the same size as var." ; exit()
        # variables: initial check
        self.x = checktab(self.x,mess="x",allownone=True,allownumber=True)
        self.y = checktab(self.y,mess="y",allownone=True,allownumber=True)
        self.z = checktab(self.z,mess="z",allownone=True,allownumber=True)
        self.t = checktab(self.t,mess="t",allownone=True,allownumber=True)
        # for the moment not var- nor file- dependent. 
        # but this could be the case.
        sx = readslices(self.x) ; sy = readslices(self.y)
        sz = readslices(self.z) ; st = readslices(self.t)
        # get methods
        mx = findmethod(sx) ; my = findmethod(sy)
        mz = findmethod(sz) ; mt = findmethod(st)
        # get number of plots to be done
        if mx == "fixed": self.nplotx = sx.size/2
        else:             self.nplotx = 1
        if my == "fixed": self.nploty = sy.size/2
        else:             self.nploty = 1
        if mz == "fixed": self.nplotz = sz.size/2
        else:             self.nplotz = 1
        if mt == "fixed": self.nplott = st.size/2
        else:             self.nplott = 1
        if self.verbose:  print "**** OK. Plots over x,y,z,t -->",self.nplotx,self.nploty,self.nplotz,self.nplott
        # create the list of onerequest() objects
        self.request = [[[[[[ \
                       onerequest() \
                       for x in range(self.nplotx)] for y in range(self.nploty)] \
                       for z in range(self.nplotz)] for t in range(self.nplott)] \
                       for j in range(self.nvin)]   for i in range(self.nfin)] 
        # loop on onerequest() objects
        for i in range(self.nfin):
         for j in range(self.nvin):
          for t in range(self.nplott):
           for z in range(self.nplotz):
            for y in range(self.nploty):
             for x in range(self.nplotx):
              obj = self.request[i][j][t][z][y][x]
              # fill in names for files and variables
              obj.verbose = self.verbose
              obj.file = self.file[i]
              obj.var = self.var[j]
              # indicate the computation method
              obj.compute = self.compute
              # open the files (the same file might be opened several times but this is cheap)
              obj.openfile()
              ### get x,y,z,t dimensions from file
              obj.getdim()
              ### get methods 
              obj.method_x = mx ; obj.method_y = my
              obj.method_z = mz ; obj.method_t = mt           
              ### get index
              obj.getindextime(dalist=st,ind=t,stride=self.stridet)
              obj.getindexvert(dalist=sz,ind=z,stride=self.stridez)
              obj.getindexhori(dalistx=sx,dalisty=sy,indx=x,indy=y,stridex=self.stridex,stridey=self.stridey)
        # change status
        self.status = "defined"

    ##############################################################################################
    # retrieve method
    # --> for each element onerequest() in the array, get field .var from .f file
    # --> see below the onerequest() class: 
    #        - only get what is needed for computing and plotting
    #        - averages etc... are computed here
    # --> RESULT: each onerequest() object has now its attribute .field filled
    # --> if one wants to perform operations on fields, this should be done after retrieve()
    ##############################################################################################
    def retrieve(self):
        self.printstatus()
        # check if things were done OK before
        if self.status != "defined": print "!! ERROR !! Please use .define() to define your pp object." ; exit()
        ## first get fields
        ## ... only what is needed is extracted from the files
        ## ... and computations are performed
        for i in range(self.nfin):
         for j in range(self.nvin):
          for t in range(self.nplott):
           for z in range(self.nplotz):
            for y in range(self.nploty):
             for x in range(self.nplotx):
              obj = self.request[i][j][t][z][y][x]
              obj.getfield()
              obj.computations()
        # change status
        self.status = "retrieved"

    ##########################################################
    # get: a shortcut method for the define + retrieve chain #
    ##########################################################
    def get(self):
        self.define()
        self.retrieve()

    ########################################
    # smooth: smooth the field in 1D or 2D #
    ########################################
    ## TBD: smooth not OK with masked array in the end of retrieve()
    def smooth(self,window):
        if self.verbose: 
            print "!! WARNING !! Performing a smoothing with a window size",window
            print "!! WARNING !! To come back to unsmoothed file, use .get() again"
        for i in range(self.nfin):
         for j in range(self.nvin):
          for t in range(self.nplott):
           for z in range(self.nplotz):
            for y in range(self.nploty):
             for x in range(self.nplotx):
              obj = self.request[i][j][t][z][y][x]
              if obj.field.ndim == 1:
                  print "!! ERROR !! 1D smoothing not supported yet because reduces array sizes."
                  exit()
                  # TBD TBD TBD
                  #obj.field = ppcompute.smooth1d(obj.field,window=window)
              elif obj.field.ndim == 2:
                  obj.field = ppcompute.smooth2d(obj.field,window=window)

    ##############################################################################################  
    # defineplot method
    # --> defineplot first defines arrays of plot objects and set each of them
    #     ... simple looping except cases where goal is not main (e.g. contour or vector)
    # --> principle: each onerequest() object gives birth to a subplot
    # --> defineplot vs. makeplot: defining plot and actually plotting it are clearly separated
    # --> THE KEY OUPUT OF defineplot IS AN ARRAY self.p OF PLOT OBJECTS
    # optional arguments
    # --> extraplot: to indicate a number of plots to be added afterwards (use self.plotin)
    # --> loadfile: to use self.p from a previously saved file
    ##############################################################################################
    def defineplot(self,extraplot=0,loadfile=None):
        # -----------------------------------------------------
        # LOAD MODE: load a self.p object. count plots from it.
        # -----------------------------------------------------
        if loadfile is not None:
            try: filehandler = open(loadfile, 'r') ; self.p = pickle.load(filehandler) 
            except IOError: print "!! ERROR !! Cannot find object file to load." ; exit()
            self.status = "definedplot" ; self.plotin = None
            self.nplot = len(self.p) ; self.howmanyplots = self.nplot
            return
        # -----------------------------------------------------
        # REGULAR MODE
        # -----------------------------------------------------
        self.printstatus()
        # check if things were done OK before
        if self.status in ["init","defined"]: 
            print "!! ERROR !! Please use .retrieve() to get fields for plots with your pp object." ; exit()
        # check self.plotin (an existing fig on which to add plots afterwards)
        if self.plotin.__class__.__name__ == "pp":
            if self.plotin.fig is None:
                self.plotin = None # this is an additional security in case 
                                   #   a pp object is given without figure opened yet.
        elif self.plotin is not None:
            print "!! ERROR !! plotin argument must be a pp object." ; exit()
        # initialize the array of subplot objects
        # either something new or attributes coming from plotin object
        if self.plotin is None:  self.p = [ ]
        else:                    self.p = self.plotin.p
        # create an array of subplot objects
        # ... in theory the order of looping can be changed without any harm
        # ... the only important thing is to keep i,j,t,z,y,x resp. for file,var,t,z,y,x
        count = 0
        for i in range(self.nfin):
         if self.filegoal[i] == "main": 
          for j in range(self.nvin):
           if self.vargoal[j] == "main":
            for t in range(self.nplott):
             for z in range(self.nplotz):
              for y in range(self.nploty):
               for x in range(self.nplotx):
                # look at dimension and append the right plot object
                obj = self.request[i][j][t][z][y][x]
                dp = obj.dimplot
                if dp == 1 or self.forcedimplot == 1:    plobj = ppplot.plot1d()
                elif dp == 2 or self.forcedimplot == 2:  plobj = ppplot.plot2d()
                elif dp == 0: print "**** OK. VALUES VALUES VALUES",obj.field
                else:         print "!! ERROR !! 3D or 4D plots not supported" ; exit()
                # load abscissa and ordinate in obj
                obj.definecoord()
                # we start to define things here before appending
                # (convenient: could be overridden by the user before makeplot)
                # ... the if loop is necessary so that we can loop above on the dp=0 case
                if dp in [1,2]:
                    # and define what to do in plobj
                    plobj.invert = obj.invert_axes
                    plobj.swap = obj.swap_axes
                    # axis labels
                    plobj.xlabel = obj.absclab ; plobj.ylabel = obj.ordilab
                    # superpose or not (this is mostly for saving purpose)
                    plobj.superpose = self.superpose
                    # get title, colormaps, labels, etc.. from var
                    plobj.var = obj.var 
                    plobj.define_from_var()
                    # generic 1D/2D: load field and coord in plot object 
                    plobj.field = obj.field    # field to be plotted
                    plobj.absc = obj.absc      # abscissa (or longitude)
                    plobj.ordi = obj.ordi      # ordinate (or latitude)
                                               # -- useless in 1D but not used anyway
                    if dp == 2:
                    # specific 2d plot stuff
                        # -- light grey background for missing values
                        if type(plobj.field).__name__ in 'MaskedArray': plobj.axisbg = '0.75'
                        # -- define if it is a map or a plot
                        plobj.mapmode = ( obj.method_x+obj.method_y == "freefree" \
                                       and "grid points" not in obj.name_x \
                                       and not self.noproj )
                    # finally append
                    self.p.append(plobj)
                    count = count + 1
        # self.nplot is number of plot to be defined in this call to defineplot()
        # (because of self.plotin this might less than length of self.p)
        self.nplot = count
        # --- superimposed contours and vectors ---
        # we have to start another loop because we need forward information
        # TBD: there is probably a more flexible way to do that
        count = 0
        for i in range(self.nfin):
         for j in range(self.nvin):
          for t in range(self.nplott):
           for z in range(self.nplotz):
            for y in range(self.nploty):
             for x in range(self.nplotx):
              goal = self.vargoal[j] + self.filegoal[i]
              obj = self.request[i][j][t][z][y][x]
              if "mainmain" in goal and obj.dimplot == 2:
                  # the plot object we consider in the loop
                  pl = self.p[count]
                  # -- see if there is a contour requested...
                  # (we use try because we might be at the end of the list)
                  found = 0
                  try:    condvar = self.vargoal[j+1]
                  except: condvar = "itisok"
                  try:    condfile = self.filegoal[i+1]
                  except: condfile = "itisok"
                  # ... get contour
                  ##########################################
                  # NB: contour is expected to be right after main otherwise it is not displayed
                  ##########################################
                  if condvar == "contour":
                      plobj.addcontour = self.request[i][j+1][t][z][y][x].field ; found += 1
                  if condfile == "contour":
                      plobj.addcontour = self.request[i+1][j][t][z][y][x].field ; found += 1
                  # see if there is a vector requested...
                  # (we use try because we might be at the end of the list)
                  try:    condvar = self.vargoal[j+found+1]+self.vargoal[j+found+2]
                  except: condvar = "itisok"
                  try:    condfile = self.filegoal[i+found+1]+self.filegoal[i+found+2]
                  except: condfile = "itisok"
                  # ... get vector and go directly to the next iteration
                  # (in some cases we would do this twice but this is cheap)
                  if "vector" in condvar:
                      plobj.addvecx = self.request[i][j+found+1][t][z][y][x].field
                      plobj.addvecy = self.request[i][j+found+2][t][z][y][x].field
                  if "vector" in condfile:
                      plobj.addvecx = self.request[i+found+1][j][t][z][y][x].field
                      plobj.addvecy = self.request[i+found+2][j][t][z][y][x].field
                  count = count + 1
        # COUNT PLOTS. if 0 just exit.
        # self.howmanyplots is self.nplot + possible extraplots 
        self.howmanyplots = self.nplot + extraplot
        if self.howmanyplots > 0:
            if self.verbose: print "**** OK. expect %i plots" % (self.howmanyplots)
        else:
            exit() # because this means that we only had 0D values !
        # final status
        self.status = "definedplot"

    ##############################################################################################
    # makeplot method
    # --> after defineplot and before makeplot, user-defined plot settings can be easily given
    #     simply by modifying the attributes of each elements of self.p
    # --> to change only one plot setting, no need to call defineplot again before makeplot
    # --> in the end, the array self.p of plot objects is saved for easy and convenient replotting
    # --> see practical examples in the folder 'examples'
    ##############################################################################################
    def makeplot(self):
        self.printstatus()
        # a few initial operations
        # ------------------------
        if "definedplot" not in self.status: 
            print "!! ERROR !! Please use .defineplot() before .makeplot() can be used with your pp object." ; exit()
        # open a figure and define subplots          
        # ---------------------------------
        if self.plotin is None:  
            # start from scratch
            self.fig = mpl.figure(figsize=(16,8))
            self.subv,self.subh = ppplot.definesubplot(self.howmanyplots,self.fig) 
            self.n = 0
            ## adapted space for labels etc
            ## ... except for ortho because there is no label anyway
            self.customplot = self.p[0].field.ndim == 2 \
                        and self.p[0].mapmode == True \
                        and self.p[0].proj not in ["ortho"]
            if self.customplot:
                margin = 0.07
                self.fig.subplots_adjust(left=margin,right=1-margin,bottom=margin,top=1-margin)
        else:
            # start from an existing figure.
            # extraplot must have been set in the call to the previous figure.
            self.fig = self.plotin.fig
            self.subv,self.subh = self.plotin.subv,self.plotin.subh
            self.n = self.plotin.n
            self.howmanyplots = self.plotin.howmanyplots
            self.customplot = self.plotin.customplot
        # LOOP on all subplots
        # NB: cannot use 'for pl in self.p' if self.plotin not None
        # --------------------
        for count in range(self.nplot):
            # the plot object we consider in the loop
            pl = self.p[self.n]
            # before making the plot, create a subplot. the first one is numbered 1 not 0.
            # ... if pl.superpose, we use only one and only figure
            # ... (and we have to be careful with not doing things several times)
            if pl.superpose:
                if self.n == 0: 
                    self.fig.add_subplot(1,1,1,axisbg=pl.axisbg) # define one subplot (still needed for user-defined font sizes)
                    sav = pl.xlabel,pl.ylabel,pl.xcoeff,pl.ycoeff,pl.title,pl.swaplab # save titles and labels
                else: 
                    pl.invert = False ; pl.lstyle = None # don't invert again axis
                    # set saved titles and labels
                    if self.plotin is None:
                        pl.xlabel,pl.ylabel,pl.xcoeff,pl.ycoeff,pl.title,pl.swaplab = sav 
                    else:
                        prev_plot = self.plotin.p[self.n-1]
                        pl.xlabel = prev_plot.xlabel
                        pl.ylabel = prev_plot.ylabel
                        pl.xcoeff = prev_plot.xcoeff
                        pl.ycoeff = prev_plot.ycoeff
                        pl.title = prev_plot.title
                        pl.swaplab = prev_plot.swaplab
            else:
                self.fig.add_subplot(self.subv,self.subh,self.n+1,axisbg=pl.axisbg)
            if self.verbose: print "**** Done subplot %i / %i " %( self.n+1,self.howmanyplots ) 
            # finally make the plot
            pl.make()
            self.n = self.n+1 
        # once completed show the plot (cannot show intermediate plotin)
        # ... added a fix (customplot=True) for the label problem in basemap
        print "**** Done step: makeplot"
        if (self.n == self.howmanyplots): 
            ppplot.save(mode=self.out,filename=self.filename,folder=self.folder,custom=self.customplot)
            mpl.close()
        # SAVE A PICKLE FILE WITH THE self.p ARRAY OF OBJECTS
        if self.verbose: print "**** Saving session in "+self.filename + ".ppobj"
        savfile = self.folder + "/" + self.filename + ".ppobj"
        try: 
            filehandler = open(savfile, 'w')
            pickle.dump(self.p, filehandler)
        except IOError: 
            print "!! WARNING !! Saved object file not written. Probably do not have permission to write here."

    ###########################################################
    # plot: a shortcut method for the defineplot + plot chain #
    ###########################################################
    def plot(self,extraplot=0):
        self.defineplot(extraplot=extraplot)
        self.makeplot()

    #######################################################
    # getplot: a shortcut method for the get + plot chain #
    #######################################################
    def getplot(self,extraplot=0):
        self.get()
        self.plot(extraplot=extraplot)

    ###################################################################
    # getdefineplot: a shortcut method for the get + defineplot chain #
    ###################################################################
    def getdefineplot(self,extraplot=0):
        self.get()
        self.defineplot(extraplot=extraplot)

    ##############################################################
    # f: operation on two pp objects being on status 'definedplot'
    # this allows for one field being function of another one
    # e.g. u.f(v) means u will be displayed as a function of v
    # ... no need to do defineplot after u.f(v), makeplot directly
    ##############################################################
    def f(self,other):
        # preamble: for this operation to work, defineplot() must have been done
        if self.status != "definedplot":
            if self.verbose: print "!! WARNING !! performing defineplot on operand"
            self.defineplot()
        if other.status != "definedplot":
            if self.verbose: print "!! WARNING !! performing defineplot on operand"
            other.defineplot()
        # check total number of plots
        if self.howmanyplots != other.howmanyplots:
               print "!! ERROR !! The two operands do not have the same number of subplots."
               exit()
        # and now operation. 
        count = 0
        while count < self.howmanyplots:
           sobj = self.p[count] ; oobj = other.p[count]
           if sobj.field.ndim !=1 or oobj.field.ndim !=1:
               if self.verbose: print "!! WARNING !! Flattening arrays because more than one-dimensional."
               sobj.field = np.ravel(sobj.field)
               oobj.field = np.ravel(oobj.field)
           sobj.absc = oobj.field
           sobj.xlabel = oobj.ylabel
           if sobj.absc.size > sobj.field.size:
               if self.verbose:
                   print "!! WARNING !! Trying to define y=f(x) with x and y not at the same size.",sobj.absc.size,sobj.field.size
                   print "!! WARNING !! Modifying x to fit y size but please check." 
               sobj.absc = sobj.absc[0:sobj.field.size]
           count = count + 1
        return self

    ###########################################################
    # copyopt: get options from e.g. a parser
    # ... allow for simple scripting and user-defined settings
    # ... must be called between defineplot and makeplot
    # REQUIRED: attributes of opt must be the same as in the pp object
    ###########################################################
    def getopt(self,opt):
        # -- if only one, or less than the number of plots --> we take the first one
        # -- if as many as number of plots --> OK, each plot has its own setting
        # (except a few cases such as trans)
        for iii in range(self.howmanyplots):
            ###
            try: self.p[iii].trans = opt.trans
            except: pass
            ###
            try: self.p[iii].div = opt.div
            except: pass
            ###
            try: self.p[iii].logy = opt.logy
            except: pass
            ###
            try: self.p[iii].colorb = opt.colorb[iii]
            except: 
                try: self.p[iii].colorb = opt.colorb[0]
                except: pass
            ###
            try: self.p[iii].title = opt.title[iii]
            except: 
                try: self.p[iii].title = opt.title[0]
                except: pass
            ###
            try: self.p[iii].xlabel = opt.xlabel[iii]
            except: 
                try: self.p[iii].xlabel = opt.xlabel[0]
                except: pass
            ###
            try: self.p[iii].ylabel = opt.ylabel[iii]
            except: 
                try: self.p[iii].ylabel = opt.ylabel[0]
                except: pass
            ###
            try: self.p[iii].lstyle = opt.lstyle[iii]
            except: 
                try: self.p[iii].lstyle = opt.lstyle[0]
                except: pass
            ###
            try: self.p[iii].color = opt.color[iii]
            except:  
                try: self.p[iii].color = opt.color[0]
                except: pass
            ###
            try: self.p[iii].marker = opt.marker[iii]
            except:  
                try: self.p[iii].marker = opt.marker[0]
                except: pass
            ###
            try: self.p[iii].proj = opt.proj[iii]
            except: 
                try: self.p[iii].proj = opt.proj[0]
                except: pass
            ###
            try: self.p[iii].back = opt.back[iii]
            except: 
                try: self.p[iii].back = opt.back[0]
                except: pass
            ###
            try: self.p[iii].area = opt.area[iii]
            except: 
                try: self.p[iii].area = opt.area[0]
                except: pass
            ###
            try: self.p[iii].blon = opt.blon[iii]
            except: 
                try: self.p[iii].blon = opt.blon[0]
                except: pass
            ###
            try: self.p[iii].blat = opt.blat[iii]
            except: 
                try: self.p[iii].blat = opt.blat[0]
                except: pass
            ###
            try: self.p[iii].vmin = opt.vmin[iii]
            except: 
                try: self.p[iii].vmin = opt.vmin[0]
                except: pass
            ###
            try: self.p[iii].vmax = opt.vmax[iii]
            except: 
                try: self.p[iii].vmax = opt.vmax[0]
                except: pass

##########################################################
### THE ONEREQUEST SUBOBJECT TO PP (ON WHICH IT LOOPS) ###
##########################################################
class onerequest():

    # default settings. mostly initialized to diagnose problem, except dimplot, nplot, verbose, swap_axes, invert_axes
    # -------------------------------
    def __init__(self):
        self.file  = '!! file: I am not set, damned !!'
        self.f     = None
        self.dim   = None
        self.var   = '!! var: I am not set, damned !!'
        self.index_x = [] ; self.index_y = [] ; self.index_z = [] ; self.index_t = []
        self.index_x2d = [] ; self.index_y2d = []
        self.method_x = '!! method_x: I am not set, damned !!'
        self.method_y = '!! method_y: I am not set, damned !!'
        self.method_z = '!! method_z: I am not set, damned !!'
        self.method_t = '!! method_t: I am not set, damned !!'
        self.field = None
        self.name_x = None ; self.name_y = None ; self.name_z = None ; self.name_t = None
        self.dim_x = None ; self.dim_y = None ; self.dim_z = None ; self.dim_t = None
        self.field_x = None ; self.field_y = None ; self.field_z = None ; self.field_t = None
        self.dimplot = 0
        self.nplot = 1
        self.absc = None ; self.ordi = None ; self.absclab = None ; self.ordilab = None
        self.verbose = True
        self.swap_axes = False ; self.invert_axes = False
        self.compute = None

    # open a file. for now it is netcdf. TBD for other formats.
    # check that self.var is inside.
    # -------------------------------
    def openfile(self):
        if not os.path.exists(self.file): print '!! ERROR !! I could not find the following file: '+self.file ; exit()
        if not os.path.isfile(self.file): print '!! ERROR !! This does not appear to be a file: '+self.file ; exit()
        self.f = netCDF4.Dataset(self.file)
        if self.verbose: print "**** OK. Opened file "+self.file
        if self.var not in self.f.variables.keys(): 
            print '!! ERROR !! File '+self.file+' does not contain variable: '+self.var
            print '..... try instead with ',self.f.variables.keys() ; exit()

    # copy attributes from another existing object
    # --------------------------------------------
    def copy(self,source):
        for k, v in vars(source).items():
            setattr(self,k,v)

    # get x,y,z,t dimensions from NETCDF file
    # TBD: user could request for a specific altitude dimension
    # TBD: staggered variables could request specific dimensions
    # -------------------------------
    def getdim(self):
          # GET SIZES OF EACH DIMENSION
          if self.verbose: print "**** OK. Found variable "+self.var
          shape = self.f.variables[self.var].shape
          self.dim = len(shape)
          if self.dim == 1:
              if self.verbose: print "**** OK. 1D field. I assume this varies with time."
              self.dim_x = 1 ; self.dim_y = 1 ; self.dim_z = 1 ; self.dim_t = shape[0]
          elif self.dim == 2:
              if self.verbose: print "**** OK. 2D field. I assume this is not-time-varying lat-lon map."
              self.dim_x = shape[1] ; self.dim_y = shape[0] ; self.dim_z = 1 ; self.dim_t = 1
          elif self.dim == 3:
              if self.verbose: print "**** OK. 3D field. I assume this is time-varying lat-lon map."
              self.dim_x = shape[2] ; self.dim_y = shape[1] ; self.dim_z = 1 ; self.dim_t = shape[0]
          elif self.dim == 4:
              if self.verbose: print "**** OK. 4D field."
              self.dim_x = shape[3] ; self.dim_y = shape[2] ; self.dim_z = shape[1] ; self.dim_t = shape[0]
          # LONGITUDE. Try preset fields. If not present set grid points axis.
          self.name_x = "nothing"
          for c in glob_listx:
            if c in self.f.variables.keys():
             self.name_x = c
          if self.name_x == "nothing":
            self.field_x = np.array(range(self.dim_x))
            self.name_x = "x grid points"
          else:
            self.field_x = self.f.variables[self.name_x]
          # LATITUDE. Try preset fields. If not present set grid points axis.
          self.name_y = "nothing"
          for c in glob_listy:
            if c in self.f.variables.keys():
             self.name_y = c
          if self.name_y == "nothing":
            self.field_y = np.array(range(self.dim_y))
            self.name_y = "y grid points"
          else:
            self.field_y = self.f.variables[self.name_y]
          # ensure that lon and lat are 2D fields
          # 1. simple 1D case (not time-varying)
          if len(self.field_x.shape)*len(self.field_y.shape) == 1:
               if self.verbose: print "**** OK. recasting lon and lat as 2D fields."  
               [self.field_x,self.field_y] = np.meshgrid(self.field_x,self.field_y)
          # 2. complex 3D case (time-varying, actually just copied over time axis)
          elif len(self.field_x.shape)*len(self.field_y.shape) == 9:
               if self.verbose: print "**** OK. reducing lon and lat as 2D fields. get rid of time."
               self.field_x = self.field_x[0,:,:]
               self.field_y = self.field_y[0,:,:]
          # if xy axis are apparently undefined, set 2D grid points axis.
          if "grid points" not in self.name_x:
            if self.field_x.all() == self.field_x[0,0]:
               print "!! WARNING !! xy axis look undefined. creating a non-dummy ones."
               self.field_x = np.array(range(self.dim_x)) ; self.name_x = "x grid points"
               self.field_y = np.array(range(self.dim_y)) ; self.name_y = "y grid points"
               [self.field_x,self.field_y] = np.meshgrid(self.field_x,self.field_y)
          if self.dim_x > 1: 
               if self.verbose: print "**** OK. x axis %4.0f values [%5.1f,%5.1f]" % (self.dim_x,self.field_x.min(),self.field_x.max())
          if self.dim_y > 1: 
               if self.verbose: print "**** OK. y axis %4.0f values [%5.1f,%5.1f]" % (self.dim_y,self.field_y.min(),self.field_y.max())
          # ALTITUDE. Try preset fields. If not present set grid points axis.
          # WARNING: how do we do if several are available?
          self.name_z = "nothing"
          for c in glob_listz:
            if c in self.f.variables.keys():
             self.name_z = c
          if self.name_z == "nothing":
            self.field_z = np.array(range(self.dim_z))
            self.name_z = "z grid points"
          else:
            self.field_z = self.f.variables[self.name_z][:] # specify dimension
                                                            # TBD: have to check that this is not a 3D field
          if self.dim_z > 1: 
               if self.verbose: print "**** OK. z axis %4.0f values [%5.1f,%5.1f]" % (self.dim_z,self.field_z.min(),self.field_z.max())
          # TIME. Try preset fields.
          self.name_t = "nothing"
          for c in glob_listt:
            if c in self.f.dimensions.keys():
             self.name_t = c
          try:
            # speed up: only get first value, last one.
            dafirst = self.f.variables[self.name_t][0]
            dalast = self.f.variables[self.name_t][self.dim_t-1]
            self.field_t = np.linspace(dafirst,dalast,num=self.dim_t)
            if dafirst == dalast: self.field_t = np.array([dafirst])
          except:
            # ... or if a problem encountered, define a simple time axis
            self.field_t = np.array(range(self.dim_t))
            self.name_t = "t grid points"
          if self.dim_t > 1: 
               if self.verbose: print "**** OK. t axis %4.0f values [%5.1f,%5.1f]" % (self.dim_t,self.field_t.min(),self.field_t.max())     

    # get list of index to be retrieved for time axis
    ### TBD: il faudrait ne prendre que les indices qui correspondent a l interieur d un plot (dans all)
    # -------------------------------
    def getindextime(self,dalist=None,ind=None,stride=1):
        if self.method_t == "free": 
            self.index_t = np.arange(0,self.dim_t,stride)
            if self.dim_t > 1:  
                self.dimplot = self.dimplot + 1 
                if self.verbose: print "**** OK. t values. all."
            else:               
                self.method_t = "fixed"
                if self.verbose: print "**** OK. no t dimension."
        elif self.method_t == "comp":
            start = np.argmin( np.abs( self.field_t - dalist[ind][0] ) )
            stop = np.argmin( np.abs( self.field_t - dalist[ind][1] ) )
            self.index_t = np.arange(start,stop,stride)
            if self.verbose: print "**** OK. t values. comp over interval ",self.field_t[start],self.field_t[stop]," nvalues=",self.index_t.size
        elif self.method_t == "fixed":
            self.index_t.append( np.argmin( np.abs( self.field_t - dalist[ind][0] ) ))
            if self.verbose: print "**** OK. t values",self.field_t[self.index_t]
        else:
            print "!! ERROR !! method "+self.method_t+" not supported"
        self.index_t = np.array(self.index_t)
             
    # get list of index to be retrieved for vertical axis
    ### TBD: il faudrait ne prendre que les indices qui correspondent a l interieur d un plot (dans all)
    # -------------------------------
    def getindexvert(self,dalist=None,ind=None,stride=1):
        if self.method_z == "free": 
            self.index_z = np.arange(0,self.dim_z,stride)
            if self.dim_z > 1:  
                self.dimplot = self.dimplot + 1
                if self.verbose: print "**** OK. z values. all."
            else:               
                self.method_z = "fixed"
                if self.verbose: print "**** OK. no z dimension."
        elif self.method_z == "comp":
            start = np.argmin( np.abs( self.field_z - dalist[ind][0] ) )
            stop = np.argmin( np.abs( self.field_z - dalist[ind][1] ) )
            self.index_z = np.arange(start,stop,stride)
            if self.verbose: print "**** OK. z values. comp over interval",self.field_z[start],self.field_z[stop]," nvalues=",self.index_z.size
        elif self.method_z == "fixed":
            self.index_z.append( np.argmin( np.abs( self.field_z - dalist[ind][0] ) ))
            if self.verbose: print "**** OK. z values",self.field_z[self.index_z]
        else:
            if self.verbose: print "!! ERROR !! method "+self.method_z+" not supported"
        self.index_z = np.array(self.index_z)

    # get list of index to be retrieved for horizontal grid
    # --> index_x and index_y are slices to be retrieved from NETCDF files
    # --> index_x2d and index_y2d are the actual (x,y) coordinates corresponding to each relevant point
    # [this is slightly more complicated because 2D arrays for lat-lon projection possibly irregular]
    # NB: to append index we use lists (the most convenient) then we convert into a numpy.array
    ### TBD: il faudrait ne prendre que les indices qui correspondent a l interieur d un plot (dans all)
    # -------------------------------
    def getindexhori(self,dalistx=None,dalisty=None,indx=None,indy=None,stridex=1,stridey=1):
        ## get what is the method over x and y axis
        test = self.method_x+self.method_y
        ## CASE 0, EASY CASES: 
        ## - LAT IS FREE (we do here what must be done whatever LON case is)
        ## - LON IS FREE (we do here what must be done whatever LAT case is)
        ## - LAT IS COMP AND LON IS FREE
        ## - LON IS COMP AND LAT IS FREE
        if self.method_x == "free" or test in ["compfree","compcomp"]:
            self.index_x = range(0,self.dim_x,stridex)
            if self.dim_x > 1:  
                if self.method_x == "free": self.dimplot = self.dimplot + 1
                if self.verbose: print "**** OK. x values. all."
            else:               
                self.method_x = "fixed"
                if self.verbose: print "**** OK. no x dimension."
        if self.method_y == "free" or test in ["freecomp","compcomp"]:
            self.index_y = range(0,self.dim_y,stridey)
            if self.dim_y > 1:  
                if self.method_y == "free": self.dimplot = self.dimplot + 1
                if self.verbose: print "**** OK. y values. all."
            else:               
                self.method_y = "fixed"
                if self.verbose: print "**** OK. no y dimension."
        ## CASE 0 above, this is just for continuity.
        if self.method_x in ["free","comp"] and self.method_y in ["free","comp"]:
            self.index_x2d = self.index_x
            self.index_y2d = self.index_y
        ## AND NOW THE LITTLE BIT MORE COMPLICATED CASES
        ## CASE 1 LAT AND LON ARE FIXED
        elif test == "fixedfixed":
            idy,idx = np.unravel_index( np.argmin( ( self.field_x - dalistx[indx][0])**2 + (self.field_y - dalisty[indy][0])**2 ), self.field_x.shape ) 
                          #TBD: pb with staggered coord
            if idx not in self.index_x:  self.index_x.append(idx)
            if idy not in self.index_y:  self.index_y.append(idy)
            self.index_x2d.append(idx)
            self.index_y2d.append(idy)
        ## CASE 2 LON IS FIXED BUT NOT LAT
        elif test in ["fixedfree","fixedcomp"]:
            # find where are requested x values for each y on the free dimension
            # NB: this does not work for non-bijective cases e.g. polar stereographic
            for iy in range(self.dim_y):
              idx = np.argmin( np.abs( self.field_x[iy,:] - dalistx[indx][0] ) )
              # if comp is requested we select only indexes which yield values between requested min and max
              storeval = (self.method_y == "comp") and (self.field_y[iy,idx] > dalisty[indy][0]) and (self.field_y[iy,idx] < dalisty[indy][1])
              storeval = storeval or (self.method_y == "free")
              if storeval:
                  if idx not in self.index_x:  self.index_x.append(idx)
                  if iy not in self.index_y and self.method_y == "comp": self.index_y.append(iy)
                  if idx not in self.index_x2d or iy not in self.index_y2d:
                    self.index_x2d.append(idx)
                    self.index_y2d.append(iy)
        ## CASE 3 LAT IS FIXED BUT NOT LON
        elif test in ["freefixed","compfixed"]:
            # find where are requested y values for each x on the free dimension
            # NB: this does not work for non-bijective cases e.g. polar stereographic
            for ix in range(self.dim_x):
              idy = np.argmin( np.abs( self.field_y[:,ix] - dalisty[indy][0] ) )
              # if comp is requested we select only indexes which yield values between requested min and max
              storeval = (self.method_x == "comp") and (self.field_x[idy,ix] > dalistx[indx][0]) and (self.field_x[idy,ix] < dalistx[indx][1])
              storeval = storeval or (self.method_x == "free")
              if storeval:
                  if idy not in self.index_y:  self.index_y.append(idy)
                  if ix not in self.index_x and self.method_x == "comp": self.index_x.append(ix)
                  if ix not in self.index_x2d or idy not in self.index_y2d:
                    self.index_x2d.append(ix)
                    self.index_y2d.append(idy)
        ## check index tab
        if len(self.index_x) == 0 or len(self.index_y) == 0:
            print "!! ERROR !! no indices found. check prescribed latitudes or longitudes" ; exit()
        ## ensure the array is a numpy array for getfield to work
        self.index_x = np.array(self.index_x)
        self.index_y = np.array(self.index_y)
        self.index_x2d = np.array(self.index_x2d)
        self.index_y2d = np.array(self.index_y2d)
        ### print extrema
        printx = self.field_x[np.ix_(self.index_y2d, self.index_x2d)]
        printy = self.field_y[np.ix_(self.index_y2d, self.index_x2d)]
        if self.verbose: 
            print "**** OK. x values (min,max).", printx.min(),printx.max()
            print "**** OK. y values (min,max).", printy.min(),printy.max()

    # get the field from the NETCDF file and perform averages
    # -------------------------------
    def getfield(self):
        ## first tell what is to be done
        if self.dimplot > 2:                       print "**** !! ERROR !! "+str(self.dimplot)+"D plots not supported!" ; exit()
        elif self.dimplot == 0 and self.verbose:   print "**** OK. 0D value requested."
        elif self.dimplot == 1 and self.verbose:   print "**** OK. 1D plot requested."
        elif self.verbose:                         print "**** OK. 2D section requested."
        # well, now get field from netcdf file
        # part below is necessary otherwise there is an index error below
        if self.index_x.size == 1: self.index_x = self.index_x[0]
        if self.index_y.size == 1: self.index_y = self.index_y[0]
        if self.index_z.size == 1: self.index_z = self.index_z[0]
        if self.index_t.size == 1: self.index_t = self.index_t[0]
        # then retrieve what is requested by user
        # each self.dim case corresponds to tests in the beginning of getdim.
        time0 = timelib.time()
        if self.verbose: print "**** OK. I am getting values from files. Please wait."
        if self.dim == 1:  
            nt = self.index_t.size ; nz = 1 ; ny = 1 ; nx = 1
            self.field = self.f.variables[self.var][self.index_t]
        elif self.dim == 2:
            nt = 1 ; nz = 1 ; ny = self.index_y.size ; nx = self.index_x.size
            self.field = self.f.variables[self.var][self.index_y,self.index_x]
        elif self.dim == 3:
            nt = self.index_t.size ; nz = 1 ; ny = self.index_y.size ; nx = self.index_x.size
            self.field = self.f.variables[self.var][self.index_t,self.index_y,self.index_x]
            # this is far faster than retrieving each term with a loop
        elif self.dim == 4:
            nt = self.index_t.size ; nz = self.index_z.size ; ny = self.index_y.size ; nx = self.index_x.size
            self.field = self.f.variables[self.var][self.index_t,self.index_z,self.index_y,self.index_x]
        else:
            print "!! ERROR !! field would have more than four dimensions ?" ; exit()
        # NB: ... always 4D array but possibly with "size 1" dimensions 
        #     ... if one dimension is missing because 1D 2D or 3D requests, make it appear again
        self.field = np.reshape(self.field,(nt,nz,ny,nx))
        if self.verbose: print "**** OK. I got %7.1e values. This took me %6.4f seconds" % (nx*ny*nz*nt,timelib.time() - time0)
        if self.verbose: print "**** OK. I got var "+self.var+" with shape",self.field.shape
        # reduce coordinates to useful points
        # ... TBD: this should be ordered in the case of non-regular projections
        if self.method_x in ["free","comp"] and self.method_y in ["free","comp"]:
          # we need 2D coordinates (free) or we get broadcast problem (comp) so we use np.ix
          self.field_x = self.field_x[np.ix_(self.index_y2d, self.index_x2d)]
          self.field_y = self.field_y[np.ix_(self.index_y2d, self.index_x2d)]
        else:
          # we are OK with 1D coordinates
          self.field_x = self.field_x[self.index_y2d, self.index_x2d]
          self.field_y = self.field_y[self.index_y2d, self.index_x2d]
        self.field_z = self.field_z[self.index_z]
        self.field_t = self.field_t[self.index_t]
        # now have to obtain the new indexes which correspond to the extracted self.field
        # for it to work with unique index, ensure that any index_* is a numpy array
        if not isinstance(self.index_x, np.ndarray): self.index_x = np.array([self.index_x])
        if not isinstance(self.index_y, np.ndarray): self.index_y = np.array([self.index_y])
        if not isinstance(self.index_z, np.ndarray): self.index_z = np.array([self.index_z])
        if not isinstance(self.index_t, np.ndarray): self.index_t = np.array([self.index_t])
        for val in self.index_x: self.index_x2d[np.where(self.index_x2d == val)] = np.where(self.index_x == val)[0]
        for val in self.index_y: self.index_y2d[np.where(self.index_y2d == val)] = np.where(self.index_y == val)[0]
        for val in self.index_z: self.index_z  [np.where(self.index_z   == val)] = np.where(self.index_z == val)[0]
        for val in self.index_t: self.index_t  [np.where(self.index_t   == val)] = np.where(self.index_t == val)[0]
               ##### VERY EXPENSIVE
               ## recast self.field with 2D horizontal arrays because we might have extracted
               ## more than what is to be actually plot or computed, in particular for comps on 2D lat,lon coordinates
               #self.field = self.field[np.ix_(self.index_t,self.index_z,self.index_y2d,self.index_x2d)]
               #(nt,nz,ny,nx) = self.field.shape        
        # extract relevant horizontal points
        # TBD: is compfree OK with computing on irregular grid?
        test = self.method_x + self.method_y
        if test in ["fixedfixed","freefree"]:
          pass
        elif test in ["fixedfree","fixedcomp"] or test in ["freefixed","compfixed"]: 
          time0 = timelib.time()  
          # prepare the loop on all relevant horizontal points
          if self.method_x in ["comp","free"]:    
              nnn = self.index_x2d.shape[0] ; what_I_am_supposed_to_do = "keepx"
          elif self.method_y in ["comp","free"]:  
              nnn = self.index_y2d.shape[0] ; what_I_am_supposed_to_do = "keepy" 
          # LOOP to extract only useful values over horizontal dimensions
          # only take diagonal terms, do not loop on all self.index_x2d*self.index_y2d
          # ... this method is fast enough, perhaps there is a faster way though
          # ... (for sure the method with np.diag is much slower)
          for iii in range(nnn):
           ix = self.index_x2d[iii] ; iy = self.index_y2d[iii]
           for iz in self.index_z:
            for it in self.index_t:
              if what_I_am_supposed_to_do == "keepx":    self.field[it,iz,0,ix] = self.field[it,iz,iy,ix]
              elif what_I_am_supposed_to_do == "keepy":  self.field[it,iz,iy,0] = self.field[it,iz,iy,ix]
          if self.verbose: print "**** OK. I got to pick the right values for your request. This took me %6.4f seconds" % (timelib.time() - time0)
          # we only keep the one value that was modified on the dimension which is not free
          if what_I_am_supposed_to_do == "keepx":     self.field = self.field[:,:,0,:] ; ny = 1 ; self.field = np.reshape(self.field,(nt,nz,ny,nx))
          elif what_I_am_supposed_to_do == "keepy":   self.field = self.field[:,:,:,0] ; nx = 1 ; self.field = np.reshape(self.field,(nt,nz,ny,nx))
        # make a mask in case there are non-NaN missing values. (what about NaN missing values?)
        # ... this is important for computations below (see ppcompute)
        masked = np.ma.masked_where(np.abs(self.field) > 1e25,self.field)
        if masked.mask.any() == True:
             if self.verbose: print "!! WARNING !! Values over +-1e25 are considered missing values."
             self.field = masked
             self.field.set_fill_value([np.NaN])

    # perform computations
    # -------------------------------
    # available: mean, max, min, meanarea
    # TB: integrals? for derivatives, define a function self.dx()
    def computations(self):  
        nt,nz,ny,nx = self.field.shape
        # treat the case of mean on fields normalized with grid mesh area
        # ... this is done in the .area() method. 
        # after that self.field contains field*area/totarea
        if "area" in self.compute: self.area()
        # now ready to compute [TBD: we would like to have e.g. mean over x,y and min over t ??]
        if self.method_t == "comp":
            if self.verbose: print "**** OK. Computing over t axis."
            if "mean" in self.compute: self.field = ppcompute.mean(self.field,axis=0)
            elif self.compute == "min": self.field = ppcompute.min(self.field,axis=0)
            elif self.compute == "max": self.field = ppcompute.max(self.field,axis=0)
            else: print "!! ERROR !! operation not supported." ; exit()
            nt = 1 ; self.field = np.reshape(self.field,(nt,nz,ny,nx))
        if self.method_z == "comp": 
            if self.verbose: print "**** OK. Computing over z axis."
            if "mean" in self.compute: self.field = ppcompute.mean(self.field,axis=1)
            elif self.compute == "min": self.field = ppcompute.min(self.field,axis=1)
            elif self.compute == "max": self.field = ppcompute.max(self.field,axis=1)
            else: print "!! ERROR !! operation not supported." ; exit()
            nz = 1 ; self.field = np.reshape(self.field,(nt,nz,ny,nx))
        if self.method_y == "comp": 
            if self.verbose: print "**** OK. Computing over y axis."
            if self.compute == "mean": self.field = ppcompute.mean(self.field,axis=2)
            elif self.compute == "min": self.field = ppcompute.min(self.field,axis=2)
            elif self.compute == "max": self.field = ppcompute.max(self.field,axis=2)
            elif self.compute == "meanarea": self.field = ppcompute.sum(self.field,axis=2)
            else: print "!! ERROR !! operation not supported." ; exit()
            ny = 1 ; self.field = np.reshape(self.field,(nt,nz,ny,nx))
            if self.field_x.ndim == 2: self.field_x = self.field_x[0,:] # TBD: this is OK for regular grid but not for irregular
        if self.method_x == "comp":
            if self.verbose: print "**** OK. Computing over x axis."
            if self.compute == "mean": self.field = ppcompute.mean(self.field,axis=3)
            elif self.compute == "min": self.field = ppcompute.min(self.field,axis=3)
            elif self.compute == "max": self.field = ppcompute.max(self.field,axis=3)
            elif self.compute == "meanarea": self.field = ppcompute.sum(self.field,axis=3)
            else: print "!! ERROR !! operation not supported." ; exit()
            nx = 1 ; self.field = np.reshape(self.field,(nt,nz,ny,nx))
            if self.field_y.ndim == 2: self.field_y = self.field_y[:,0] # TBD: this is OK for regular grid but not for irregular
        # remove all dimensions with size 1 to prepare plot (and check the resulting dimension with dimplot)
        self.field = np.squeeze(self.field)
        if self.field.ndim != self.dimplot: 
            print "!! ERROR !! Problem: self.field is different than plot dimensions", self.field.ndim, self.dimplot ; exit()
        if self.verbose: 
            print "**** OK. Final shape for "+self.var+" after averaging and squeezing",self.field.shape
    
    # get areas for computations and ponderate field by those areas
    # -------------------------------------------------------------
    def area(self):
        if self.verbose: print "**** OK. Get area array for computations."
        # create a request object for area
        # ... and copy known attributes from self
        aire = onerequest()
        aire.copy(self)
        # get area field name
        aire.var = "nothing"
        for c in glob_listarea:
         if c in aire.f.variables.keys():
            aire.var = c
        # do not try to calculate areas automatically 
        if aire.var == "nothing":
            print "!! ERROR !! area variable not found... needs to be added in set_ppclass.txt?"
            exit()
        # define area request dimensions
        aire.getdim()
        # ensure this is a 2D horizontal request and define indexes
        #    ... areas are not supposed to vary with time and height
        aire.method_x = "free" ; aire.method_y = "free"
        aire.getindexhori() ; aire.dimplot = 2
        aire.method_z = "fixed" ; aire.field_z = np.array([0]) ; aire.index_z = np.array([0])
        aire.method_t = "fixed" ; aire.field_t = np.array([0]) ; aire.index_t = np.array([0])
        # read the 2D area array in netCDF file
        aire.getfield()
        # normalize by total area
        aire.field = np.squeeze(aire.field)
        totarea = ppcompute.sum(ppcompute.sum(aire.field,axis=1),axis=0)
        if self.verbose: print "**** OK. Total area is: ",totarea
        aire.field = aire.field / totarea
        # reduce with self horizontal indexes
        if "fixed" in self.method_x+self.method_y:
            aire.field = aire.field[self.index_y,self.index_x]
        # tile area array over self t and z axis so that area field could be multiplied with self.field
        aire.field = np.tile(aire.field,(self.index_t.size,self.index_z.size,1,1))
        if self.field.shape != aire.field.shape:
            print "!! ERROR !! Problem in area(). Check array shapes." ; exit()
        else:
            self.field = self.field*aire.field

    # define coordinates for plot
    # -------------------------------
    def definecoord(self):
        I_got_abs = False ; I_got_ord = False
        # here is the thing. time is usually taken as an abscissa so we start with time.
        if self.method_t ==  "free": 
            self.absc = self.field_t ; self.absclab = self.name_t 
            I_got_abs = True
        # then we usually have x as an abscissa.
        if self.method_x == "free":
            if I_got_abs:  
                self.ordi = self.field_x ; self.ordilab = self.name_x
                I_got_ord = True
            else:          
                self.absc = self.field_x ; self.absclab = self.name_x
                I_got_abs = True
        # ... or we have y
        if self.method_y == "free":
            if I_got_abs:   
                self.ordi = self.field_y ; self.ordilab = self.name_y
                I_got_ord = True
            else:          
                self.absc = self.field_y ; self.absclab = self.name_y
                I_got_abs = True
        # ... and we end with z because it is usually not an abscissa (profiles).
        if self.method_z == "free":
            if self.field_z[0] > self.field_z[1]:
                self.invert_axes = True # the axis will be turned upside-down
            if I_got_abs:  
                self.ordi = self.field_z ; self.ordilab = self.name_z
                I_got_ord = True
            else:
                self.absc = self.field_z ; self.absclab = self.name_z
                I_got_abs = True
                self.swap_axes = True # says that altitude is not supposed to remain as an abscissa
        if I_got_abs and self.verbose: print "**** OK. abscissa:",self.absclab, self.absc.shape
        if I_got_ord and self.verbose: print "**** OK. ordinate:",self.ordilab, self.ordi.shape

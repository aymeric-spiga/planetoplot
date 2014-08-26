#######################
##### PLANETOPLOT #####
#######################

### A. Spiga     -- LMD -- 06~09/2011 -- General building and mapping capabilities
### T. Navarro   -- LMD -- 10~11/2011 -- Improved use for GCM and added sections + 1Dplot capabilities 
### A. Colaitis  -- LMD --    11/2011 -- Mostly minor improvements and inter-plot operation capabilities + zrecast interpolation for gcm
### A. Spiga     -- LMD -- 11~12/2011 -- Extended multivar subplot capabilities + cosmetic changes + general cleaning and tests
### A. Colaitis  -- LMD --    12/2011 -- Added movie capability [mencoder must be installed]
### A. Spiga     -- LMD --    12/2011 -- Added HTML animated page capability + general tests of consistency [winds, etc...] + consistent generic movie loop
### J. Leconte   -- LMD --    02/2012 -- Added area weighted averaging. Compatibility with terrestrial gcm.
### A. Spiga     -- LMD --    03/2012 -- Cleaning and improved comments 
### T. Navarro   -- LMD --    04/2012 -- Added capabilities (e.g. histograms for difference maps)
### A. Colaitis  -- LMD -- 05~06/2012 -- Added capabilities for analysis of mesoscale files (wind speed, etc...)
### A. Spiga     -- LMD -- 04~07/2012 -- Added larger support of files + planets. Corrected a few bugs. Cleaning and improved comments
### A. Colaitis  -- LMD --    08/2012 -- Added functionalities for further analysis: spectra, hodographs, etc... plus improved flexibility (e.g. grid)

def planetoplot (namefiles,\
           level=0,\
           vertmode=0,\
           proj=None,\
           back=None,\
           target=None,
           stride=3,\
           var=None,\
           clb=None,\
           winds=False,\
           addchar=None,\
           vmin=None,\
           vmax=None,\
           tile=False,\
           zoom=None,\
           display=True,\
           hole=False,\
           save="gui",\
           anomaly=False,\
           var2=None,\
           ndiv=10,\
           mult=1.,\
           add=0.,\
           zetitle=["fill"],\
           slon=None,\
           slat=None,\
           svert=None,\
           stime=None,\
           outputname=None,\
           resolution=200,\
           ope=None,\
           fileref=None,\
           minop=0.,\
           maxop=0.,\
           titleref="fill",\
           invert_y=False,\
           xaxis=[None,None],\
           yaxis=[None,None],\
           ylog=False,\
           xlog=False,\
           yintegral=False,\
           blat=None,\
           blon=None,\
           tsat=False,\
           flagnolow=False,\
           mrate=None,\
           mquality=False,\
           trans=1,\
           zarea=None,\
           axtime=None,\
           redope=None,\
           seevar=False,\
           xlab=None,\
           ylab=None,\
           lbls=None,\
           lstyle=None,\
           cross=None,\
           markdevil=False,\
           facwind=1.,\
           trycol=False,\
           streamflag=False,\
           nocolorb=False,\
           analysis=None,\
           monster=False):

    ####################################################################################################################
    ### Colorbars http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps?action=AttachFile&do=get&target=colormaps3.png

    #################################
    ### Load librairies and functions
    import netCDF4

    import matplotlib as mpl
    import matplotlib.pyplot
    import matplotlib.cm
    from mpl_toolkits.basemap import cm
    if mrate is not None: from videosink import VideoSink

    import numpy as np
    from numpy.core.defchararray import find
    import scipy
    if analysis in ['laplace']: import scipy.ndimage.laplace as laplace

    import itertools
    import os

    import myplot
    from mymath import deg,max,min,mean,writeascii,fig2data,fig2img
    from times import sol2ls
    #from singlet import singlet


#########################
### PRELIMINARY STUFF ###
#########################
### we say hello
    print "********************************************"
    print "********** WELCOME TO PLANETOPLOT **********"
    print "********************************************"
    if monster: 
        print "********* SPEED MODE. EXPERIMENTAL. ********"
        print "********************************************"
### we ensure namefiles and var are arrays and we get their size
    if not isinstance(namefiles, np.ndarray): namefiles = [namefiles]
    if not isinstance(var, np.ndarray):       var = [var]
    zelenfile = len(namefiles) ; zelenvar = len(var)
### we check colorbars options (including trycol i.e. the user wants to try a set of colorbars)
    if trycol: clb = ["Greys","Blues","YlOrRd","jet","spectral","hot","RdBu","RdYlBu","Paired"] ; zetitle = clb ; var = [var[0]]*9 ; zelenvar = len(var)
    if clb is None:            clb = ["def"]*zelenvar
    elif len(clb) < zelenvar:  clb = [clb[0]]*zelenvar ; print "WARNING: less color than vars! setting all to 1st value."
### we initialize a few variables
    initime=-1 ; sslon = None ; sslat = None
    k = 0 ; firstfile = True ; count = 0
### we perform sanity checks and correct for insufficient information
    if slon is not None: sslon = np.zeros([1,2])
    if slat is not None: sslat = np.zeros([1,2])
    if redope is not None and winds: winds=False ; print "WARNING: no winds with redope. setting winds to False."
### we avoid specific cases not yet implemented
    if mrate is not None and zelenvar > 1: myplot.errormess("multivar not allowed in movies. should be fixed soon!")
### we prepare number of plot fields [zelen] and number of plot instances [numplot] according to user choices
### --> we support multifile and multivar plots : files + vars separated by commas are on the same figure
    nlon, nlat, nvert, ntime, mapmode, nslices = myplot.determineplot(slon, slat, svert, stime, redope)
    zelen = zelenfile*zelenvar
    if (nslices > 1 and monster): myplot.errormess("multislice + monster not supported yet. to be done soon")
### we have a special mode obtained by -p noproj in which lat/lon plots are just flat 2D plots
    if proj == "noproj": mapmode = 0
### we correct number of plot fields for possible operation (substract, etc...)
    if ope is not None:
        if fileref is not None:       zelen = 3*zelenvar*zelenfile
        elif "var" in ope:            zelen = zelen + 1
    numplot = zelen*nslices
    print "********** FILES, SLICES, VARS, TOTAL PLOTS: ", zelenfile, nslices, zelenvar, numplot
    print "********** MAPMODE: ", mapmode
### we define the arrays for plot fields -- which will be placed within multiplot loops
    all_var  = [[]]*zelen ; all_var2  = [[]]*zelen
    all_title = [[]]*zelen ; all_varname = [[]]*zelen ; all_namefile = [[]]*zelen ; all_colorb = [[]]*zelen
    all_time = [[]]*zelen ; all_vert = [[]]*zelen ; all_lat = [[]]*zelen ; all_lon = [[]]*zelen
    all_windu = [[]]*zelen ; all_windv = [[]]*zelen
    plot_x = [[]]*zelen ; plot_y = [[]]*zelen ; multiplot = [[]]*zelen
### tool (should be move to mymath)
    getVar = lambda searchList, ind: [searchList[i] for i in ind]
#############################
### LOOP OVER PLOT FIELDS ###
#############################
   
    for nnn in range(zelenfile):
     for vvv in range(zelenvar): 

    ### we load each NETCDF objects in namefiles
      namefile = namefiles[nnn] 
      nc  = netCDF4.Dataset(namefile)
    ### we explore the variables in the file
      varinfile = nc.variables.keys()
      if seevar: print varinfile ; exit() 
    ### we define the type of file we have (gcm, meso, etc...)
      typefile = myplot.whatkindfile(nc)
      if firstfile:                 typefile0 = typefile
      elif typefile != typefile0:   myplot.errormess("Not the same kind of files !", [typefile0, typefile])
    ### we care for input file being 1D simulations
      is1d=999 
      if "longitude" in nc.dimensions and "latitude" in nc.dimensions: is1d = len(nc.dimensions["longitude"])*len(nc.dimensions["latitude"])
      elif "lon" in nc.dimensions and "lat" in nc.dimensions: is1d = len(nc.dimensions["lon"])*len(nc.dimensions["lat"])
      if typefile in ['gcm','earthgcm'] and is1d == 1:       mapmode=0 ; winds=False
    ### we create default vert and time prescriptions if not here in case mapping mode is on (lat/lon maps)
      if redope is None and mapmode == 1:
          if svert is None:  svert = myplot.readslices(str(level)) ; nvert=1
          if stime is None and mrate is None:
             stime = myplot.readslices(str(0)) ; ntime=1 ## this is a default choice
             print "WELL... nothing about time axis. I took default: first time reference stored in file."
    ### we get the names of variables to be read. in case only one available, we choose this one.
    ### (we take out of the test specific names e.g. UV is not in the file but used to ask a wind speed computation)
      varname = myplot.select_getfield(zvarname=var[vvv],znc=nc,ztypefile=typefile,mode='check')
    ### we get the names of wind variables to be read (if any)
      if winds:                                                    
         [uchar,vchar,metwind] = myplot.getwinddef(nc)             
         if uchar == 'not found': winds = False
    ### we tell the user that either no var or no wind is not acceptable
      if not varname and not winds: myplot.errormess("please set at least winds or var",printvar=nc.variables)
    ### we get the coordinates lat/lon to be used
      [lon2d,lat2d] = myplot.getcoorddef(nc)
    ### we get an adapted map projection if none is provided by the user
      if proj == None:   proj = myplot.getproj(nc)   
    ### we define plot boundaries given projection or user choices
      if firstfile:
         if proj in ["npstere","spstere"]: [wlon,wlat] = myplot.polarinterv(lon2d,lat2d)
         elif proj in ["lcc","laea"]:      [wlon,wlat] = myplot.wrfinterv(lon2d,lat2d)
         else:                             [wlon,wlat] = myplot.simplinterv(lon2d,lat2d)
         if zoom:                          [wlon,wlat] = myplot.zoomset(wlon,wlat,zoom)
         elif zarea is not None:           [wlon,wlat] = myplot.latinterv(area=zarea)

#############################################################
############ WE LOAD 4D DIMENSIONS : x, y, z, t #############
#############################################################

      ###should be available for GCM or simple files ?
      ##if ope in ["cat"] and nnn > 0:    count = time[-1] + 1  ## so that a cat is possible with simple subscripts

    ### TYPE 1 : GCM files or simple files
      if typefile in ["gcm","earthgcm","ecmwf"]:
      ### this is needed for continuity 
          if slon is not None: sslon = slon  
          if slat is not None: sslat = slat  
      ### we define lat/lon vectors. we get what was done in getcoorddef.
          lat = lat2d[:,0] ; lon = lon2d[0,:]
      ### we define areas. this is needed for calculate means and weight with area. this is not compulsory (see reduce_field).
          if "aire" in nc.variables:      area = nc.variables["aire"][:,:]
          ### --> add a line here if your reference is not present
          else:                           area = None
      ### we define altitude vector. either it is referenced or it is guessed based on last variable's dimensions.
          if "altitude" in nc.variables:   vert = nc.variables["altitude"][:]
          elif "Alt" in nc.variables:      vert = nc.variables["Alt"][:]
          elif "lev" in nc.variables:      vert = nc.variables["lev"][:]
          elif "presnivs" in nc.variables: vert = nc.variables["presnivs"][:]
          ### --> add a line here if your reference is not present
          else: 
              dadim = myplot.getdimfromvar(nc) ; print "No altitude found. Try to build a simple axis.",dadim
              if   len(dadim) == 4:  print "-- 4D field. Assume z is dim 2." ; vert = np.arange(dadim[-3])
              elif len(dadim) == 3:  print "-- 3D field. Assume z is dim 1." ; vert = [0.]
              else:                  vert = [0.]
      ### we define time vector. either it is referenced or it is guessed based on last variable's dimensions.
          if "Time" in nc.variables:            letime = "Time"
          elif "time_counter" in nc.variables:  letime = "time_counter"
          elif "time" in nc.variables:          letime = "time"
          ### --> add a line here if your reference is not present
          else:                                 letime = None
          if letime is not None:          
              if monster:
                  ### very long simulation... we just retrieve 3 values for time
                  timelength = len(nc.dimensions[letime])
                  dafirst = nc.variables[letime][0]
                  dalast = nc.variables[letime][timelength-1]
                  step = nc.variables[letime][1] - dafirst
                  time = np.arange(dafirst,dalast,step)
                  print "Time (start,stop,nvalues)",dafirst,dalast,timelength
              else:
                  ### if the time record is not too long, what follows is pretty quick
                  time = nc.variables[letime][:]
          else: 
              print "No time found. Try to build a simple axis. Assume t is dim 1." 
              dadim = myplot.getdimfromvar(nc)
              if   len(dadim) == 4:  time = np.arange(dadim[-4])
              elif len(dadim) == 3:  time = np.arange(dadim[-3])
              else:                  time = [0.] #myplot.errormess("no time axis found.")
          ### (SPECIFIC?)
          if "time_counter" in nc.variables:  time = time / 86400. #### convert from s to days
          ### (SPECIFIC. convert to Ls for Martian GCM files.)
          if axtime in ["ls"]:
              print "converting to Ls ..."
              for iii in range(len(time)):
                time[iii] = sol2ls(time[iii])
                if iii > 0:
                  while abs(time[iii]-time[iii-1]) > 300: time[iii] = time[iii]+360
          ### (a case where the user would like to set local times. e.g. : 1D plot with no longitude reference)
          if axtime in ["lt"]:
              if initime == -1: initime=input("Please type initial local time:")
              time = (initime+time*24)%24 ; print "LOCAL TIMES.... ", time
          ### (simply ask for subscript)
          if axtime in ["ind"]:
              dadim = myplot.getdimfromvar(nc)
              if   len(dadim) == 4:  time = np.arange(dadim[-4])
              elif len(dadim) == 3:  time = np.arange(dadim[-3])
              elif len(dadim) == 2:  time = np.arange(dadim[-2]) 

    ### TYPE 2 : MESOSCALE FILES
      elif typefile in ['meso','geo']:
          ### area not active with mesoscale files
          area = None 
          ### HACK TO GET THE CORRECT LAT/LON FROM MESOSCALE FILES WITH 2D LAT/LON ARRAYS
          ### principle: calculate correct indices then repopulate slon and slat
          if slon is not None or slat is not None:
              show_topo_map_user = ((save == 'png') and (typefile == 'meso') and ("HGT" in varinfile) and display)
              #show_topo_map_user = False
              if firstfile and show_topo_map_user:  iwantawhereplot = nc     #show a topo map with a cross on the chosen point
              else:                                 iwantawhereplot = None   #do not show anything, just select indices
              numlon = 1 ; numlat = 1 
              if slon is not None:   numlon = slon.shape[1]    
              if slat is not None:   numlat = slat.shape[1]
              indices = np.ones([numlon,numlat,2]) ; vlon = None ; vlat = None
              for iii in range(numlon):  
               for jjj in range(numlat):
                 if slon is not None:  vlon = slon[0][iii]  ### note: slon[:][0] does not work
                 if slat is not None:  vlat = slat[0][jjj]  ### note: slon[:][0] does not work
                 indices[iii,jjj,:] = myplot.bidimfind(lon2d,lat2d,vlon,vlat,file=iwantawhereplot)  
                 lonp,latp = ( lon2d[indices[iii,jjj,0],indices[iii,jjj,1]] , lat2d[indices[iii,jjj,0],indices[iii,jjj,1]] )
              ### possible bug here if several --lat
              for iii in range(numlon):
               for jjj in range(numlat):
                 if slon is not None: sslon[0][iii] = indices[iii,0,1] #...this is idx
                 if slat is not None: sslat[0][jjj] = indices[0,jjj,0] #...this is idy
              lonp,latp = ( lon2d[indices[0,0,0],indices[0,0,1]] , lat2d[indices[0,0,0],indices[0,0,1]] )
          ### we get rid of boundary relaxation zone for plots. important to do that now and not before.
          if (typefile in ['meso'] and mapmode == 1):
             if '9999' not in getattr(nc,'START_DATE'): lon2d = myplot.dumpbdy(lon2d,6) ; lat2d = myplot.dumpbdy(lat2d,6)  
          ### we read the keyword for vertical dimension. we take care for vertical staggering.
          if varname in ['PHTOT','W']:    vertdim='BOTTOM-TOP_PATCH_END_STAG'
          else:                           vertdim='BOTTOM-TOP_PATCH_END_UNSTAG'
          if (var2 is not None and var2 not in ['PHTOT','W']): dumped_vert_stag=True ; vertdim='BOTTOM-TOP_PATCH_END_UNSTAG'
          else:                                                dumped_vert_stag=False
          ### we read the keyword for horizontal dimensions. we take care for horizontal staggering.
          if varname in ['V']:  latdim='SOUTH-NORTH_PATCH_END_STAG'
          else:                 latdim='SOUTH-NORTH_PATCH_END_UNSTAG'
          if varname in ['U']:  londim='WEST-EAST_PATCH_END_STAG'
          else:                 londim='WEST-EAST_PATCH_END_UNSTAG'
          lon = np.arange(0,getattr(nc,londim),1) ; lat = np.arange(0,getattr(nc,latdim),1)
          ### we define the time axis and take care of various specificities (lt, ls, sol) or issues (concatenation)
          if axtime in ["ls","sol"]:
              lstab, soltab, lttab = myplot.getlschar ( namefile, getaxis = True )
              if axtime == "ls":      time = lstab 
              elif axtime == "sol":   time = soltab
          else:
              if ope in ["cat"] and nnn > 0: 	count = time[-1] + 1  ## so that a cat is possible with simple subscripts
              else: 				count = 0
              if "Times" in nc.dimensions:   time = count + np.arange(0,len(nc.dimensions["Times"]),1)
              elif "Time" in nc.dimensions:  time = count + np.arange(0,len(nc.dimensions["Time"]),1)
              else:                          time = count + np.arange(0,1,1)
          if axtime in ["lt"]:
              ftime = np.zeros(len(time))
              for i in range(len(time)): ftime[i] = myplot.localtime ( time[i], slon , namefile )
              time=ftime ; time = myplot.check_localtime(time)
              print "LOCAL TIMES.... ", time
          ### we define the vertical axis (or lack thereof) and cover possibilities for it to be altitude, pressure, geopotential. quite SPECIFIC.
          if typefile in ['geo']:   vert = [0.] ; stime = myplot.readslices(str(0))
          else:
              if vertmode is None:  vertmode=0
              if vertmode == 0:     
                  if "vert" in nc.variables: vert = nc.variables["vert"][:]/1000. ; vertmode = 1
                  else:                      vert = np.arange(0,getattr(nc,vertdim),1)
              elif vertmode == -1:  vert = nc.variables["PHTOT"][0,:,0,0]/3.72 ; vert = np.array(vert[0:len(vert)-1]) #; print vert
              elif vertmode == 1 or vertmode == 2:  vert = nc.variables["vert"][:]        ## pressure in Pa
              else:                                 vert = nc.variables["vert"][:]/1000.  ## altitude in km
####################################################################
############ END of WE LOAD 4D DIMENSIONS : x, y, z, t #############
####################################################################


    ### we fill the arrays of varname, namefile, time, colorbar at the current step considered (NB: why use both k and nnn ?)
      all_varname[k] = varname
      all_namefile[k] = namefile
      all_colorb[k] = clb[vvv]

##############################################################################
##### LARGE FILE. WE'D BETTER NOT FILL THE MEMORY WITH THE STUPID WHOLE ARRAY!
      if monster:
        ####################################################################
        ## get all indexes to be taken into account for this subplot and then reduce field
        ## We plot 1) all lon slices 2) all lat slices 3) all vert slices 4) all time slices and then go to the next slice
        if ope is not None:
            if fileref is not None:      index_f = (k//(nlon*nlat*nvert*ntime))%(3*zelenfile)  ## OK only 1 var,  see test in the beginning
            elif "var" in ope:           index_f = (k//(nlon*nlat*nvert*ntime))%(zelenvar+1)        ## OK only 1 file, see test in the beginning
            elif "cat" in ope:           index_f = 0
        elif not (k == 0):
           #if zelenfile > 1 and zelenvar == 1 and which == "unidim": pass ## TROUVER UN MOYEN POUR unidim
           if zelenfile > 1 and zelenvar == 1: pass
           else: yeah = zelenfile*zelenvar ; index_f = (k//(nlon*nlat*nvert*ntime))%yeah
        else: yeah = zelenfile*zelenvar ; index_f = (k//(nlon*nlat*nvert*ntime))%yeah

        ilon  = myplot.getsindex(sslon,k%nlon,lon)
        ilat  = myplot.getsindex(sslat,(k//nlon)%nlat,lat)
        ivert = myplot.getsindex(svert,(k//(nlon*nlat))%nvert,vert)
 
        if mrate is not None:                 itime = None
        else:                                 itime = myplot.getsindex(stime,(k//(nlon*nlat*nvert))%ntime,time)
        ltst = None
        if typefile in ['meso'] and itime is not None and len(itime) < 2: ltst = myplot.localtime ( itime, slon , all_namefile[index_f] )
 
        if ilat  is not None:  print "********** LAT : INDEX",ilat[0],  ilat[-1],  "VALUE",lat[ilat[0]],   lat[ilat[-1]]
        else:                  ilat = range(len(lat))
        if ilon  is not None:  print "********** LON : INDEX",ilon[0],  ilon[-1],  "VALUE",lon[ilon[0]],   lon[ilon[-1]]
        else:                  ilon = range(len(lon))
        if ivert is not None:  print "********** VERT: INDEX",ivert[0], ivert[-1], "VALUE",vert[ivert[0]], vert[ivert[-1]]
        else:                  ivert = range(len(vert))
        if itime is not None:  print "********** TIME: INDEX",itime[0], itime[-1], "VALUE",time[itime[0]], time[itime[-1]]
        else:                  itime = range(len(time))

        all_var[k] = myplot.getfieldred(nc,all_varname[k],ilon,ilat,ivert,itime)
        if var2:   all_var2[k] = myplot.getfieldred(nc,var2,ilon,ilat,ivert,itime)
        if winds:  
            all_windu[k] = myplot.getfield(nc,uchar,ilon,ilat,ivert,itime) 
            all_windv[k] = myplot.getfield(nc,vchar,ilon,ilat,ivert,itime)
        plot_x[k] = None ; plot_y[k] = None
        
        all_time[k] = time[itime]
        all_vert[k] = vert[ivert]
        all_lat[k]  = lat[ilat]
        all_lon[k]  = lon[ilon]

      else:
        ## regular stuff: not large array.
        all_time[k] = time
        all_vert[k] = vert
        all_lat[k]  = lat
        all_lon[k]  = lon
        if var2:  all_var2[k], plot_x[k], plot_y[k] = myplot.select_getfield(zvarname=var2,znc=nc,ztypefile=typefile,\
                         mode='getvar',ztsat=tsat,ylon=all_lon[k],ylat=all_lat[k],yalt=all_vert[k],ytime=all_time[k],analysis=analysis)
        if winds: all_windu[k] = myplot.getfield(nc,uchar) ; all_windv[k] = myplot.getfield(nc,vchar)
        ### we fill the arrays of fields to be plotted at the current step considered
        all_var[k], plot_x[k], plot_y[k] = myplot.select_getfield(zvarname=all_varname[k],znc=nc,ztypefile=typefile,\
                         mode='getvar',ztsat=tsat,ylon=all_lon[k],ylat=all_lat[k],yalt=all_vert[k],ytime=all_time[k],analysis=analysis)
####################################################################

      # last line of for namefile in namefiles
      print "**** GOT SUBDATA:",k," NAMEFILE:",namefile," VAR:",varname, var2 ; k += 1 ; firstfile = False

    ### note to self : stopped comment rewriting here.

    ##################################
    ### Operation on files (I) with _var
    if ope is not None and "var" in ope:
         print "********** OPERATION: ",ope
         if zelenfile > 1: myplot.errormess("for this operation... please set only one file !")
         if zelenvar > 2:       myplot.errormess("not sure this works for more than 2 vars... please check.")
         if   "div_var" in ope: all_var[k] = all_var[k-2] / all_var[k-1] ; insert = '_div_'
         elif "mul_var" in ope: all_var[k] = all_var[k-2] * all_var[k-1] ; insert = '_mul_'
         elif "add_var" in ope: all_var[k] = all_var[k-2] + all_var[k-1] ; insert = '_add_'
         elif "sub_var" in ope: all_var[k] = all_var[k-2] - all_var[k-1] ; insert = '_sub_'
         else:                    myplot.errormess(ope+" : non-implemented operation. Check pp.py --help")
         plot_x[k] = None ; plot_y[k] = None
         all_time[k] = all_time[k-1] ; all_vert[k] = all_vert[k-1] ; all_lat[k] = all_lat[k-1] ; all_lon[k] = all_lon[k-1] ; all_namefile[k] = all_namefile[k-1]
         all_varname[k] = all_varname[k-2] + insert + all_varname[k-1]
         if len(clb) >= zelen: all_colorb[k] = clb[-1]   # last additional user-defined color is for operation plot 
         else:                 all_colorb[k] = all_colorb[k-1]  # if no additional user-defined color... set same as var
         ### only the operation plot. do not mention colorb so that it is user-defined?
         if "only" in ope:
             numplot = 1 ; all_var[0] = all_var[k]
             all_time[0] = all_time[k] ; all_vert[0] = all_vert[k] ; all_lat[0] = all_lat[k] ; all_lon[0] = all_lon[k] ; all_namefile[0] = all_namefile[k]
             all_varname[0] = all_varname[k-2] + insert + all_varname[k-1]

    if ope is not None and "var" not in ope:
       print "********** OPERATION: ",ope
       if zelenvar > 1: myplot.errormess("for this operation... please set only one var !")
       if fileref is None and ope not in ['cat']: myplot.errormess("fileref is missing!")
       ################################## this is terribly complex for the main program file. got to do something about this
       ### Operation on files (II) without _var
       # we re-iterate on the plots to set operation subplots to make it compatible with multifile (using the same ref file)
       # (k+1)%3==0 is the index of operation plots
       # (k+2)%3==0 is the index of reference plots
       # (k+3)%3==0 is the index of first plots
       opefirstpass=True
       for k in np.arange(zelen):
               if ope in ["-","+","-%","-_only","+_only","-%_only","-_histo"]:
                  ncref = netCDF4.Dataset(fileref)

                  if opefirstpass: ## first plots
                     for ll in np.arange(zelenfile):
                        print "SETTING FIRST PLOT"
                        myplot.alltransfer(3*ll,ll,all_varname,all_time,all_vert,all_lat,all_lon,all_namefile,all_var2,all_colorb,all_var)
                        if plot_y[ll] is not None: plot_y[3*ll] = plot_y[ll] ; plot_x[3*ll] = plot_x[ll]
                        else: plot_y[3*ll] = None ; plot_x[3*ll] = None
                        opefirstpass=False

                  if (k+2)%3==0: ## reference plots
                        print "SETTING REFERENCE PLOT"
                        myplot.alltransfer(k,k-1,all_varname,all_time,all_vert,all_lat,all_lon,all_namefile,all_var2,all_colorb,all_var)
                        ### all_var is actually overwritten here
                        all_var[k], plot_x[k], plot_y[k] = myplot.select_getfield(zvarname=all_varname[k-1],znc=ncref,ztypefile=typefile,mode='getvar',ztsat=tsat,ylon=all_lon[k],ylat=all_lat[k],yalt=all_vert[k],ytime=all_time[k],analysis=analysis)
                        if winds: all_windu[k] = myplot.getfield(ncref,uchar) ; all_windv[k] = myplot.getfield(ncref,vchar)

                  if (k+1)%3==0: ## operation plots
                     print "SETTING OPERATION PLOT"
                     all_varname[k] = all_varname[k-1] ; all_time[k] = all_time[k-1] ; all_vert[k] = all_vert[k-1] ; all_lat[k] = all_lat[k-1] ; all_lon[k] = all_lon[k-1] ; all_namefile[k] = all_namefile[k-1] ; all_var2[k] = all_var2[k-1]
                     if ope in ["-","-_only","-_histo"]:
                         all_var[k]= all_var[k-2] - all_var[k-1]
                         if plot_y[k-1] is not None and plot_y[k-2] is not None: plot_y[k] = plot_y[k-2] - plot_y[k-1]
                         if plot_y[k-2] is None: plot_y[k] = None; plot_x[k] = None
                     elif ope in ["+","+_only"]:   
                         all_var[k]= all_var[k-2] + all_var[k-1]
                         if plot_y[k-1] is not None and plot_y[k-2] is not None: plot_y[k] = plot_y[k-2] + plot_y[k-1]
                         if plot_y[k-2] is None: plot_y[k] = None; plot_x[k] = None
                     elif ope in ["-%","-%_only"]:
                         masked = np.ma.masked_where(all_var[k-1] == 0,all_var[k-1])
                         masked.set_fill_value([np.NaN])
                         all_var[k]= 100.*(all_var[k-2] - masked)/masked
                         if plot_y[k-1] is not None and plot_y[k-2] is not None: 
                            masked = np.ma.masked_where(plot_y[k-1] == 0,plot_y[k-1])
                            masked.set_fill_value([np.NaN])
                            plot_y[k]= 100.*(plot_y[k-2] - masked)/masked
                         if plot_y[k-2] is None: plot_y[k] = None; plot_x[k] = None
                     if len(clb) >= zelen: all_colorb[k] = clb[-1]
                     else: all_colorb[k] = "RdBu_r" # if no additional user-defined color... set a good default one
                     if winds: all_windu[k] = all_windu[k-2]-all_windu[k-1] ; all_windv[k] = all_windv[k-2] - all_windv[k-1]

               elif ope in ["cat"]:
                  tabtime = all_time[0];tab = all_var[0];k = 1
                  if var2: tab2 = all_var2[0]
                  while k != zelenfile and len(all_time[k]) != 0:
                      if var2: tab2 = np.append(tab2,all_var2[k],axis=0) 
                      tabtime = np.append(tabtime,all_time[k]) ; tab = np.append(tab,all_var[k],axis=0) ; k += 1
                  all_time[0] = np.array(tabtime) ; all_var[0] = np.array(tab) ; numplot = 1
                  if var2: all_var2[0] = np.array(tab2)
               else: myplot.errormess(ope+" : non-implemented operation. Check pp.py --help")
       if "only" in ope:
           numplot = 1 ; all_var[0] = all_var[k]
           all_time[0] = all_time[k] ; all_vert[0] = all_vert[k] ; all_lat[0] = all_lat[k] ; all_lon[0] = all_lon[k] ; all_namefile[0] = all_namefile[k] ; plot_x[0]=plot_x[k] ; plot_y[0]=plot_y[k]
           all_varname[0] = all_varname[k]

    ##################################
    ### Open a figure and set subplots
    fig = mpl.pyplot.figure()
    subv,subh = myplot.definesubplot( numplot, fig, ipreferline = (mapmode == 1) ) 
    if ope in ['-','-%','-_histo'] and zelenfile ==1 : subv,subh = 2,2
    elif ope in ['-','-%'] and zelenfile>1 : subv, subh = zelenfile, 3
 
    #################################
    ### Time loop for plotting device
    nplot = 1 ; error = False ; firstpass = True 
    if lstyle is not None: linecycler = itertools.cycle(lstyle)
    else: linecycler = itertools.cycle(["-","--","-.",":"])
    print "********************************************"
    while error is False:
     
       print "********** PLOT", nplot, " OF ",numplot
       if nplot > numplot: break

       ####################################################################
       ## get all indexes to be taken into account for this subplot and then reduce field
       ## We plot 1) all lon slices 2) all lat slices 3) all vert slices 4) all time slices and then go to the next slice
       if ope is not None:
           if fileref is not None:      index_f = ((nplot-1)//(nlon*nlat*nvert*ntime))%(3*zelenfile)  ## OK only 1 var,  see test in the beginning
           elif "var" in ope:           index_f = ((nplot-1)//(nlon*nlat*nvert*ntime))%(zelenvar+1)        ## OK only 1 file, see test in the beginning
           elif "cat" in ope:           index_f = 0
       elif not firstpass:
          if zelenfile > 1 and zelenvar == 1 and which == "unidim": pass
          else: yeah = zelenfile*zelenvar ; index_f = ((nplot-1)//(nlon*nlat*nvert*ntime))%yeah
       else: yeah = zelenfile*zelenvar ; index_f = ((nplot-1)//(nlon*nlat*nvert*ntime))%yeah
       time = all_time[index_f] ; vert = all_vert[index_f] ; lat = all_lat[index_f] ; lon = all_lon[index_f]
       indexlon  = myplot.getsindex(sslon,(nplot-1)%nlon,lon)
       indexlat  = myplot.getsindex(sslat,((nplot-1)//nlon)%nlat,lat)
       indexvert = myplot.getsindex(svert,((nplot-1)//(nlon*nlat))%nvert,vert)
       plotx=plot_x[index_f] ; ploty=plot_y[index_f]
       if mrate is not None:                 indextime = None 
       else:                                 indextime = myplot.getsindex(stime,((nplot-1)//(nlon*nlat*nvert))%ntime,time)
       ltst = None 
       if typefile in ['meso'] and indextime is not None and len(indextime) < 2: ltst = myplot.localtime ( indextime, slon , all_namefile[index_f] )

       if not monster:
         if indexlat  is not None: print "********** LAT : INDEX",indexlat[0],  indexlat[-1],  "VALUE",lat[indexlat[0]],   lat[indexlat[-1]]
         if indexlon  is not None: print "********** LON : INDEX",indexlon[0],  indexlon[-1],  "VALUE",lon[indexlon[0]],   lon[indexlon[-1]]
         if indexvert is not None: print "********** VERT: INDEX",indexvert[0], indexvert[-1], "VALUE",vert[indexvert[0]], vert[indexvert[-1]]
         if indextime is not None: print "********** TIME: INDEX",indextime[0], indextime[-1], "VALUE",time[indextime[0]], time[indextime[-1]]

       ##var = nc.variables["phisinit"][:,:]
       ##contourf(np.transpose(var),30,cmap = get_cmap(name="Greys_r") ) ; mpl.pyplot.axis('off') ; plot(indexlat,indexlon,'mx',mew=4.0,ms=20.0)
       ##mpl.pyplot.show()
       ##exit()
       #truc = True 
       #truc = False
       #if truc: indexvert = None
       ####################################################################
       ########## REDUCE FIELDS
       ####################################################################
       error = False
       varname = all_varname[index_f]
       which=''
       if varname:   ### what is shaded.
           what_I_plot, error = myplot.reducefield( all_var[index_f], d4=indextime, d1=indexlon, d2=indexlat, d3=indexvert, \
                                             yint=yintegral, alt=vert, anomaly=anomaly, redope=redope, mesharea=area, unidim=is1d)
           if add != 0.:      what_I_plot = what_I_plot + add
           if mult != 2718.:  what_I_plot = what_I_plot*mult 
           else:              what_I_plot = np.log10(what_I_plot) ; print "log plot"

       if var2:      ### what is contoured.
           what_I_plot_contour, error = myplot.reducefield( all_var2[index_f], d4=indextime, d1=indexlon, d2=indexlat , d3=indexvert, \
                                                     yint=yintegral, alt=vert, redope=redope )
       if winds:     ### what is plot as vectors.
           vecx, error = myplot.reducefield( all_windu[index_f], d4=indextime, d3=indexvert, yint=yintegral, alt=vert)
           vecy, error = myplot.reducefield( all_windv[index_f], d4=indextime, d3=indexvert, yint=yintegral, alt=vert)
           if varname in [uchar,vchar]: what_I_plot = np.sqrt( np.square(vecx) + np.square(vecy) ) ; varname = "wind"
 
       if plotx is not None:
          plotx, error = myplot.reducefield( plotx, d4=indextime, d1=indexlon, d2=indexlat, d3=indexvert, \
                                             yint=yintegral, alt=vert, anomaly=anomaly, redope=redope, mesharea=area, unidim=is1d)
          ploty, error = myplot.reducefield( ploty, d4=indextime, d1=indexlon, d2=indexlat, d3=indexvert, \
                                             yint=yintegral, alt=vert, anomaly=anomaly, redope=redope, mesharea=area, unidim=is1d)
          which='xy'
       #####################################################################
       #if truc:
       #   nx = what_I_plot.shape[2] ; ny = what_I_plot.shape[1] ; nz = what_I_plot.shape[0] 
       #   for k in range(nz): print k,' over ',nz ; what_I_plot[k,:,:] = what_I_plot[k,:,:] / myplot.smooth(what_I_plot[k,:,:],12)
       #   for iii in range(nx):
       #    for jjj in range(ny):
       #     deviation = what_I_plot[:,jjj,iii] ; mx = max(deviation) ; mn = min(deviation)
       #     if iii > 6 and iii < nx-6 and jjj > 6 and jjj < ny-6:   what_I_plot[0,jjj,iii],rel = singlet(deviation,vert/1000.)  ### z must be in km
       #     else:                                                   what_I_plot[0,jjj,iii]     = 0.
       #     if np.abs(what_I_plot[0,jjj,iii]) > 1.5: 
       #         print iii,jjj,what_I_plot[0,jjj,iii],int(abs(1.-mx)*100.),int(abs(1.-mn)*100.)
       #         mpl.pyplot.plot(rel)
       #   mpl.pyplot.show()
       #   anomaly = True ### pour avoir les bons reglages plots
       #   what_I_plot = what_I_plot[0,:,:]  
       #####################################################################

       ####################################################################
       ### General plot settings
       changesubplot = (numplot > 1) and (len(what_I_plot.shape) != 1) and (which != "xy")  ## default for 1D plots: superimposed. to be reworked for better flexibility.
       if changesubplot: mpl.pyplot.subplot(subv,subh,nplot) #; mpl.pyplot.subplots_adjust(wspace=0,hspace=0)
       colorb = all_colorb[index_f]
       ####################################################################
       if error:
               myplot.errormess("There is an error in reducing field !")
       else:
               ticks = ndiv + 1 
               fvar = varname 
               if anomaly: fvar = 'anomaly'
               ###
               if mapmode == 0:    ### could this be moved inside imov loop ?
                   itime=indextime
                   if len(what_I_plot.shape) == 3: itime=[0]
                   m = None ; x = None ; y = None
                   latyeah = lat ; lonyeah = lon
                   if typefile in ['meso']:
                       # now that the section is determined we can set the real lat
                       # ... or for now, a temptative one.
                       milieux = int(lat2d.shape[1]/2.)
                       milieuy = int(lat2d.shape[0]/2.)
                       if slon is not None or proj == "noproj": latyeah = lat2d[:,milieux]
                       if slat is not None or proj == "noproj": lonyeah = lon2d[milieuy,:]
                   what_I_plot, x, y = myplot.define_axis(lonyeah,latyeah,vert,time,indexlon,indexlat,indexvert,\
                         itime,what_I_plot, len(all_var[index_f].shape),vertmode,redope)
               ###
               if analysis in ['laplace']: what_I_plot = laplace(what_I_plot)
               ###
               if (fileref is not None) and ((index_f+1)%3 == 0):    zevmin, zevmax = myplot.calculate_bounds(what_I_plot,vmin=minop,vmax=maxop)
               else:                                                   zevmin, zevmax = myplot.calculate_bounds(what_I_plot,vmin=vmin,vmax=vmax)
               #if (fileref is not None) and (index_f == numplot-1):    colorb = "RdBu_r"
               if colorb in ["def","nobar","onebar"]:                  palette = mpl.cm.get_cmap(name=myplot.defcolorb(fvar.upper()))
               elif colorb == "relief":                                palette = cm.GMT_relief
               elif colorb == "haxby":                                 palette = cm.GMT_haxby
               else:                                                   palette = mpl.cm.get_cmap(name=colorb)
               #palette = cm.GMT_split
               #palette = cm.GMT_globe
               ##### 1. ELIMINATE 0D or >3D CASES
               if len(what_I_plot.shape) == 0:   
                 print "VALUE VALUE VALUE VALUE ::: ",what_I_plot
                 save = 'donothing'
               elif len(what_I_plot.shape) >= 4:
                 print "WARNING!!! ",len(what_I_plot.shape),"-D PLOT NOT SUPPORTED !!! dimensions: ",what_I_plot.shape
                 myplot.errormess("Are you sure you did not forget to prescribe a dimension ?")
               ##### 2. HANDLE simple 1D/2D field and movies of 1D/2D fields
               else:
                 if mrate is not None: iend=len(time)-1
                 else:                 iend=0
                 imov = 0
                 if analysis in ['density','histo','fft','histodensity']: which="xy"
                 elif len(what_I_plot.shape) == 3:
                    if var2 and which == '':               which = "contour" ## have to start with contours rather than shading
                    elif which == '':                      which = "regular"
                    if mrate is None:      myplot.errormess("3D field. Use --rate RATE for movie or specify --time TIME. Exit.")
                 elif len(what_I_plot.shape) == 2:
                    if var2 and which == '':               which = "contour" ## have to start with contours rather than shading
                    elif which == '':                      which = "regular"
                    if mrate is not None and which == '':  which = "unidim"
                 elif len(what_I_plot.shape) == 1 and which == '' :
                    which = "unidim"
                    if what_I_plot.shape[-1] == 1:      print "VALUE VALUE VALUE VALUE ::: ", what_I_plot[0] ; save = 'donothing'
##                 if which == "unidim" and zelenfile > 1: numplot = 1 # this case is similar to several vars from one file 
                 ##### IMOV LOOP #### IMOV LOOP
                 while imov <= iend:
                    print "-> frame ",imov+1, which
                    if which == "regular":   
                        if mrate is None:                                   what_I_plot_frame = what_I_plot
                        else:                                               what_I_plot_frame = what_I_plot[imov,:,:]
                        if winds:
                            if mrate is None:                                   vecx_frame = vecx ; vecy_frame = vecy
                            else:                                               vecx_frame = vecx[imov,:,:] ; vecy_frame = vecy[imov,:,:]
                    elif which == "contour":  
                        if mrate is None or what_I_plot_contour.ndim < 3:   what_I_plot_frame = what_I_plot_contour
                        else:                                               what_I_plot_frame = what_I_plot_contour[imov,:,:]
                    elif which == "unidim":
                        if mrate is None:                                   what_I_plot_frame = what_I_plot
                        else:                                               what_I_plot_frame = what_I_plot[:,imov]  ## because swapaxes 
                    #if mrate is not None:     
                    if mapmode == 1: 
                        m = myplot.define_proj(proj,wlon,wlat,back=back,blat=blat,blon=blon)  ## this is dirty, defined above but out of imov loop
                        x, y = m(lon2d, lat2d)                                         ## this is dirty, defined above but out of imov loop
                    if (typefile in ['meso'] and mapmode == 1):
                       if '9999' not in getattr(nc,'START_DATE'): what_I_plot_frame = myplot.dumpbdy(what_I_plot_frame,6,condition=True)
#                   if typefile in ['mesoideal']:    what_I_plot_frame = myplot.dumpbdy(what_I_plot_frame,0,stag='W',condition=dumped_vert_stag)

                    if which == "unidim":
                        if lbls is not None: lbl=lbls[nplot-1] 
#                        if lbls is not None: lbl=lbls[index_f]
                        else:
                           lbl = ""
                           if indexlat is not None:  lbl = lbl + " ix" + str(indexlat[0])
                           if indexlon is not None:  lbl = lbl + " iy" + str(indexlon[0])
                           if indexvert is not None: lbl = lbl + " iz" + str(indexvert[0])
                           if indextime is not None: lbl = lbl + " it" + str(indextime[0])
                           if lbl == "": lbl = all_namefile[index_f]

                        if mrate is not None: x = y  ## because swapaxes...
                        #what_I_plot_frame = np.diff(what_I_plot_frame, n=1) ; x = x[1:]
                      
                        zeline = next(linecycler) ## "-" for simple lines
                        if tile:      zemarker = 'x'
                        else:         zemarker = None 
                        this_is_a_regular_plot = (indexvert is not None) or (indextime is None) or (indexlat is None) or (indexlon is None)
                        if this_is_a_regular_plot:   mpl.pyplot.plot(x,what_I_plot_frame,zeline,label=lbl,marker=zemarker)  ## vertical profile
                        else:                        mpl.pyplot.plot(what_I_plot_frame,x,zeline,label=lbl,marker=zemarker)  ## regular plot
                        mpl.pyplot.grid(True)
                        if nplot > 1: mpl.pyplot.legend(loc='best')
                        if indextime is None and axtime is not None and xlab is None:    mpl.pyplot.xlabel(axtime.upper()) ## define the right label
                        if save == 'txt':  writeascii(np.transpose([x,np.transpose(what_I_plot)]),'profile'+str(nplot*1000+imov)+'.txt')
                        if axtime == "lt" and indextime is None:
                            ax = mpl.pyplot.gca()
                            # set ticks where your images will be
                            ax.get_xaxis().set_ticks(np.arange(0,48,2))
                            # rename tick labels
                            ax.get_xaxis().set_ticklabels(["0","2","4","6","8","10","12","14","16","18","20","22",\
                                                       "0","2","4","6","8","10","12","14","16","18","20","22"])
                            ## rebound everyone
                            ax.set_xbound(lower=min(x), upper=max(x))
                    elif which == "xy":
                        if lbls is not None: lbl=lbls[index_f]
                        else: lbl=None
                        if analysis is not None:
                           if analysis == 'histo': 
                               if zelen == 1:
                                  mpl.pyplot.hist(ploty.flatten(),bins=ndiv,normed=True, alpha=0.67, facecolor = 'green', label = lbls)
                                  mpl.pyplot.legend()
                                  mpl.pyplot.grid(True)
                                  mpl.pyplot.title(zetitle)
                               else:
                                  multiplot[index_f]=ploty.flatten()
                                  if index_f == zelen-1: 
                                     if ope is not None: multiplot = getVar(multiplot,3*np.arange((index_f+1)/3)+2) ## we only compute histograms for the operation plots
                                     mpl.pyplot.hist(multiplot,bins=ndiv,normed=True, alpha=0.75, label = lbls) 
                                     mpl.pyplot.legend()
                                     mpl.pyplot.grid(True)
                                     mpl.pyplot.title(zetitle)
                           elif analysis in ['density','histodensity']:
                                  if ope is not None and (index_f+1)%3 !=0: pass
                                  else:
                                     plotx = np.linspace(min(ploty.flatten()),max(ploty.flatten()),1000)
                                     density = scipy.stats.gaussian_kde(ploty.flatten())
   #                                  density.covariance_factor = lambda : .25  # adjust the covariance factor to change the bandwidth if needed
   #                                  density._compute_covariance()
                                        # display the mean and variance of the kde:
                                     sample = density.resample(size=20000)
                                     n, (smin, smax), sm, sv, ss, sk = scipy.stats.describe(sample[0])
                                     mpl.pyplot.plot(plotx,density(plotx), label = lbl+'\nmean: '+str(sm)[0:5]+'   std: '+str(np.sqrt(sv))[0:5]+'\nskewness: '+str(ss)[0:5]+'   kurtosis: '+str(sk)[0:5])
                                     if analysis == 'histodensity':  # plot both histo and density (to assess the rightness of the kernel density estimate for exemple) and display the estimated variance
                                        mpl.pyplot.hist(ploty.flatten(),bins=ndiv,normed=True, alpha=0.30, label = lbl)
                                     if index_f == zelen-1: mpl.pyplot.legend() ; mpl.pyplot.title(zetitle)
                           else:
                              mpl.pyplot.plot(plotx,ploty,label = lbl)
                              if index_f == zelen-1: mpl.pyplot.legend() ; mpl.pyplot.title(zetitle)
                              mpl.pyplot.grid(True)
                        else:
                           mpl.pyplot.plot(plotx,ploty,label = lbl)
                           if index_f == zelen-1: mpl.pyplot.legend() ; mpl.pyplot.title(zetitle)
                        if varname == 'hodograph':
                            a=0
                            for ii in np.arange(len(time)): 
                               if a%6 == 0: mpl.pyplot.text(plotx[ii],ploty[ii],time[ii]) 
                               a=a+1
                            mpl.pyplot.grid(True)

                    elif which == "regular":
                    
                        # plot stream lines if there is a stream file and a vert/lat slice. Might not work with movies ??
                        if streamflag and sslat is None and svert is None:
                             streamfile = all_namefile[index_f].replace('.nc','_stream.nc')
                             teststream = os.path.exists(streamfile)
                             if teststream:
                                print 'INFO: Using stream file',streamfile, 'for stream lines'
                                ncstream = netCDF4.Dataset(streamfile)
                                psi = myplot.getfield(ncstream,'psi')
                                psi = psi[0,:,:,0] # time and longitude are dummy dimensions
                                if psi.shape[1] != len(x) or psi.shape[0] != len(y):
                                    myplot.errormess('stream file does not fit! Dimensions: '+str(psi.shape)+' '+str(x.shape)+' '+str(y.shape))
                                zelevels = np.arange(-1.e10,1.e10,1.e9)
                                zemin = np.min(abs(zelevels))
                                zemax = np.max(abs(zelevels))
                                zewidth  =  (abs(zelevels)-zemin)*(5.- 0.5)/(zemax - zemin) + 0.5 # linewidth ranges from 5 to 0.5
                                cs = mpl.pyplot.contour( x,y,psi, zelevels, colors='k', linewidths = zewidth)
                                mpl.pyplot.clabel(cs, inline=True, fontsize = 4.*mpl.pyplot.rcParams['font.size']/5., fmt="%1.1e")
                             else:
                                print 'WARNING: STREAM FILE',streamfile, 'DOES NOT EXIST !'
                             
                        if hole:         what_I_plot_frame = myplot.hole_bounds(what_I_plot_frame,zevmin,zevmax)
                        else:            what_I_plot_frame = myplot.bounds(what_I_plot_frame,zevmin,zevmax)
                        if flagnolow:    what_I_plot_frame = myplot.nolow(what_I_plot_frame)
                        if not tile:
                            #zelevels = np.linspace(zevmin*(1. + 1.e-7),zevmax*(1. - 1.e-7)) #,num=20)
                            zelevels = np.linspace(zevmin,zevmax,num=ticks)
                            #what_I_plot_frame = myplot.smooth(what_I_plot_frame,100)
                            if mapmode == 1:       m.contourf( x, y, what_I_plot_frame, zelevels, cmap = palette, alpha=trans)
                            elif mapmode == 0:     mpl.pyplot.contourf( x, y, what_I_plot_frame, zelevels, cmap = palette, alpha=trans)
                        else:
                            if mapmode == 1:       m.pcolor( x, y, what_I_plot_frame, cmap = palette, vmin=zevmin, vmax=zevmax, alpha=trans)
                            elif mapmode == 0:     mpl.pyplot.pcolor( x, y, what_I_plot_frame, cmap = palette, vmin=zevmin, vmax=zevmax, alpha=trans)

                        if (cross is not None or markdevil) and mapmode == 1:
                            if cross is not None: 
                               howmuch = np.array(cross).shape[0]
                               for ttt in range(howmuch):
                                  idx,idy=m(cross[ttt][0],cross[ttt][1])
                                  mpl.pyplot.plot([idx],[idy],'+k',mew=0.5,ms=18)
                            elif markdevil:
                                idx,idy=myplot.find_devil(nc,indextime)
                                idx,idy=x[idx,idy],y[idx,idy]
                                mpl.pyplot.plot([idx],[idy],'+k',mew=0.5,ms=18)

                        if not nocolorb:
                          if colorb not in ['nobar','onebar']:
                            if (fileref is not None) and ((index_f+1)%3 == 0):   daformat = "%.3f" 
                            elif mult != 1:                                        daformat = "%.1f"
                            else:                                                  daformat = myplot.fmtvar(fvar.upper())
                            if proj in ['moll']:  zeorientation="horizontal" ; zepad = 0.07 ; zefrac = 0.1 #zepad=0.05
                            elif proj in ['cyl']: zeorientation="vertical" ; zepad = 0.03 ; zefrac = 0.023
                            else:                 zeorientation="vertical" ; zepad = 0.03 ; zefrac = 0.05
                            if mapmode == 0:      zefrac = 0.1
                            zecb = mpl.pyplot.colorbar( fraction=zefrac,pad=zepad,format=daformat,orientation=zeorientation,\
                                      ticks=np.linspace(zevmin,zevmax,num=min([ticks/2+1,21])),extend='neither',spacing='proportional' ) 
                            #if zeorientation == "horizontal": 
                            #    daticks = zecb.ax.get_xticks()
                            #    zecb.ax.set_xticklabels(daticks,rotation='vertical')
                            if zeorientation == "horizontal" and zetitle[0] != "fill": zecb.ax.set_xlabel(zetitle[index_f]) ; zetitle[index_f]=""
                        if winds:
                            if typefile in ['meso']:
                                if '9999' not in getattr(nc,'START_DATE') : [vecx_frame,vecy_frame] = [myplot.dumpbdy(vecx_frame,6,stag=uchar,condition=True), myplot.dumpbdy(vecy_frame,6,stag=vchar,condition=True)]
                                key = True
                                if fvar in ['UV','uv','uvmet']: key = False
                            elif typefile in ['gcm']:
                                key = False
                            if metwind and mapmode == 1:   [vecx_frame,vecy_frame] = m.rotate_vector(vecx_frame, vecy_frame, lon2d, lat2d)
                            if trans != 0.0:   colorvec = myplot.definecolorvec(colorb) 
                            else:              colorvec = myplot.definecolorvec(back) 
                            myplot.vectorfield(vecx_frame, vecy_frame, x, y, stride=stride, csmooth=2,\
                                             #scale=15., factor=300., color=colorvec, key=key)
                                             scale=20., factor=250./facwind, color=colorvec, key=key)
                                                              #200.         ## or csmooth=stride
                        ### THIS IS A QUITE SPECIFIC PIECE (does not work for mesoscale files)
                        if ope == '-_histo' and nplot == numplot: # this should work as long as ope is '-' guarantees 3 plots for 4 panels without contour
                            mpl.pyplot.subplot(subv,subh,nplot+1)
                            mpl.pyplot.rcParams["legend.fontsize"] = 'xx-large'
                            if indexlat is None:
                                latmin = -50.; latmax = 50. # latitude range for histogram of difference
                                zeindexlat = (lat<latmax)*(lat>latmin)
                                if typefile in ['meso']: zeindexlat = 10
                                # this follows the define_axis logic in myplot.py:
                                if indextime is None or indexlon is None: what_I_plot_frame = what_I_plot_frame[zeindexlat,:]
                                else: what_I_plot_frame = what_I_plot_frame[:,zeindexlat]
                            toplot = np.ravel(what_I_plot_frame[np.isnan(what_I_plot_frame)==False])
                            zebins = np.linspace(zevmin,zevmax,num=30)
                            mpl.pyplot.hist(toplot,bins=zebins,histtype='step',linewidth=2,color='k',normed=True)
                            zestd = np.std(toplot);zemean = np.mean(toplot)
                            zebins = np.linspace(zevmin,zevmax,num=300)
                            zegauss = (1./(zestd * np.sqrt(2 * np.pi)) * np.exp( - (zebins - zemean)**2 / (2 * zestd**2) ) )
                            mpl.pyplot.plot(zebins, zegauss, linewidth=1, color='r',label="mean: "+str(zemean)[0:5]+"\nstd: "+str(zestd)[0:5])
                            mpl.pyplot.legend(loc=0,frameon=False)
                            mpl.pyplot.subplot(subv,subh,nplot) # go back to last plot for title of contour difference
                        if ope is not None and "only" not in ope: mpl.pyplot.title("fig(1) "+ope+" fig(2)")
                        elif ope is not None and "only" in ope: mpl.pyplot.title("fig(1) "+ope[0]+" fig(2)")
                            
                    elif which == "contour":
                        #mpl.pyplot.rcParams['contour.negative_linestyle'] = 'solid' # no dashed line for negative values
                        zevminc, zevmaxc = myplot.calculate_bounds(what_I_plot_frame, vmin=min(what_I_plot_frame), vmax=max(what_I_plot_frame))
                        zelevels = np.linspace(zevminc,zevmaxc,ticks/2) #20)
                        ### another dirty specific stuff in the wall
                        if var2 == 'HGT':        zelevels = np.arange(-10000.,30000.,250.) #1000.)
                        elif var2 == 'tpot':     zelevels = np.arange(270,370,5)
                        elif var2 in ['tk','temp']:     zelevels = np.arange(0.,1000.,5.)#(150,250,5)
                        elif var2 == 'wstar':    zelevels = np.arange(0,10,1.0)
                        elif var2 == 'zmax_th':  zelevels = np.arange(0,10,2.0) ; what_I_plot_frame = what_I_plot_frame / 1000.
                        elif var2 in ['u']:      zelevels = np.arange(-400.,400.,5.)
                        ###
                        if mapmode == 0:   
                            what_I_plot_frame, x, y = myplot.define_axis( lonyeah,latyeah,vert,time,indexlon,indexlat,indexvert,\
                                                              itime,what_I_plot_frame, len(all_var2[index_f].shape),vertmode,redope)
                            ## is this needed? only if len(all_var2[index_f].shape) != len(all_var[index_f].shape)
                            cs = mpl.pyplot.contour( x,y,what_I_plot_frame, zelevels, colors='k', linewidths = 0.33)#, alpha=0.5, linestyles='solid')
                            ##cs = mpl.pyplot.contour( x,y,what_I_plot_frame, zelevels, colors='w', linewidths = 0.5)#, alpha=0.5, linestyles='solid')
                            #mpl.pyplot.clabel(cs, inline=1, fontsize = 4.*rcParams['font.size']/5., fmt=fmtvar(var2.upper()))
                        elif mapmode == 1:  
                            cs = m.contour( x,y,what_I_plot_frame, zelevels, colors='k', linewidths = 0.33)#, alpha=0.5, linestyles='solid')
                            #mpl.pyplot.clabel(cs, inline=0, fontsize = mpl.pyplot.rcParams['font.size'], fmt="%.0f") #myplot.fmtvar(var2.upper()))
                    if which in ["regular","unidim","xy"]:

                        if nplot > 1 and which in ["unidim","xy"]:
                           pass  ## because we superimpose nplot instances
                        else:
                           # Axis directives for movie frames [including the first one).
                           zxmin, zxmax = xaxis ; zymin, zymax = yaxis
                           if zxmin is not None: mpl.pyplot.xlim(xmin=zxmin)
                           if zxmax is not None: mpl.pyplot.xlim(xmax=zxmax)
                           if zxmin is not None and zxmax is not None:
                               ax = mpl.pyplot.gca()
                               ax.get_xaxis().set_ticks(np.linspace(xaxis[0],xaxis[1],ticks/2+1)) 
                           if zymin is not None: mpl.pyplot.ylim(ymin=zymin)
                           if zymax is not None: mpl.pyplot.ylim(ymax=zymax)
                           if ylog and not xlog:      mpl.pyplot.semilogy()
                           if xlog and not ylog:      mpl.pyplot.semilogx()
                           if xlog and ylog: 
                                mpl.pyplot.xscale('log') 
                                mpl.pyplot.yscale('log')
                           if invert_y:  ax = mpl.pyplot.gca() ; ax.set_ylim(ax.get_ylim()[::-1])
                           if xlab is not None: mpl.pyplot.xlabel(xlab)
                           if ylab is not None: mpl.pyplot.ylabel(ylab)

                        if mrate is not None:
                           ### THIS IS A MENCODER MOVIE
                           if mrate > 0:
                             figframe=mpl.pyplot.gcf()
                             if mquality:   figframe.set_dpi(600.)
                             else:          figframe.set_dpi(200.)
                             mframe=fig2img(figframe)
                             if imov == 0:
                                moviename='movie' ;W,H = figframe.canvas.get_width_height()
                                video = VideoSink((H,W), moviename, rate=mrate, byteorder="rgba")
                             video.run(mframe) ; mpl.pyplot.close()
                             if imov == iend: video.close()                            
                           ### THIS IS A WEBPAGE MOVIE 
                           else:
                             nameframe = "image"+str(1000+imov)
                             myplot.makeplotres(nameframe,res=100.,disp=False) ; mpl.pyplot.close()
                             if imov == 0: myfile = open("zepics", 'w')
                             myfile.write("modImages["+str(imov)+"] = '"+nameframe+"_100.png';"+ '\n')
                             if imov == iend:
                                 myfile.write("first_image = 0;"+ '\n')
                                 myfile.write("last_image = "+str(iend)+";"+ '\n')
                                 myfile.close()
                        if var2 and which == "regular":  which = "contour"
                        imov = imov+1
                    elif which == "contour":
                        which = "regular"

       ### Next subplot
       zevarname = varname
       if redope is not None: zevarname = zevarname + "_" + redope 
       basename = myplot.getname(var=zevarname,var2=var2,winds=winds,anomaly=anomaly)
       if len(what_I_plot.shape) > 3:
           basename = basename + myplot.getstralt(nc,level) 
       if mrate is not None: basename = "movie_" + basename 
       if typefile in ['meso']:
            if sslon is not None: basename = basename + "_lon_" + str(int(round(lonp)))
            if sslat is not None: basename = basename + "_lat_" + str(int(round(latp)))
            plottitle = basename
            ### dans le nouveau systeme time=ls,sol,lt cette ligne pourrait ne servir a rien (ou deplacer au dessus)
            if addchar and indextime is not None:   [addchar,gogol,gogol2] = myplot.getlschar ( all_namefile[index_f] )  ;  plottitle = plottitle + addchar
            ### en fait redope is None doit etre remplace par : n'est ni maxt ni mint
            if redope is None and ltst is not None and ( (mapmode == 0) or (proj in ["lcc","laea","merc","nsper"]) ):  plottitle = plottitle + "_LT" + str(ltst)
       else:
            if fileref is not None:
               if (index_f+1)%3 == 0:     plottitle = basename+' '+"fig(1) "+ope+" fig(2)"
               elif (index_f+2)%3 == 0:   plottitle = basename+' '+fileref
               else:                        plottitle = basename+' '+all_namefile[index_f]
            else:                            plottitle = basename+' '+all_namefile[index_f]
       if mult != 1:                         plottitle = '{:.0e}'.format(mult) + "*" + plottitle
       if zetitle[0] != "fill":                 
          if index_f < len(zetitle): plottitle = zetitle[index_f]
          else: plottitle = zetitle[0]
          if titleref is "fill":             titleref=zetitle[index_f]
          if fileref is not None and which != "xy":
             if (index_f+2)%3 == 0:        plottitle = titleref
             if (index_f+1)%3 == 0:        plottitle = "fig(1) "+ope+" fig(2)"
#       if indexlon is not None:      plottitle = plottitle + " lon: " + str(min(lon[indexlon])) +" "+ str(max(lon[indexlon]))
#       if indexlat is not None:      plottitle = plottitle + " lat: " + str(min(lat[indexlat])) +" "+ str(max(lat[indexlat]))
#       if indexvert is not None:     plottitle = plottitle + " vert: " + str(min(vert[indexvert])) +" "+ str(max(vert[indexvert]))
#       if indextime is not None:     plottitle = plottitle + " time: " + str(min(time[indextime])) +" "+ str(max(time[indextime]))
       if colorb != "onebar": mpl.pyplot.title( plottitle )
       if nplot >= numplot: error = True
       nplot += 1

       if zelenfile > 1 and zelenvar == 1 and which == "unidim": index_f=index_f+1
       firstpass=False

    if colorb == "onebar":
        cax = mpl.pyplot.axes([0.1, 0.2, 0.8, 0.03]) # a ameliorer
        zecb = mpl.pyplot.colorbar(cax=cax, orientation="horizontal", format=myplot.fmtvar(fvar.upper()),\
                 ticks=np.linspace(zevmin,zevmax,num=min([ticks/2+1,21])),extend='neither',spacing='proportional')
        if zetitle[0] != "fill": zecb.ax.set_xlabel(zetitle[index_f]) ; zetitle[index_f]=""


    ##########################################################################
    ### Save the figure in a file in the data folder or an user-defined folder
    if outputname is None:
       if typefile in ['meso']:   prefix = myplot.getprefix(nc)
       elif typefile in ['gcm']:            prefix = 'LMD_GCM_'
       else:                                prefix = ''
    ###
       zeplot = prefix + basename 
       if zoom:            zeplot = zeplot + "zoom"+str(abs(zoom))
       if addchar:         zeplot = zeplot + addchar
       if numplot <= 0:    zeplot = zeplot + "_LT"+str(abs(numplot))
    ###
       if not target:      zeplot = namefile[0:find(namefile,'wrfout')] + zeplot
       else:               zeplot = target + "/" + zeplot  
    ###
    else:
       zeplot=outputname

    if mrate is None:
        pad_inches_value = 0.35
        if wlon[1]-wlon[0] < 2.: pad_inches_value = 0.5  # LOCAL MODE (small values)
        print "********** SAVE ", save
        if save == 'png': 
            if display: myplot.makeplotres('.thumb',res=50.,pad_inches_value=pad_inches_value) #,erase=True)  ## a miniature
            myplot.makeplotres(zeplot,res=resolution,pad_inches_value=pad_inches_value,disp=False)
        elif save in ['eps','svg','pdf']:     myplot.makeplotres(zeplot,pad_inches_value=pad_inches_value,disp=False,ext=save)
        elif save == 'gui':                   mpl.pyplot.show()
        elif save == 'donothing':             pass
        elif save == 'txt':                   print "Saved results in txt file." 
        else: 
            print "INFO: save mode not supported. using gui instead."
            mpl.pyplot.show()

    ###################################
    #### Getting more out of this video -- PROBLEMS WITH CREATED VIDEOS
    #
    if mrate is not None and save != "html":
        print "Re-encoding movie.. first pass"
        video.first_pass(filename=zeplot,quality=mquality,rate=mrate)
    #    print "Re-encoding movie.. second pass"
    #    video.second_pass(filename=moviename,quality=mquality,rate=mrate)   

    ###############
    ### Now the end
    return zeplot

#!/usr/bin/env python

### A. Colaitis

##
#   This routine transforms a diagfi.nc file into a diagfi_MCS.nc file where
#   the fields are directly comparable to those contained in MCS data, i.e. 
#   fields are re-binned at times over the ranges specified in the MCS file.
###

###########################################################################################
###########################################################################################
### What is below relate to running the file as a command line executable (very convenient)
if __name__ == "__main__":
   import sys
   from optparse import OptionParser    ### to be replaced by argparse
   from netCDF4 import Dataset
   from os import system,path
   from times import sol2ls
   import numpy as np
   from mymath import find_nearest 
   from myplot import getfield,separatenames
   from make_netcdf import make_gcm_netcdf
   from gcm_transformations import call_zrecast,call_hrecast
   parser = OptionParser() 

   #############################
   ### Options
   parser.add_option('-f', '--file',   action='store',dest='file',     type="string",  default=None,  help='[NEEDED] filename.')
   parser.add_option('-m', '--mfile',  action='store',dest='mcsfile',     type="string",  default=None,  help='[NEEDED] filename for MCS comparison.')
   parser.add_option('-v', '--var',    action='append',dest='var',      type="string",  default=None,  help='[NEEDED] Variables to process. (coma-separated list. aps and bps are always included.)')
   parser.add_option('-x', action='store_false',dest='recast',     default=True,  help='Force aps and bps to be included in output file (usefull if your file is already recasted along z) [True]')
   parser.add_option('-i', '--zrecast', action='store_true', dest='zrecast', default=False, help='Cast zrecast.e on diagfi file with MCS pressure levels. Will pass this operation is recasted file is already present, unless --override is specified. [False]')
   parser.add_option('-H', '--hrecast', action='store_true', dest='hrecast', default=False, help='Cast hrecast.e on diagfi file on MCS lat/lon grid. Will pass this operation is recasted file is already present, unless --override is specified. [False]')
   parser.add_option('--override', action='store_true', dest='override', default=False, help='Force zrecast.e to act even if recasted file is already present(will erase previous recasted file) [False]')
   parser.add_option('--ditch', action='store_true', dest='ditch', default=False, help='Ditch recasted file when interpolation is performed. [False]')
   parser.add_option('--latreverse', action='store_true', dest='latreverse', default=False, help='Reverse the latitude axis in output diagfi. [False]')

   #############################
   ### Get options and variables
   (opt,args) = parser.parse_args()

   #############################
   ### Load and check data

   if opt.var is None:
      print "You must specify at least a field to process with -v."
      exit()

   # Hrecast & Zrecast

   varznames=separatenames(opt.var[0])
   filename=opt.file

#   if opt.hrecast:
#      if (path.exists(filename[0:len(filename)-3]+"_h.nc") and (not opt.override)):
#         print "--> "+filename[0:len(filename)-3]+"_h.nc"
#         print "Recasted file is already there, skipping interpolation. [use --override to force interpolation]"
#         filename=filename[0:len(filename)-3]+"_h.nc"
#      else:
#         print "--> "+filename[0:len(filename)-3]+"_h.nc"
#         filename=call_hrecast (  input_name      = [filename], \
#                    fields  = varznames, \
#                    predefined = 'mcs')[0]

   if opt.zrecast:
      if (path.exists(filename[0:len(filename)-3]+"_P.nc") and (not opt.override)):
         print "--> "+filename[0:len(filename)-3]+"_P.nc"
         print "Recasted file is already there, skipping interpolation. [use --override to force interpolation]"
         filename=filename[0:len(filename)-3]+"_P.nc"
      else:
         print "--> "+filename[0:len(filename)-3]+"_P.nc"
         filename=call_zrecast (  interp_mode   = 2, \
                    input_name      = [filename], \
                    fields  = varznames, \
                    predefined = 'mcs')[0]

   if opt.hrecast:
      if (path.exists(filename[0:len(filename)-3]+"_h.nc") and (not opt.override)):
         print "--> "+filename[0:len(filename)-3]+"_h.nc"
         print "Recasted file is already there, skipping interpolation. [use --override to force interpolation]"
         filename=filename[0:len(filename)-3]+"_h.nc"
      else:
         print "--> "+filename[0:len(filename)-3]+"_h.nc"
         filename=call_hrecast (  input_name      = [filename], \
                    fields  = varznames, \
                    predefined = 'mcs')[0]


   # Files

   print "--> Loading diagfi dataset."

   nc=Dataset(filename)
   ncmcs=Dataset(opt.mcsfile)

   # Dimensions

   lon=nc.variables["longitude"][:]
   lat=nc.variables["latitude"][:]
   alt=nc.variables["altitude"][:]
   time=nc.variables["Time"][:] # in fraction of sols
   if "controle" in nc.variables:
      controle=nc.variables["controle"][:]
      day_ini=controle[3]%669
   else:
      if opt.zrecast:
         nccontrol=Dataset(opt.file)
         if "controle" in nccontrol.variables:
            controle=nccontrol.variables["controle"][:]
            day_ini=controle[3]%669
         else:
            print "Error: could not find controle variable in diagfi."
            day_ini=input("Please type initial sol number:")%669
      else:
         print "Error: could not find controle variable in diagfi."
         day_ini=input("Please type initial sol number:")%669
   time[:]=time[:]+day_ini
   nx=len(lon)
   ny=len(lat)
   nz=len(alt)
   nt=len(time)
   lstime=sol2ls(time)

   # MCS

   print "--> Loading and preparing MCS dataset."

   dtimemintmp=ncmcs.variables["dtimemin"][:,:,:]
   dtimemaxtmp=ncmcs.variables["dtimemax"][:,:,:]
   ntimemintmp=ncmcs.variables["ntimemin"][:,:,:]
   ntimemaxtmp=ncmcs.variables["ntimemax"][:,:,:]
   lonmcs=ncmcs.variables["longitude"][:]
   latmcs=ncmcs.variables["latitude"][:]
   timemcs=ncmcs.variables["time"][:]%360 # IN LS

   dtimemin=np.ma.masked_where(dtimemintmp < 0.,dtimemintmp)
   dtimemin.set_fill_value([np.NaN])
   dtimemax=np.ma.masked_where(dtimemaxtmp < 0.,dtimemaxtmp)
   dtimemax.set_fill_value([np.NaN])
   ntimemin=np.ma.masked_where(ntimemintmp < 0.,ntimemintmp)
   ntimemin.set_fill_value([np.NaN])
   ntimemax=np.ma.masked_where(ntimemaxtmp < 0.,ntimemaxtmp)
   ntimemax.set_fill_value([np.NaN])

   # Variables to treat

   print "--> Preparing diagfi dataset."

   varz=[]
   n=0
   for zn in varznames:
       load=getfield(nc,zn)
       load=np.ma.masked_where(load < -1.e-20,load)
       load.set_fill_value([np.NaN])
       load=load.filled()
       load=np.ma.masked_invalid(load)
       load.set_fill_value([np.NaN])
       load=load.filled()
       varz.append(load)
       load=0.
       print "Found: "+zn+" with dimensions: "
       print np.array(varz[n]).shape
       n=n+1

   nzvar=len(varz)
   dimensions={}
   vv=0
   for var in varz:
       a=len(np.array(var).shape)
       if a == 3: dimensions[vv]=a
       elif a == 4: dimensions[vv]=a
       else:
          print "Warning, only 3d and 4d variables are supported for time-recasting"
          exit()
       vv=vv+1

   # Variables to save but not treated (only along z, or phisinit-like files)

   aps=nc.variables["aps"][:]
   bps=nc.variables["bps"][:]
   fullnames=["aps","bps"]
   for name in varznames:
       fullnames.append("d"+name)
       fullnames.append("n"+name)
   print "Will output: "
   if opt.recast: print fullnames[2:]
   else: print fullnames
   #############################
   ### Building
   #############################

   ### We loop over chunks of gcm data corresponding to MCS time dimension
   ### Bin sizes for mcs data is 5 degrees ls centered on value
   varday=np.zeros([len(timemcs),nz,ny,nx])
   varnight=np.zeros([len(timemcs),nz,ny,nx])
   vardayout=np.zeros([nzvar,len(timemcs),nz,ny,nx])
   varnightout=np.zeros([nzvar,len(timemcs),nz,ny,nx])
   vardayout=np.ma.masked_invalid(vardayout)
   varnightout=np.ma.masked_invalid(varnightout)
   i=0
   for ls in timemcs:
       lsstart=ls-2.5
       lsstop=ls+2.5
       istart=find_nearest(lstime,lsstart,strict=True)
       istop=find_nearest(lstime,lsstop,strict=True)
       varchk=0
       if ((istart is np.NaN) or (istop is np.NaN)):
          vardayout[:,i,:,:,:]=np.NaN
          varnightout[:,i,:,:,:]=np.NaN
          print "Time interval skipped. Ls MCS: (",lsstart,';',lsstop,')',"// Ls diagfi: (",lstime.min(),';',lstime.max(),')'
          i=i+1
          continue
       print "--->> Processing Data. Ls MCS: (",lsstart,';',lsstop,')',"// Ls diagfi: (",lstime.min(),';',lstime.max(),')'
       # warning, python's convention is that the second index of array[a:b] is the array index of element after the one being picked last.
       # for that reason, array[0:0] is nan and array[0:1] is only one value. Hence, len(array[a:b+1]) is b-a+1 and not b-a+2
       print "     .initialisation."
       
       
       varchk=np.zeros([nzvar,istop-istart+1,nz,ny,nx])
       vv=0
       for variable in varz:
           if dimensions[vv] is 3:
              varchk[vv,:,0,:,:]=variable[istart:istop+1,:,:]
           else:
              varchk[vv,:,:,:,:]=variable[istart:istop+1,:,:,:]
           vv=vv+1
       varchk=np.ma.masked_invalid(varchk)
       varchk.set_fill_value([np.NaN])
       varchktime=time[istart:istop+1]
       ndays=np.floor(varchktime[len(varchktime)-1])-np.floor(varchktime[0])
       dtmichk=dtimemin[i,:,:]
       dtmachk=dtimemax[i,:,:]
       ntmichk=ntimemin[i,:,:]
       ntmachk=ntimemax[i,:,:]
       dtmichk.set_fill_value([np.NaN])
       dtmachk.set_fill_value([np.NaN])
       ntmichk.set_fill_value([np.NaN])
       ntmachk.set_fill_value([np.NaN])
       dtmichk=dtmichk.filled()
       dtmachk=dtmachk.filled()
       ntmichk=ntmichk.filled()
       ntmachk=ntmachk.filled()

   ### We iterate for each day in the chunk, on each grid point we find
   ### the closest corresponding MCS grid point and the index of the 
   ### time in the chunk closest to the time in the closest MCS grid point.
   ### (yea it's complicated)

       vartmpnight=np.zeros([nzvar,ndays,nz,ny,nx])
       vartmpday=np.zeros([nzvar,ndays,nz,ny,nx])
       vartmpnight=np.ma.masked_invalid(vartmpnight)
       vartmpday=np.ma.masked_invalid(vartmpday)
       vartmpnight.set_fill_value([np.NaN])
       vartmpday.set_fill_value([np.NaN])

       nd=0
       print "     .time indices MCS grid -> diagfi grid."
       while nd < ndays:

          daystart=find_nearest(varchktime-varchktime[0],nd)
          daystop=find_nearest(varchktime-varchktime[0],nd+1)
#          varchktime_lon=np.zeros([daystop-daystart+1,len(lon)])
          ix=0
          for x in lon:

             varchktime_lon = 24.*(varchktime[daystart:daystop+1]-varchktime[daystart]) + x/15.

             iy=0
             for y in lat:
                niy=find_nearest(latmcs,y)
                nix=find_nearest(lonmcs,x)
                nitdtmichk=find_nearest(varchktime_lon,dtmichk[niy,nix])
                nitdtmachk=find_nearest(varchktime_lon,dtmachk[niy,nix])
                nitntmichk=find_nearest(varchktime_lon,ntmichk[niy,nix])
                nitntmachk=find_nearest(varchktime_lon,ntmachk[niy,nix])
                for vv in np.arange(nzvar):
                   if ((nitdtmichk is np.NaN) or (nitdtmachk is np.NaN)):
                       vartmpday[vv,nd,:,iy,ix]=np.NaN
                   elif nitdtmichk > nitdtmachk:
                       vartmpday[vv,nd,:,iy,ix]=(np.ma.mean(varchk[vv,daystart+nitdtmichk:daystop+1,:,iy,ix],axis=0)+np.ma.mean(varchk[vv,daystart:daystart+nitdtmachk+1,:,iy,ix],axis=0))/2.
                   else:
                       vartmpday[vv,nd,:,iy,ix]=np.ma.mean(varchk[vv,daystart+nitdtmichk:daystart+nitdtmachk+1,:,iy,ix],axis=0)
                   if ((nitntmichk is np.NaN) or (nitntmachk is np.NaN)):
                       vartmpnight[vv,nd,:,iy,ix]=np.NaN
                   elif nitntmichk > nitntmachk:
                       vartmpnight[vv,nd,:,iy,ix]=(np.ma.mean(varchk[vv,daystart+nitntmichk:daystop+1,:,iy,ix],axis=0)+np.ma.mean(varchk[vv,daystart:daystart+nitntmachk+1,:,iy,ix],axis=0))/2.
                   else:                                                            
                       vartmpnight[vv,nd,:,iy,ix]=np.ma.mean(varchk[vv,daystart+nitntmichk:daystart+nitntmachk+1,:,iy,ix],axis=0)
                iy=iy+1
             ix=ix+1
          nd=nd+1

       print "     .creating bins."

       vartmpdaymasked=np.ma.masked_invalid(vartmpday)
       vartmpnightmasked=np.ma.masked_invalid(vartmpnight)
       vartmpdaymasked.set_fill_value([np.NaN])
       vartmpnightmasked.set_fill_value([np.NaN])
       for vv in np.arange(nzvar):
          vardayout[vv,i,:,:,:]=np.ma.mean(vartmpdaymasked[vv,:,:,:,:],axis=0)
          varnightout[vv,i,:,:,:]=np.ma.mean(vartmpnightmasked[vv,:,:,:,:],axis=0)
          print "          ."+varznames[vv]+".done"
       i=i+1

   print "--->> Preparing Data for ncdf. Missing value is NaN."

   vardayout=np.ma.masked_invalid(vardayout)
   varnightout=np.ma.masked_invalid(varnightout)
   vardayout.set_fill_value([np.NaN])
   varnightout.set_fill_value([np.NaN])

   if opt.latreverse:
       vardayout[:,:,:,:,:]=vardayout[:,:,:,::-1,:]
       varnightout[:,:,:,:,:]=varnightout[:,:,:,::-1,:]

   all=[aps,bps]
   for vv in np.arange(nzvar):
       if dimensions[vv] == 3:
          all.append(vardayout[vv,:,0,:,:].filled())
          all.append(varnightout[vv,:,0,:,:].filled())
       elif dimensions[vv] == 4:
          all.append(vardayout[vv,:,:,:,:].filled())
          all.append(varnightout[vv,:,:,:,:].filled())

   if opt.recast:
      all=all[2:]
      fullnames=fullnames[2:] 

   if opt.latreverse:
      lat=lat[::-1]

   make_gcm_netcdf (zfilename=filename[0:len(filename)-3]+"_MCS.nc", \
                        zdescription="Temperatures from diagfi reworked to match MCS format", \
                        zlon=lon, \
                        zlat=lat, \
                        zalt=alt, \
                        ztime=timemcs, \
                        zvariables=all, \
                        znames=fullnames)
   if opt.zrecast and opt.ditch:
      print "removing interpolated file"
      system("rm -f "+opt.file[0:len(opt.file)-3]+"_P.nc")

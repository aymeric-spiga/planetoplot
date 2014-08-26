### A. Colaitis -- LMD -- 11/11/2011

def make_gcm_netcdf (zfilename="filename.nc", \
                     zdescription=None, \
                     zlon=None, \
                     zlat=None, \
                     zalt=None, \
                     ztime=None, \
                     zvariables=None, \
                     znames=None):

   from os import system
   from netCDF4 import Dataset
   import numpy as np

#####################
## Run some checks ##
#####################

# Number of variables and number of variables' names

   nvar = len(zvariables)
   nnames = len(znames)

   if nnames < nvar:
      print "not enough variables name specified, using automatic labbeling"
      i=0
      zznames=[[]]*nvar
      for zvar in zvariables:
          zznames[i]="var"+str(i+1)
          i=i+1
   
   elif nnames > nvar:
      print "Warning, more variables than names are specified."
      i=0
      zznames=[[]]*nvar
      for zvar in zvariables:
          zznames[i]=znames[i]
          i=i+1
   else:
      zznames=znames[:]
   
   # Dimensions
   
   nx,ny,nz,nt = 0,0,0,0
   
   if zlon is not None:
      nx=len(zlon)
   if zlat is not None:
      ny=len(zlat)
   if zalt is not None:
      nz=len(zalt)
   if ztime is not None:
      nt=len(ztime)
   
   zdims={}
   zdims['longitude']=nx
   zdims['latitude']=ny
   zdims['altitude']=nz
   zdims['Time']=nt
   print zdims 
   # Find which variable uses which dimensions
   
   i=0
   zvarcarac={}
   for zvar in zvariables:
       zvardim=np.array(zvar).shape
       ndim=len(zvardim)
       zzvardim=[[]]*ndim
       j=0
       for dim in zvardim:
           if dim not in zdims.values():
              print "WARNING -----------------------------"
              print "Dimensions given to subroutine do not match variables dimensions :"
              print "Dimensions: ",zdims
              print "Variable: ",zznames[i],", dimensions: ",zvardim
              print " PROGRAM EXIT"
              exit()
           else:
              a=get_key(zdims,dim)
              if len(a) is not 1:
                 if j is 0:                ##this should solve most conflicts with Time
                    zzvardim[j]=a[1]
                 elif j is 3:
                    zzvardim[j]=a[1]
                 else:
                    zzvardim[j]=a[0]
              else:
                 zzvardim[j]=a[0]
              j=j+1
       zvarcarac[zznames[i]]=zzvardim
       i=i+1

   print "creating "+zfilename

   print zvarcarac  
 
   #########################
   ## Clean previous file ##
   #########################
   
   system("rm -f "+zfilename)
   
   ## Open file
   file = Dataset(zfilename, 'w', format='NETCDF3_64BIT') #netcdf4 does not work with ncview and ncdump yet
   if zdescription is not None:
      file.description = zdescription
   ## Dimensions
   if ztime is not None:
      file.createDimension('Time', None)
   if zalt is not None:
      file.createDimension('altitude',nz)
   if zlat is not None:
       file.createDimension('latitude',ny)
   if zlon is not None:
      file.createDimension('longitude',nx)
   ## Variables for dimensions & Data assignment
   if ztime is not None:
      times = file.createVariable('Time', 'f', ('Time',))
      times[:] = ztime[:]
   if zalt is not None:
      altitudes = file.createVariable('altitude', 'f', ('altitude',))
      altitudes[:] = zalt[:]
   if zlat is not None:
      latitudes = file.createVariable('latitude', 'f', ('latitude',))
      latitudes[:] = zlat[:]
   if zlon is not None:
      longitudes = file.createVariable('longitude', 'f', ('longitude',))
      longitudes[:] = zlon[:]
   ## Other Variables Creation & Data assignment
   i=0
   for name in zznames:
       za = file.createVariable(name, 'f', tuple(zvarcarac[name]))
       if len(zvarcarac[name]) is 1:
          za[:] = zvariables[i][:]
       elif len(zvarcarac[name]) is 2:
          za[:,:] = zvariables[i][:,:]
       elif len(zvarcarac[name]) is 3:
          za[:,:,:] = zvariables[i][:,:,:]
       elif len(zvarcarac[name]) is 4:
          za[:,:,:,:] = zvariables[i][:,:,:,:]
       i=i+1
   ## close file
   file.close()

   print "closing "+zfilename

   return

def find_key(dic, val):
    """return the key of dictionary dic given the value"""
    return [k for k, v in dic.iteritems() if v == val][0]

def get_key(self, value):
    """find the key(s) as a list given a value"""   
    return [item[0] for item in self.items() if item[1] == value]

###############################################
## PLANETOPLOT                               ##
## --> PPCOMPUTE                             ##
###############################################
# python built-in librairies
import os
# added librairies
import numpy as np
import math  as m
###############################################

## first a useful function to find settings in a folder in PYTHONPATH
def findset(whereset,string="planetoplot"):
    # ... can set a local set of settings in .settings
    cwd = os.getcwd()
    islocal = os.path.isdir(cwd+"/.settings")
    if islocal:
      whereset = cwd+"/.settings"
    else:
      # ... set a default whereset if it was set to None
      # ... default is in the planetoplot folder
      if whereset is None:
        for path in os.environ['PYTHONPATH'].split(os.pathsep):
            if string in path: whereset = path + "/../settings/"
        if whereset is None:
            print "!! WARNING !! "+ string + " not in $PYTHONPATH"
            print "--> set to current directory instead"
            whereset = "./"
    # ... if the last / is missing put it
    if whereset[-1] != "/": whereset = whereset + "/"
    return whereset

##########################
#### MAX MEAN MIN SUM ####
#####################################
#### WITH SUPPORT FOR NaN VALUES ####
#####################################

## compute log10 of field
## author DB + JVO
def log10(field): 
        if field is None: return None
        if type(field).__name__=='MaskedArray':
              field.set_fill_value(np.NaN)
              return np.ma.array(np.log10(field))
#        elif (np.isnan(field) and (type(field).__name__ not in 'MaskedArray')):
#              return np.ma.masked_invalid(np.log10(field))
        else: return np.log10(field)

## compute max
## compute min
## author AS + AC
def min(field,axis=None): 
        if field is None: return None
        if type(field).__name__=='MaskedArray':
              field.set_fill_value(np.NaN)
              return np.ma.array(field).min(axis=axis)
        elif (np.isnan(np.sum(field)) and (type(field).__name__ not in 'MaskedArray')):
              return np.ma.masked_invalid(field).min(axis=axis)
        else: return np.array(field).min(axis=axis)

## compute max
## author AS + AC
def max(field,axis=None):
        if field is None: return None
        if type(field).__name__=='MaskedArray':
              field.set_fill_value(np.NaN)
              return np.ma.array(field).max(axis=axis)
        elif (np.isnan(np.sum(field)) and (type(field).__name__ not in 'MaskedArray')):
              return np.ma.masked_invalid(field).max(axis=axis) 
        else: return np.array(field).max(axis=axis)

## compute mean
## author AS + AC
def mean(field,axis=None):
        if field is None: return None
        else: 
           if type(field).__name__=='MaskedArray':
              field.set_fill_value(np.NaN)
              zout=np.ma.array(field).mean(axis=axis,dtype=np.float64)
              if axis is not None:
                 zout.set_fill_value(np.NaN)
                 return zout.filled()
              else:return zout
           elif (np.isnan(np.sum(field)) and (type(field).__name__ not in 'MaskedArray')):
              zout=np.ma.masked_invalid(field).mean(axis=axis,dtype=np.float64)
              if axis is not None:
                 zout.set_fill_value([np.NaN])
                 return zout.filled()
              else:return zout
           else: 
              return np.array(field).mean(axis=axis,dtype=np.float64)

## compute sum
## author AS + AC
def sum(field,axis=None):
        if field is None: return None
        else:
           if type(field).__name__=='MaskedArray':
              field.set_fill_value(np.NaN)
              zout=np.ma.array(field).sum(axis=axis)
              if axis is not None:
                 zout.set_fill_value(np.NaN)
                 return zout.filled()
              else:return zout
           elif (np.isnan(np.sum(field)) and (type(field).__name__ not in 'MaskedArray')):
              zout=np.ma.masked_invalid(field).sum(axis=axis)
              if axis is not None:
                 zout.set_fill_value([np.NaN])
                 return zout.filled()
              else:return zout
           else:
              return np.array(field).sum(axis=axis)

## compute mean over bins
## author AS
def meanbin(y,x,bins):
    bins.sort() # sorting is needed for binning
    meanvalues = []
    for iii in range(len(bins)):
        ## GET VALUES OVER RELEVANT INTERVALS
        if iii == 0:
          ind = x<bins[iii+1] ; ys = y[ind]
        elif iii == len(bins)-1:
          ind = x>=bins[iii-1] ; ys = y[ind]
        else:
          ind = x>=bins[iii-1] ; interm = x[ind] ; intermf = y[ind]
          ind = interm<bins[iii+1] ; xs = interm[ind] ; ys = intermf[ind]
        ## COMPUTE MEAN and TREAT NAN CASE
        meanvalues.append(mean(ys))
    ## RETURN A NUMPY ARRAY
    meanvalues = np.array(meanvalues)
    return meanvalues

## compute rolling mean
## author DB
## --
## sfield,scoord = rollingmean(field,coord)
## --
def rollingmean(field,coord,axis=None,n=None):
    ret = np.cumsum(field, axis=axis, dtype=float)
    rec = np.cumsum(coord, dtype=float)
    rec[n:] = rec[n:] - rec[:-n]
    if axis == 0:
       ret[n:,:,:,:] = ret[n:,:,:,:] - ret[:-n,:,:,:]
       return ret[n:,:,:,:] / n, rec[n:] / n
    if axis == 1:
       ret[:,n:,:,:] = ret[:,n:,:,:] - ret[:,:-n,:,:]
       return ret[:,n:,:,:] / n, rec[n:] / n
    if axis == 2:
       ret[:,:,n:,:] = ret[:,:,n:,:] - ret[:,:,:-n,:]
       return ret[:,:,n:,:] / n, rec[n:] / n
    if axis == 3:
       ret[:,:,:,n:] = ret[:,:,:,n:] - ret[:,:,:,:-n]
       return ret[:,:,:,n:] / n, rec[n:] / n

## compute perturbation to mean
## -- the whole dimension must exist!
## author AS
def perturbation(field,axis=None,mm=None):
    # calculate mean (averaged dim is reduced)
    if mm is None:
      mm = mean(field,axis=axis)
    else:
      if mm.ndim != 3: print "input mean must be a 3-dim array" ; exit() 
      else: print "mean to compute pert is provided!"
    # include back the reduced dim
    if field.ndim == 4:
      if axis == 0:   mm = mm[np.newaxis,:,:,:]
      elif axis == 1: mm = mm[:,np.newaxis,:,:]
      elif axis == 2: mm = mm[:,:,np.newaxis,:]
      elif axis == 3: mm = mm[:,:,:,np.newaxis]
    else:
      print "not supported. yet!"
      exit()
    # repeat the calculated mean on this dimension
    mm = np.repeat(mm,field.shape[axis],axis=axis)
    # compute perturbations
    field = field - mm
    return field

## get leading exponent
## author AS
def get_exponent(value):
    return np.floor(np.log10(value))

## get normalizing factor corr. to leading exponent
## author AS
def get_norm(value):
    return 10.**get_exponent(value)

################
#### SMOOTH ####
################
### TBD: works with missing values

def smooth2diter(field,n=1):
        count = 0
        result = field
        while count < n:
            result = smooth2d(result, window=2)
            count = count + 1
        return result

## Author: AS. uses gauss_kern and blur_image.
def smooth2d(field, window=10):
	## actually blur_image could work with different coeff on x and y
        if True in np.isnan(field):
            print "!! ERROR !! Smooth is a disaster with missing values. This will be fixed."
            exit()
	if window > 1:	result = blur_image(field,int(window))
	else:		result = field
	return result

## FROM COOKBOOK http://www.scipy.org/Cookbook/SignalSmooth
def gauss_kern(size, sizey=None):
    	# Returns a normalized 2D gauss kernel array for convolutions
    	size = int(size)
    	if not sizey:
	        sizey = size
	else:
	        sizey = int(sizey)
	x, y = np.mgrid[-size:size+1, -sizey:sizey+1]
	g = np.exp(-(x**2/float(size)+y**2/float(sizey)))
	return g / g.sum()

## FROM COOKBOOK http://www.scipy.org/Cookbook/SignalSmooth
def blur_image(im, n, ny=None) :
	# blurs the image by convolving with a gaussian kernel of typical size n. 
	# The optional keyword argument ny allows for a different size in the y direction.
    	g = gauss_kern(n, sizey=ny)
	try:
          import scipy.signal as sp_signal
	  improc = sp_signal.convolve(im, g, mode='same')
	except:
	  print "blur_image needs scipy. please add scipy lib." ; exit()
    	return improc

## FROM COOKBOOK http://www.scipy.org/Cookbook/SignalSmooth
def smooth1d(vec,window=11,window_type='hanning'):
    """smooth the data using a window with requested size.
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    input:
        x: the input signal 
        window: the dimension of the smoothing window; should be an odd integer
        window_type: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.
    output:
        the smoothed signal
    example:
    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    see also: 
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
    TODO: the window parameter could be the window itself if an array instead of a string   
    """
    y = None ; w = None ; s = None ; x = None
    if True in np.isnan(vec):
        print "!! ERROR !! Smooth is a disaster with missing values. This will be fixed."
        exit()    
    x = np.array(vec)
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
    if x.size < window:
        print x.size,window
        raise ValueError, "Input vector needs to be bigger than window size."
    if window<3:
        return x
    if not window_type in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
    s=np.r_[x[window-1:0:-1],x,x[-1:-window:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window,'d')
    else:
        w=eval('np.'+window_type+'(window)')
    y=np.convolve(w/w.sum(),s,mode='valid')
    return y

########################
#### TIME CONVERTER ####
########################

#  mars_sol2ls
#  author T. Navarro
#  -----------------
#  convert a given martian day number (sol)
#  into corresponding solar longitude, Ls (in degr.),
#  where sol=0=Ls=0 is the
#  northern hemisphere spring equinox.
#  -----------------
def mars_sol2ls(soltabin,forcecontinuity=True):
      year_day  = 668.6
      peri_day  = 485.35
      e_elips   = 0.09340
      radtodeg  = 57.2957795130823
      timeperi  = 1.90258341759902
      if type(soltabin).__name__ in ['int','float','float32','float64']:
        soltab=[soltabin]
        solout=np.zeros([1])
      else:
        soltab=soltabin
        solout=np.zeros([len(soltab)])
      i=0
      for sol in soltab:
         zz=(sol-peri_day)/year_day
         zanom=2.*np.pi*(zz-np.floor(zz))
         xref=np.abs(zanom)
#  The equation zx0 - e * sin (zx0) = xref, solved by Newton
         zx0=xref+e_elips*m.sin(xref)
         iter=0
         while iter <= 10:
            iter=iter+1
            zdx=-(zx0-e_elips*m.sin(zx0)-xref)/(1.-e_elips*m.cos(zx0))
            if(np.abs(zdx) <= (1.e-7)):
              continue
            zx0=zx0+zdx
         zx0=zx0+zdx
         if(zanom < 0.): zx0=-zx0
# compute true anomaly zteta, now that eccentric anomaly zx0 is known
         zteta=2.*m.atan(m.sqrt((1.+e_elips)/(1.-e_elips))*m.tan(zx0/2.))
# compute Ls
         ls=zteta-timeperi
         if(ls < 0.): ls=ls+2.*np.pi
         if(ls > 2.*np.pi): ls=ls-2.*np.pi
# convert Ls in deg.
         ls=radtodeg*ls
         solout[i]=ls
         i=i+1
      if forcecontinuity:
         for iii in range(len(soltab)-1):
           while solout[iii+1] - solout[iii] > 180. :  solout[iii+1] = solout[iii+1] - 360.
           while solout[iii] - solout[iii+1] > 180. :  solout[iii+1] = solout[iii+1] + 360.
      return solout

# mars_date
# author A. Spiga
# ------------------------------
# get Ls, sol, utc from a string with format yyyy-mm-dd_hh:00:00
# -- argument timechar is a vector of such strings indicating dates
# -- example: timechar = nc.variables['Times'][:] in mesoscale files
# NB: uses mars_sol2ls function above
# ------------------------------
def mars_date(timechar):
    # some preliminary information
    days_in_month = [61, 66, 66, 65, 60, 54, 50, 46, 47, 47, 51, 56]
    plus_in_month = [ 0, 61,127,193,258,318,372,422,468,515,562,613]
    # get utc and sol from strings
    utc = [] ; sol = []
    for zetime in timechar:
        dautc = float(zetime[11]+zetime[12]) + float(zetime[14]+zetime[15])/37.
        dasol = dautc / 24.
        dasol = dasol + plus_in_month[ int(zetime[5]+zetime[6])-1 ] + int(zetime[8]+zetime[9])
        dasol = dasol - 1 ##les sols GCM commencent a 0
        utc.append(dautc)
        sol.append(dasol)
    sol = np.array(sol)
    utc = np.array(utc)
    # get ls from sol
    ls = mars_sol2ls(sol)
    return ls, sol, utc

# timecorrect
# author A. Spiga
# -----------------------------
# ensure time axis is monotonic
# correct negative values
# -----------------------------
def timecorrect(time):
    for ind in range(len(time)-1):
        if time[ind] < 0.: time[ind] = time[ind] + 24.
        if time[ind+1] < time[ind]: time[ind+1] = time[ind+1] + 24.
    return time

#########################
#### FLOW DIAGNOSTIC ####
#########################

def deriv1d(field,coord=None): 
  ## dfdx = deriv(f,x)
  df = np.gradient(field)
  if coord is not None:
    dx = np.gradient(coord,edge_order=2)
    return df/dx
  else:
    return df

def dxdy(lon,lat,coeff=None,lonlat=False):
  ####################
  # fill 2D array w/ derivatives of two coordinates
  ####################
  # lon: coord1
  # lat: coord2
  # coeff: e.g. pi/180 for deg --> rad
  # lonlat: True if lon is longitude and lat is latitude
  ####################
  if coeff is not None:
    # convert factor (e.g. deg to rad)
    la = lon*coeff
    ph = lat*coeff
  else:
    la = lon
    ph = lat
  # ensure 2D arrays for lon,lat
  if lon.ndim == 1 and lat.ndim == 1:
    [lar,phr] = np.meshgrid(la,ph)
  elif lon.ndim == 2 and lat.ndim == 2:
    lar,phr = la,ph
  else:
    print "ppcompute. there is a problem with coordinates rank. stop."
    exit()
  # compute normalized gradients for lat/lon grid
  # -- cartesian differential coordinates
  dump,dx = np.gradient(lar,edge_order=2)
  dy,dump = np.gradient(phr,edge_order=2)
  # treat the particular case of lon/lat
  if lonlat:
    dx = dx*np.cos(phr)
  return dx,dy
 
def deriv2d(u,coordx,coordy,coeff=None,lonlat=False,fac=1.):
  ###
  ### equivalent to  
  ###   x,y = np.meshgrid(coordx,coordy)
  ###   dx,dy = dxdy(coordx,coordy)
  ###   du_dy,du_dx = np.gradient(u,dy,dx,edge_order=2)
  ### though slightly more flexible
  ###
  # compute gradients -- with normalized coordinates
  du_y,du_x = np.gradient(u,edge_order=2)
  # compute cartesian increments
  dx,dy = dxdy(coordx,coordy,coeff=coeff,lonlat=lonlat)
  # eventually compute cartesian derivatives
  du_dx = du_x/(fac*dx)
  du_dy = du_y/(fac*dy)
  # output
  return du_dx, du_dy 

def divort(u,v,lon,lat,rad):
  ####################
  # compute divergence and vorticity
  # div,vorti = divort(u,v,lon,lat,rad)
  ####################
  # u: zonal wind
  # v: meridional wind
  # lon: longitude
  # lat: latitude
  # rad: radius of the planet
  ####################
  du_dx,du_dy = deriv2d(u,lon,lat,coeff=np.pi/180.,lonlat=True,fac=rad)
  dv_dx,dv_dy = deriv2d(v,lon,lat,coeff=np.pi/180.,lonlat=True,fac=rad)
  vorti = dv_dx - du_dy
  div = du_dx + dv_dy 
  return div,vorti

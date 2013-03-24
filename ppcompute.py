###############################################
## PLANETOPLOT                               ##
## --> PPCOMPUTE                             ##
###############################################
# python built-in librairies
import os
# added librairies
import numpy as np
import scipy.signal as sp_signal
###############################################

## first a useful function to find settings in a folder in PYTHONPATH
def findset(whereset,string="planetoplot_v2"):
    # ... set a default whereset if it was set to None
    # ... default is in the planetoplot_v2 folder
    if whereset is None:
        for path in os.environ['PYTHONPATH'].split(os.pathsep):
            if string in path: whereset = path
        if whereset is None: 
            print "!! ERROR !! "+ string + "not in $PYTHONPATH"
            print "--> either put it in $PYTHONPATH or change whereset"
            exit()
    # ... if the last / is missing put it
    if whereset[-1] != "/": whereset = whereset + "/"
    return whereset

##########################
#### MAX MEAN MIN SUM ####
#####################################
#### WITH SUPPORT FOR NaN VALUES ####
#####################################

## compute min
## author AS + AC
def min (field,axis=None): 
        if field is None: return None
        if type(field).__name__=='MaskedArray':
              field.set_fill_value(np.NaN)
              return np.ma.array(field).min(axis=axis)
        elif (np.isnan(np.sum(field)) and (type(field).__name__ not in 'MaskedArray')):
              return np.ma.masked_invalid(field).min(axis=axis)
        else: return np.array(field).min(axis=axis)

## compute max
## author AS + AC
def max (field,axis=None):
        if field is None: return None
        if type(field).__name__=='MaskedArray':
              field.set_fill_value(np.NaN)
              return np.ma.array(field).max(axis=axis)
        elif (np.isnan(np.sum(field)) and (type(field).__name__ not in 'MaskedArray')):
              return np.ma.masked_invalid(field).max(axis=axis) 
        else: return np.array(field).max(axis=axis)

## compute mean
## author AS + AC
def mean (field,axis=None):
        if field is None: return None
        else: 
           if type(field).__name__=='MaskedArray':
              field.set_fill_value(np.NaN)
              zout=np.ma.array(field).mean(axis=axis)
              if axis is not None:
                 zout.set_fill_value(np.NaN)
                 return zout.filled()
              else:return zout
           elif (np.isnan(np.sum(field)) and (type(field).__name__ not in 'MaskedArray')):
              zout=np.ma.masked_invalid(field).mean(axis=axis)
              if axis is not None:
                 zout.set_fill_value([np.NaN])
                 return zout.filled()
              else:return zout
           else: 
              return np.array(field).mean(axis=axis)

## compute sum
## author AS + AC
def sum (field,axis=None):
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

################
#### SMOOTH ####
################
### TBD: works with missing values

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
    	improc = sp_signal.convolve(im, g, mode='same')
    	return improc

## FROM COOKBOOK http://www.scipy.org/Cookbook/SignalSmooth
def smooth1d(x,window=11,window_type='hanning'):
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
    if True in np.isnan(field):
        print "!! ERROR !! Smooth is a disaster with missing values. This will be fixed."
        exit()    
    x = np.array(x)
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
    if x.size < window:
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

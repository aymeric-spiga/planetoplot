#! /usr/bin/env python

import pickle
import planets
from ppclass import pp
import numpy as np
import ppplot
from scipy import fftpack
import sys
from optparse import OptionParser ### TBR by argparse
parser = OptionParser()
parser.add_option('-f','--file',action='store',dest='file',type="string",default="extract.nc",help="file")
parser.add_option('-v','--var',action='store',dest='var',type="string",default="ps",help="var")
parser.add_option('-y','--lat',action='store',dest='y',type="string",default="0.",help="y (lat)")
parser.add_option('-z','--vert',action='store',dest='z',type="string",default="0.",help="z (vert)")
parser.add_option('-d','--dt',action='store',dest='dt',type="float",default=None,help="data points each ** time unit (year)")
parser.add_option('-s','--sigma',action='store',dest='sigma',type="float",default=None,help="how many sigma around mean for colorbar")
parser.add_option('-c','--colorbar',action='store',dest='cb',type='string',default="CMRmap",help="colorbar")
parser.add_option('-o','--output',action='store',dest='output',type='string',default=None,help="name of png output (gui if None)")
parser.add_option('--noreldis',action='store_true',dest='noreldis',default=False,help="no relation dispersion plotted")
(opt,args) = parser.parse_args()

##aymeric@aymeric-laptop:~/TMPDIR/fft$ testsat.py -v v -f yorgl.nc -y 20. -z 490
##aymeric@aymeric-laptop:~/TMPDIR/fft$ testsat.py -v v -f yorgl.nc -y 3.75 -z 118

## testsat.py -f ~/Remote/saturn_128x96x64_guided.nc -y 0. -z 5e4 -v temp
## testsat.py -f ~/Remote/saturn_128x96x64_guided.nc -y 0. -z 5e4 -v u
## testsat.py -f ~/Remote/saturn_128x96x64_guided.nc -y 0. -z 5e4 -v v

## testsat.py -y 0. -z 5e4 -v temp -f ~/Remote/saturn_0.5deg_dynamico.nc -d 0.0008 -s 15

###################################
##
## 1. SPECTRA FROM SIMULATIONS
##
###################################

## FIELD
tab = pp(file=opt.file,var=opt.var,y=opt.y,z=opt.z,verbose=True).getf()

## MAKE X=LON Y=TIME
tab = np.transpose(tab)

## FIELD SIZES
nx,ny = tab.shape
limxmax = (nx-1)/4
limxmin = -limxmax

## SAMPLE SPACES
dx = 1./nx   # data points each dx planet --> result in zonal wavenumber
if opt.dt is None:
 dy = 1./250. # data points each dy time units (e.g. sols) --> result in (time units)^-1
else:
 dy = opt.dt

## 2D FFT
## http://docs.scipy.org/doc/scipy/reference/fftpack.html
spec = fftpack.fft2(tab)

## FREQUENCY COORDINATES
specx = fftpack.fftfreq(nx,d=dx)
specy = fftpack.fftfreq(ny,d=dy)

## REMOVE DC COMPONENT 
spec[:,0] = 0. #np.nan
spec[0,:] = 0.

## MOVE 0 TO CENTER WITH NEG --> POS
spec = fftpack.fftshift(spec)
specx = fftpack.fftshift(specx)
specy = fftpack.fftshift(specy)
#right westward/eastward orientation
specx = specx[::-1]

## POWER SPECTRA
spec = np.abs(spec)**2
#spec = np.log10(spec)

## PLOT
if opt.sigma is not None:
  ppplot.how_many_sigma = opt.sigma  # the 3sig max-min is not wide enough
                                     # -- increase this to make max/min more contrasted
ppplot.cline = 0.35
ppplot.negative(howblack="0.05")

p = ppplot.plot2d()
p.f = spec
p.x = specx
p.y = specy
p.colorbar = opt.cb
p.axisbg = 0.75
p.div = 50
p.xmin = limxmin
p.xmax = limxmax
p.nxticks = 8
#p.nyticks = 1./dy/2.
p.ylabel = r"frequency $\sigma$ (year$^{-1}$)"
p.xlabel = r"$\leftarrow$ Westward | wavenumber $s$ | Eastward $\rightarrow$"
p.showcb = False
p.logy = True
p.ymax = 1000. #10.*25000.*dy
p.title = "Zonal-Time FFT for "+opt.var
if (not opt.noreldis):
  p.title = p.title + " + RW / KW dispersion relations"
else:
  foo = np.logspace(-3.,-1.,10)*25000.
  p.ymin = np.min(foo)
p.make()

###################################
##
## 2. DISPERSION RELATIONSHIP
##
###################################

if (not opt.noreldis):

  ####################################
  lz = 60000. # vue dans la simu
  ####################################
  nutab = [+1,+2,+3]
  nutab = [+1,+2,+3,+4,+5]
  nutab = [-2,-1,+1,+2]
  ####################################
  #T0=pp(file=opt.file,var="temp",y=opt.y,z=opt.z,x=0,t="0,10000").getf()
  ##--environ 120K, OK.
  ###################################
  n2tab = []
  #n2tab.append(0.3e-5) # SL nat2011
  n2tab.append(1.0e-5) # LL grl2008
  #n2tab.append(planets.Saturn.N2())
  #n2tab.append(planets.Saturn.N2(dTdz=-0.7e-3)) # proche LL
  ####################################

  specx = np.linspace(limxmin,limxmax,100)
  specy = np.logspace(-3.,-1.,100)
  s,sigma = np.meshgrid(specx,specy)

  p.f = sigma*np.nan # trick to make it transparent
  p.clev = [0.]
  p.ccol = "white"
  p.x = specx
  p.y = 25000.*specy
  p.invert = True
  for nnn in nutab:
   for n2n2 in n2tab:
     p.c = planets.Saturn.dispeqw(s,sigma,nu=nnn,lz=lz,N2=n2n2)
     p.make()

if opt.output is None:
  ppplot.show()
else:
  ppplot.save(mode="png",filename=opt.output)

####################################
# save a .sh file with the command #
####################################
command = ""
for arg in sys.argv: command = command + arg + ' '
if opt.output is not None:
  try:
    f = open(opt.output+'.sh', 'w')
    f.write(command)
  except IOError:
    print "!! WARNING !! not saved. Probably do not have permission to write here."




#! /usr/bin/env python

################
################
## SPECTRA.PY ##
################
################

## examples
## ** mars' thermal tides in a file with outputs every 0.0833 sols
##     ./spectra.py ps100.nc -v ps -d 0.0833 -u sol -y 20 --sigma 10 --ymin 0.5 --ymax 8 --xmin -8 --xmax 8 -o mars_tides
## ** saturn waves in a file with outputs every 5 sols in a year of ~25000 sols
##     ./spectra.py outremap.nc  --dt 0.0002 --unit year -y 45 --sigma 5 --logy

## modules
import planets
from ppclass import pp
import numpy as np
import ppplot
from scipy import fftpack
import sys

##############
## PREAMBLE ##
##############

## script options
from optparse import OptionParser ### TBR by argparse
parser = OptionParser()
parser.add_option('-v','--var',action='store',dest='var',type="string",default="ps",help="var (default: ps)")
parser.add_option('-y','--lat',action='store',dest='y',type="string",default="0.",help="y (lat, default: 0)")
parser.add_option('-z','--vert',action='store',dest='z',type="string",default="0.",help="z (vert, default: 0)")
parser.add_option('-u','--unit',action='store',dest='unit',type='string',default="time unit",help="time unit (spectra in UNIT^-1, default: time unit)")
parser.add_option('-d','--dt',action='store',dest='dt',type="float",default=None,help="in FILE, one data point each time UNIT (default: 1)")
parser.add_option('-o','--output',action='store',dest='output',type='string',default=None,help="name of png output (gui if None)")
parser.add_option('--reldis',action='store_true',dest='reldis',default=False,help="add dispersion relationship")
parser.add_option('--log',action='store_true',dest='log',default=False,help="set log field")

## get planetoplot-like options
parser = ppplot.opt(parser) # common options for plots
parser = ppplot.opt2d(parser) # common options for plots
(opt,args) = parser.parse_args()

## input file
infile = args[0]

## customize default behaviour for this script
if opt.colorbar is None: 
  opt.colorbar = "CMRmap"
else: 
  opt.colorbar = opt.colorbar[0]
#
if opt.ylabel is None: 
  opt.ylabel = r"frequency $\sigma$ (cycles per "+opt.unit+")"
#
if opt.xlabel is None: 
  opt.xlabel = r"$\leftarrow$ Westward | wavenumber $s$ | Eastward $\rightarrow$"
#
if opt.title is None: 
  opt.title = "2D Zonal-Time FFT for "+opt.var
#
if (opt.reldis): 
  opt.title = opt.title + " + RW / KW dispersion relations"
else: 
  opt.title = opt.title + " at latitude "+str(opt.y)

## hardwired settings for this script
opt.div = 50
opt.showcb = False
opt.back = 'black'
ppplot.cline = 0.35
#ppplot.negative(howblack="0.05")

## a few specific fixes
if opt.ymin is not None: opt.ymin = opt.ymin[0]
if opt.ymax is not None: opt.ymax = opt.ymax[0]
if opt.xmin is not None: opt.xmin = opt.xmin[0]
if opt.xmax is not None: opt.xmax = opt.xmax[0]
if opt.logy and (opt.ymin is None): opt.ymin = 1e-6
if opt.logy and (opt.ymax is None): opt.ymax = 1e6
if opt.ymin is None: opt.ymin = 0 # remove symmetric negative frequency

##############################
## SPECTRA FROM SIMULATIONS ##
##############################

## FIELD
tab,x,y,z,t = pp(file=infile,var=opt.var,y=opt.y,z=opt.z,verbose=True).getfd()

## MAKE X=LON Y=TIME
tab = np.transpose(tab)

## FIELD SIZES
nx,nt = tab.shape

### TREAT PROBLEM of EVEN NUMBER OF LONGITUDES (which is a problem with fftshift) 
### ... ADD A MIRROR POINT IN THE END
if (nx % 2 == 0):
  tab2 = np.zeros([nx+1,nt])
  tab2[0:nx-1,:] = tab[0:nx-1,:]
  tab2[nx,:] = tab[0,:]
  tab = tab2
  nx = nx + 1
     
## SAMPLE SPACES
# data points each dx planet --> result in zonal wavenumber
dx = 1./nx
# data points each opt.dt time units --> result in (time unit)^-1
if opt.dt is None: dt = 1.
else: dt = opt.dt

## PERFORM 2D FFT
## http://docs.scipy.org/doc/scipy/reference/fftpack.html
spec = fftpack.fft2(tab)

## FREQUENCY COORDINATES
specx = fftpack.fftfreq(nx,d=dx)
spect = fftpack.fftfreq(nt,d=dt)

## REMOVE DC COMPONENT 
if not opt.log:
  spec[:,0] = 0. #np.nan
  spec[0,:] = 0.

## MOVE 0 TO CENTER WITH NEG --> POS
spec = fftpack.fftshift(spec)
specx = fftpack.fftshift(specx)
spect = fftpack.fftshift(spect)

## GET THE RIGHT W/E ORIENTATION
if x[0] < x[-1]:
  #right westward/eastward orientation
  specx = specx[::-1]
else:
  #reverted axis (e.g. 180 --> -180, Venus): nothing to do
  pass

## POWER SPECTRA
spec = np.abs(spec)**2

## MAKE IT LOG
if opt.log: spec = np.log10(spec)

## HORIZONTAL WAVENUMBER LIMITS
limxmax = (nx-1)/4
limxmin = -limxmax

## PLOT
p = ppplot.plot2d()
p.transopt(opt) # transfer plot options
p.f = spec
p.x = specx
p.y = spect
p.make()

##################################################################################################
##################################################################################################
##################################################################################################

##aymeric@aymeric-laptop:~/TMPDIR/fft$ testsat.py -v v -f yorgl.nc -y 20. -z 490
##aymeric@aymeric-laptop:~/TMPDIR/fft$ testsat.py -v v -f yorgl.nc -y 3.75 -z 118
## testsat.py -f ~/Remote/saturn_128x96x64_guided.nc -y 0. -z 5e4 -v temp
## testsat.py -f ~/Remote/saturn_128x96x64_guided.nc -y 0. -z 5e4 -v u
## testsat.py -f ~/Remote/saturn_128x96x64_guided.nc -y 0. -z 5e4 -v v
## testsat.py -y 0. -z 5e4 -v temp -f ~/Remote/saturn_0.5deg_dtnamico.nc -d 0.0008 -s 15 --reldis --logy
## spectra.py -f /home/aymeric/Big_Data/DATAPLOT/diagfired.nc -v ps -d 0.00015 

###################################
##
## 2. DISPERSION RELATIONSHIP
##
###################################

if (opt.reldis):

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
  spect = np.logspace(-3.,-1.,100)
  s,sigma = np.meshgrid(specx,spect)

  p.f = sigma*np.nan # trick to make it transparent
  p.clev = [0.]
  p.ccol = "white"
  p.ccol = "red"
  p.x = specx
  p.y = 25000.*spect
  p.invert = True
  for nnn in nutab:
   for n2n2 in n2tab:
     p.c = planets.Saturn.dispeqw(s,sigma,nu=nnn,lz=lz,N2=n2n2)
     p.make()

if opt.output is None: ppplot.show()
else: ppplot.save(mode="png",filename=opt.output)

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




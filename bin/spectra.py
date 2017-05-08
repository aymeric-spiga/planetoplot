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
from ppclass import pp
import numpy as np
import ppplot
from scipy import fftpack
from scipy.ndimage.measurements import maximum_position
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
parser.add_option('-d','--dt',action='store',dest='dt',type="float",default=1.,help="in FILE, one data point each time UNIT (default: 1)")
parser.add_option('-o','--output',action='store',dest='output',type='string',default=None,help="name of png output (gui if None)")
parser.add_option('--reldis',action='store_true',dest='reldis',default=False,help="add dispersion relationship")
parser.add_option('--log',action='store_true',dest='log',default=False,help="set log field")
parser.add_option('--period',action='store_true',dest='period',default=False,help="show period (in UNIT) instead of frequency")
parser.add_option('--ndom',action='store',dest='ndom',type="int",default=10,help="print info for the NDOM dominant modes (default:10)")
parser.add_option('--noplot',action='store_true',dest='noplot',default=False,help="do not plot anything, just output dominant modes")

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
 if not opt.period:
  opt.ylabel = r"frequency $\sigma$ ($^{\circ}$ "+opt.unit+"$^{-1}$)"
 else:
  opt.ylabel = r"period ("+opt.unit+")"
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
if opt.logy and (opt.ymin is None): opt.ymin = 1e-6
if opt.logy and (opt.ymax is None): opt.ymax = 1e6
if opt.ymin is None: 
 if not opt.period:
  opt.ymin = 0 # remove symmetric negative frequency
 else:
  opt.ymin = 2.*opt.dt # minimum possible period

##############################
## SPECTRA FROM SIMULATIONS ##
##############################

## FIELD
vb = False
tab,x,y,z,t = pp(file=infile,var=opt.var,y=opt.y,z=opt.z,verbose=vb).getfd() # pert_x ne change rien

## symm / anti-symm components, see Wheeler 1999
#tabtropN = pp(file=infile,var=opt.var,y=+5.,z=opt.z,verbose=vb).getf()
#tabtropS = pp(file=infile,var=opt.var,y=-5.,z=opt.z,verbose=vb).getf()
#tab = 0.5*tabtropN + 0.5*tabtropS # symmetric
##tab = 0.5*tabtropN - 0.5*tabtropS # antisymmetric

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
lowerperiod = 3.*opt.dt # Nyquist rate + 1
higherperiod = opt.dt*float(nt-1)/2. # half size of sample

## PERFORM 2D FFT
## http://docs.scipy.org/doc/scipy/reference/fftpack.html
spec = fftpack.fft2(tab)

## FREQUENCY COORDINATES
specx = fftpack.fftfreq(nx,d=dx)
spect = fftpack.fftfreq(nt,d=opt.dt)

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
if opt.xmin is None: opt.xmin = limxmin
if opt.xmax is None: opt.xmax = limxmax

## RETAIN ONLY POSITIVE FREQUENCIES
## -- (and remove period=sample_size)
mm = min(spect[spect>0])
w = spect > mm
spect = spect[w]
spec = spec[:,w]

## SEARCH FOR DOMINANT MODES
# -- initialize output
if opt.output is not None:
  txtfile = opt.output+".txt"
else:
  txtfile = "spectra.txt"
fifi = open(txtfile, "w")
fifi.write("---------------------------------------\n")
fifi.write("%2s %4s %8s %8s %8s\n" % ("n","WN","dg/"+opt.unit,opt.unit,"log(A)"))
fifi.write("---------------------------------------\n")
# -- initialize while loop
search = np.empty_like(spec) ; search[:,:] = spec[:,:]
zelab = search > 0 # (all elements)
itit = 1 
# -- while loop
while itit <= opt.ndom:
  # -- find dominant mode
  ij = maximum_position(search,labels=zelab)
  dominant_wn = specx[ij[0]]
  dominant_fq = spect[ij[1]]
  spower = search[ij]
  if (1./dominant_fq) < lowerperiod: reliable = "x"
  else: reliable = "o"
  # -- print result
  if reliable == "o":
    fifi.write("%2i %4.0f %8.1f %8.1f %8.1f %4s\n" % (itit,dominant_wn,360.*dominant_fq,1./dominant_fq,np.log10(spower),reliable))
  # -- iterate
  zelab = search < spower # remove maximum found
  search[ij[0],:] = -9999. # remove wavenumber found (otherwise loop could find other maxima for this wn)
  itit += 1
# -- close output
fifi.write("---------------------------------------\n")
fifi.close()
# -- print results
print(open(txtfile, "r").read())

## COMPUTE FREQUENCY/PERIOD AXIS
lowerperiod = 2.5*opt.dt # Nyquist rate + 1
higherperiod = opt.dt*float(nt-1)/2. # half size of sample
if not opt.period:
  # frequency: longitude degree per UNIT
  spect = 360.*spect
  if opt.ymax is None: opt.ymax = 360./lowerperiod
  elif opt.ymin is None: opt.ymin = 360./higherperiod
else:
  # period: UNIT
  spect = 1./(spect)
  if opt.ymin is None: opt.ymin = lowerperiod
  elif opt.ymax is None: opt.ymax = higherperiod

## PLOT
if not opt.noplot:
  p = ppplot.plot2d()
  p.transopt(opt) # transfer plot options
  p.f = spec
  p.x = specx
  p.y = spect
  p.make()

###################################
##
## 2. DISPERSION RELATIONSHIP
##
###################################

if (opt.reldis):

  ## planet
  try: 
    import planets
  except: 
    print "please install the module planets"
  mypl = planets.Saturn
  #mypl = planets.Jupiter

  ####################################
  lz = 60000. # vue dans la simu?
  lz = 2.*mypl.H() # a kind of generic choice
  ####################################
  nutab = [+1,+2,+3]
  nutab = [+1,+2,+3,+4,+5]
  nutab = [-3,-2,-1,+1,+2,+3]
  #nutab = [-5,-4,-3,-2,-1,0,+1,+2,+3,+4,+5]
  #nutab = [-2,-1,0,+1,+2]
  nutab = [-1,0,+1]
  #lz = 3000.
  ####################################
  #T0=pp(file=opt.file,var="temp",y=opt.y,z=opt.z,x=0,t="0,10000").getf()
  ##--environ 120K, OK.
  ###################################
  n2tab = []
  #n2tab.append(0.3e-5) # SL nat2011
  n2tab.append(1.0e-5) # LL grl2008
  #n2tab.append(mypl.N2()) # simple
  #n2tab.append(mypl.N2(dTdz=-0.7e-3)) # proche LL
  ####################################
  hache = 100.
  #hache = 10.
  #hache = 1.

  # ensure number of points 
  # -- is enough for smooth lines
  # -- is not too high for efficiency
  n = 100
  specx = np.linspace(limxmin,limxmax,2*n)
  spect = np.linspace(spect.min(),spect.max(),n)
  spect = spect / 360. # convert back from deglon/unit to cycle/unit
  s,sigma = np.meshgrid(specx,spect)

  # a few general settings for plots
  p.f = sigma*np.nan # trick to make it transparent
  p.clev = [0.]
  p.clab = False
  p.x = specx

  ## COMPUTE FREQUENCY/PERIOD AXIS
  if not opt.period:
    # frequency: longitude degree per UNIT
    p.y = 360.*spect
  else:
    # period: UNIT
    p.y = 1./(spect)

  ## COMPUTE dispersion relationship for all modes
  for nnn in nutab:
   for n2n2 in n2tab:
     if nnn == 0: p.ccol = "cyan"
     elif nnn > 0: p.ccol = "magenta"
     else: p.ccol = "red"
     #p.c = mypl.dispeqw(s,sigma,nu=nnn,lz=lz,N2=n2n2)
     p.c = mypl.dispeqw(s,sigma,nu=nnn,h=hache)
     p.make()

### SHOW or SAVE PLOT
if not opt.noplot:
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

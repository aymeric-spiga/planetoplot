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
parser.add_option('-k','--kind',action='store',dest='kind',type='string',default="local",help="kind of diagnostic: 'local' at lat [D], 'sym' or 'antisym' at +/- lat")
parser.add_option('-e','--height',action='store',dest='height',type='float',default=10.,help="equivalent height [m, default: 10m]")
parser.add_option('--savef',action='store',dest='savef',type="string",default="png",help="choose save format of outputs")
parser.add_option('--reldis',action='store_true',dest='reldis',default=False,help="add dispersion relationship")
parser.add_option('--log',action='store_true',dest='log',default=False,help="set log field")
parser.add_option('--period',action='store_true',dest='period',default=False,help="show period (in UNIT) instead of frequency")
parser.add_option('--ndom',action='store',dest='ndom',type="int",default=10,help="print info for the NDOM dominant modes (default:10)")
parser.add_option('--noplot',action='store_true',dest='noplot',default=False,help="do not plot anything, just output dominant modes")
parser.add_option('--nutab',action='store',dest='nutab',type="int",default=1,help="dispersion relation to include (0,1,2 ; default:1)")
parser.add_option('--cycle',action='store_true',dest='cycle',default=False,help="use unit cycle instead of degree for longitude")
parser.add_option('--tmin',action='store',dest='tmin',default=None,type="float",help="minimum value on time axis")
parser.add_option('--tmax',action='store',dest='tmax',default=None,type="float",help="maximum value on time axis")
parser.add_option('--zonalwind',action='store',dest='zonalwind',default=None,type="float",help="zonal wind to compute intrinsic frequency (default is None)")

## get planetoplot-like options
parser = ppplot.opt(parser) # common options for plots
parser = ppplot.opt2d(parser) # common options for plots
(opt,args) = parser.parse_args()

## input file
infile = args[0]

## customize default behaviour for this script
if opt.cycle:
  fac = 1.
  unitlon = "cy"
else:
  fac = 360.
  unitlon = "^{\circ}"
#
if opt.colorbar is None: 
  opt.colorbar = "CMRmap"
else: 
  opt.colorbar = opt.colorbar[0]
#
if opt.ylabel is None: 
 if not opt.period:
  opt.ylabel = r"frequency $\sigma$ ($"+unitlon+"$ "+opt.unit+"$^{-1}$)"
 else:
  opt.ylabel = r"period ("+opt.unit+")"
#
if opt.xlabel is None: 
  opt.xlabel = r"$\leftarrow$ Westward | wavenumber $s$ | Eastward $\rightarrow$"
#
if opt.title is None: 
  opt.title = "2D Zonal-Time FFT for "+opt.var
#
if opt.kind == "local":
  opt.title = opt.title + " at latitude "+str(opt.y)
elif opt.kind == "sym":
  opt.title = opt.title + " symmetric, lat +/-"+str(opt.y)
elif opt.kind == "antisym":
  opt.title = opt.title + " antisymmetric, lat +/-"+str(opt.y)

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

## function for Doppler shift
if opt.zonalwind is not None:
    print "ASSUMING ZONAL WIND:",opt.zonalwind
    import planets
    mypl = planets.Saturn
    # transform u in degree per day
    # 2*pi*radius/360 m <> 1 deg
    m_per_deg = 2.*np.pi*mypl.a / 360.
    s_per_day = mypl.day
    utransform = opt.zonalwind * s_per_day / m_per_deg
    ####
    def intrinsic(sigma,s):
        sigmahat = sigma - np.sign(s)*utransform
        if np.sign(sigmahat) == -1:
            print "reversal !",sigma,sigmahat
            shat = -s
            sigmahat = -sigmahat
        else:
            shat = s
        return sigmahat,shat
else:
    utransform = None

##############################
## SPECTRA FROM SIMULATIONS ##
##############################

## FIELD
vb = False
if opt.kind == "local":
   tab,x,y,z,t = pp(file=infile,var=opt.var,y=opt.y,z=opt.z,verbose=vb).getfd() # pert_x ne change rien
else:
   # symm / anti-symm components, see Wheeler & Kiladis 1999
   dalat = np.abs(np.float(opt.y))
   tabtropN,x,y,z,t = pp(file=infile,var=opt.var,y=+dalat,z=opt.z,verbose=vb).getfd()
   tabtropS         = pp(file=infile,var=opt.var,y=-dalat,z=opt.z,verbose=vb).getf()
   if opt.kind == "sym":
      tab = 0.5*tabtropN + 0.5*tabtropS # symmetric
   elif opt.kind == "antisym":
      tab = 0.5*tabtropN - 0.5*tabtropS # antisymmetric

## MIN and MAX VALUE on TIME
if opt.tmin is not None and opt.tmax is not None:
    ind = (t >= opt.tmin)*(t <= opt.tmax)
    tab = tab[ind]
    print t[ind]

## MAKE X=LON Y=TIME
tab = np.transpose(tab)

## FIELD SIZES
nx,nt = tab.shape

### TREAT PROBLEM of EVEN NUMBER OF LONGITUDES (which is a problem with fftshift) 
### ... ADD A MIRROR POINT IN THE END
#if (nx % 2 == 0):
#  tab2 = np.zeros([nx+1,nt])
#  tab2[0:nx-1,:] = tab[0:nx-1,:]
#  tab2[nx,:] = tab[0,:]
#  tab = tab2
#  nx = nx + 1
if (nx % 2 == 0):
  tab2 = np.zeros([nx-1,nt])
  tab2[0:nx-2,:] = tab[0:nx-2,:]
  tab = tab2
  nx = nx - 1     

## SAMPLE SPACES
# data points each dx planet --> result in zonal wavenumber
dx = 1./nx
# data points each opt.dt time units --> result in (time unit)^-1
lowerperiod = 4.*opt.dt # Nyquist rate
if opt.ymax is not None: 
  lowerperiod = fac/opt.ymax
higherperiod = opt.dt*float(nt-1)/2.2 # half size of sample (and 10% not to be at limit)

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

### RETAIN ONLY POSITIVE FREQUENCIES
### -- (and remove period=sample_size)
mm = min(spect[spect>0])
w = spect > mm
spect = spect[w]
#spec = spec[:,w] #before
### have to collect all power (complex --> real)
### ... avoid spectral leakage
### 90 deg per day semble etre la limite...
spec = spec[:,w] + spec[::-1,w[::-1]]
# https://gist.github.com/endolith/236567
# https://calebmadrigal.com/fourier-transform-notes/
# https://dsp.stackexchange.com/questions/3466/amplitude-of-the-signal-in-frequency-domain-different-from-time-domain
# https://dsp.stackexchange.com/questions/4825/why-is-the-fft-mirrored
# https://electronics.stackexchange.com/questions/12407/what-is-the-relation-between-fft-length-and-frequency-resolution

## SEARCH FOR DOMINANT MODES
# -- initialize output
if opt.output is not None:
  txtfile = opt.output+".txt"
else:
  txtfile = "spectra.txt"
fifi = open(txtfile, "w")
fifi.write(opt.title+"\n")
fifi.write("---------------------------------------\n")
if opt.zonalwind is not None:
    fifi.write("%4s & %8s & %8s & %8s \\\\ \hline \n" % ("$s$","$\hat{sigma} \, ("+unitlon+"$/"+opt.unit+")","intr. period ("+opt.unit+")","log(SP)"))
else:
    fifi.write("%4s & %8s & %8s & %8s \\\\ \hline \n" % ("$s$",     "$\sigma \, ("+unitlon+"$/"+opt.unit+")",      "period ("+opt.unit+")","log(SP)"))
#fifi.write("---------------------------------------\n")
# -- initialize while loop
search = np.empty_like(spec) ; search[:,:] = spec[:,:]
zelab = search > 0 # (all elements)
itit = 1 
spowermax = -9999.
# -- while loop
while itit <= opt.ndom:
  # -- find dominant mode
  ij = maximum_position(search,labels=zelab)
  dominant_wn = specx[ij[0]]
  dominant_fq = spect[ij[1]]
  spower = search[ij]
  # -- break if a given orders of magnitude difference in power
  if spower < spowermax/100.:
    break
  if (1./dominant_fq) < lowerperiod: reliable = "x"
  else: reliable = "o"
  # -- print result
  if reliable == "o":

    # -- correct by zonalwind if needed
    if opt.zonalwind is not None:
        dominant_fq_int, dominant_wn_int = intrinsic(360.*dominant_fq,dominant_wn)
        dominant_fq_int = dominant_fq_int/360.
        if (1./dominant_fq_int) > higherperiod: # pour les stationnaires
            reliable = "x"
        else:
            fifi.write("%+4.0f & %8.1f & %8.1f & %8.1f \\\\ \n" % (dominant_wn_int,fac*dominant_fq_int,1./dominant_fq_int,np.log10(spower)))
    else:
        fifi.write("%+4.0f & %8.1f & %8.1f & %8.1f \\\\ \n" % (dominant_wn    ,fac*dominant_fq    ,1./dominant_fq    ,np.log10(spower)))

    search[ij[0],:] = -9999. # remove wavenumber found (otherwise loop could find other maxima for this wn)
    if spower > spowermax: spowermax = spower
  # -- iterate
  zelab = search < spower # remove maximum found
  itit += 1

# -- close output
fifi.write("---------------------------------------\n")
fifi.close()
# -- print results
print(open(txtfile, "r").read())

## COMPUTE FREQUENCY/PERIOD AXIS
if not opt.period:
  # frequency: longitude degree per UNIT
  spect = fac*spect
  if opt.ymax is None: opt.ymax = fac/lowerperiod
  elif opt.ymin is None: opt.ymin = fac/higherperiod
else:
  # period: UNIT
  spect = 1./(spect)
  if opt.ymin is None: opt.ymin = lowerperiod
  elif opt.ymax is None: opt.ymax = higherperiod

## PLOT
if not opt.noplot:
  p = ppplot.plot2d()
  p.transopt(opt) # transfer plot options
  #p.f = spec/np.max(spec)
  p.f = spec/1e6
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

  #####################################
  #nutab = [-5,-4,-3,-2,-1,0,+1,+2,+3,+4,+5]
  #nutab = [-2,-1,0,+1,+2]
  nutab = [-1,0,+1]
  #nutab = [-1,+1]
  #####################################

  if opt.nutab == 0: nutab = [0]
  elif opt.nutab == 1: nutab = [-1,+1]
  elif opt.nutab == 2: nutab = [-2,-1,0,+1,+2]

  hache = opt.height

  hachetab = [hache]
  hachetab = [1000.,2000.,5000.,10000.,20000.,50000.]
  hachetab = [5000.,10000.,20000.,50000.]

  for hache in hachetab:
    term1 = 1./(4.*(mypl.H()**2))
    #T0=pp(file=opt.file,var="temp",y=opt.y,z=opt.z,x=0,t="0,10000").getf() ##--environ 120K, OK.
    term2 = mypl.N2() / (mypl.g*hache)
    #0.3e-5 SL nat2011 // 1.0e-5 LL grl2008 // avec dTdz=-0.7e-3 # proche LL
    m = np.sqrt(term2-term1)
    lz = 2.*np.pi / m
    c = np.sqrt(mypl.g*hache)
    print "EQUIVALENT HEIGHT", hache
    print "vertical wavelength [km] ", lz / 1000.
    print "phase speed KW [m/s] ", c
    print "equat Rossby rad [deg] ", np.sqrt(c/mypl.beta()) / 1e6


  # ensure number of points 
  # -- is enough for smooth lines
  # -- is not too high for efficiency
  n = 500
  specx = np.linspace(limxmin,limxmax,2*n)
  spect = np.linspace(spect.min(),spect.max(),n)
  spect = spect / fac # convert back from deglon/unit to cycle/unit
  s,sigma = np.meshgrid(specx,spect)

  # a few general settings for plots
  p.f = sigma*np.nan # trick to make it transparent
  p.clev = [0.]
  p.clab = False
  p.x = specx

  ## COMPUTE FREQUENCY/PERIOD AXIS
  if not opt.period:
    # frequency: longitude degree per UNIT
    p.y = fac*spect
  else:
    # period: UNIT
    p.y = 1./(spect)



  ## COMPUTE dispersion relationship for all modes
  p.colorbar = "viridis"
  for nnn in nutab:
#   if nnn == 0: p.ccol = "cyan"
#   elif nnn > 0: p.ccol = "magenta"
#   else: p.ccol = "red"
   for hhh in hachetab:
     if hhh == hachetab[0]: p.ccol = "blue" 
     if hhh == hachetab[1]: p.ccol = "purple"
     if hhh == hachetab[2]: p.ccol = "magenta"
     if hhh == hachetab[3]: p.ccol = "red"
     p.c = mypl.dispeqw(s,sigma,nu=nnn,h=hhh)
     p.make()


### SHOW or SAVE PLOT
if not opt.noplot:
  if opt.output is None: ppplot.show()
  else: ppplot.save(mode=opt.savef,filename=opt.output)

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

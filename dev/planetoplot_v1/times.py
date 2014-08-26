def ls2sol(lstabin):

#  Returns solar longitude, Ls (in deg.), from day number (in sol),
#  where sol=0=Ls=0 at the northern hemisphere spring equinox
     import numpy as np
     import math as m

     year_day=668.6
     peri_day=485.35
     timeperi=1.90258341759902
     e_elips=0.0934
     pi=3.14159265358979
     degrad=57.2957795130823

     if type(lstabin).__name__ in ['int','float']:
        lstab=[lstabin]
        lsout=np.zeros([1])
     else: 
        lstab=lstabin
        lsout=np.zeros([len(lstab)])

     i=0
     for ls in lstab:
         if (np.abs(ls) < 1.0e-5):
            if (ls >= 0.0):
                ls2sol = 0.0
            else:
                ls2sol = year_day
         else:

            zteta = ls/degrad + timeperi
            zx0 = 2.0*m.atan(m.tan(0.5*zteta)/m.sqrt((1.+e_elips)/(1.-e_elips)))
            xref = zx0-e_elips*m.sin(zx0)
            zz = xref/(2.*pi)
            ls2sol = zz*year_day + peri_day
            if (ls2sol < 0.0): ls2sol = ls2sol + year_day
            if (ls2sol >= year_day): ls2sol = ls2sol - year_day
         
         lsout[i]=ls2sol
         i=i+1

     return lsout

def sol2ls(soltabin):

#  convert a given martian day number (sol)
#  into corresponding solar longitude, Ls (in degr.),
#  where sol=0=Ls=0 is the
#  northern hemisphere spring equinox.
      import numpy as np
      import math as m

      year_day=668.6
      peri_day=485.35
      e_elips=0.09340
      radtodeg=57.2957795130823
      timeperi=1.90258341759902

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
 
      return solout

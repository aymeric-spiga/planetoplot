cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real function sol2ls(sol)

c  convert a given martian day number (sol)
c  into corresponding solar longitude, Ls (in degr.),
c  where sol=0=Ls=0 is the
c  northern hemisphere spring equinox.

      implicit none

c     input 
      real    sol
c     output 
      real    ls

c     local variables
      double precision    year_day,peri_day,timeperi,e_elips
      double precision    pi,radtodeg
c number of martian days (sols) in a martian year
      parameter (year_day=668.6)
c perihelion date (in sols)
      parameter (peri_day=485.35)
c orbital eccentricity
      parameter (e_elips=0.09340)
      parameter (pi=3.14159265358979)
c  radtodeg: 180/pi
      parameter (radtodeg=57.2957795130823)
c  timeperi: 2*pi*( 1 - Ls(perihelion)/ 360 ); Ls(perihelion)=250.99
      parameter (timeperi=1.90258341759902)

      double precision    zanom,xref,zx0,zdx,zteta,zz
c  xref: mean anomaly, zx0: eccentric anomaly, zteta: true anomaly      
      integer iter

      zz=(sol-peri_day)/year_day
      zanom=2.*pi*(zz-nint(zz))
      xref=dabs(zanom)

c  The equation zx0 - e * sin (zx0) = xref, solved by Newton
      zx0=xref+e_elips*dsin(xref)
      do 110 iter=1,10
         zdx=-(zx0-e_elips*dsin(zx0)-xref)/(1.-e_elips*dcos(zx0))
         if(dabs(zdx).le.(1.d-7)) then
           goto 120
         endif
         zx0=zx0+zdx
  110 continue
  120 continue
      zx0=zx0+zdx
      if(zanom.lt.0.) zx0=-zx0

c compute true anomaly zteta, now that eccentric anomaly zx0 is known
      zteta=2.*datan(dsqrt((1.+e_elips)/(1.-e_elips))*dtan(zx0/2.))

c compute Ls
      ls=real(zteta-timeperi)
      if(ls.lt.0.) ls=ls+2.*real(pi)
      if(ls.gt.2.*pi) ls=ls-2.*real(pi)
c convert Ls in deg.
      ls=real(radtodeg)*ls

      sol2ls = ls
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real function ls2sol(ls)

c  Returns solar longitude, Ls (in deg.), from day number (in sol),
c  where sol=0=Ls=0 at the northern hemisphere spring equinox

      implicit none

c  Arguments:
      real ls

c  Local:
      double precision xref,zx0,zteta,zz
      double precision peri_day,timeperi,e_elips
      double precision pi,degrad 
      double precision,parameter :: year_day=668.6 
      parameter (peri_day=485.35) 
      parameter (timeperi=1.90258341759902)
      parameter (e_elips=0.0934)  
      parameter (pi=3.14159265358979)
      parameter (degrad=57.2957795130823)

      if (abs(ls).lt.1.0e-5) then
         if (ls.ge.0.0) then
            ls2sol = 0.0
         else
            ls2sol = real(year_day)
         end if
         return
      end if

      zteta = ls/degrad + timeperi
      zx0 = 2.0*datan(dtan(0.5*zteta)/dsqrt((1.+e_elips)/(1.-e_elips)))
      xref = zx0-e_elips*dsin(zx0)
      zz = xref/(2.*pi)
      ls2sol = real(zz*year_day + peri_day)
      if (ls2sol.lt.0.0) ls2sol = ls2sol + real(year_day)
      if (ls2sol.ge.year_day) ls2sol = ls2sol - real(year_day)

      return
      end

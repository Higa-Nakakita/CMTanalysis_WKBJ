c=======================================================================
      real*8 function geocen(arg)
c Returns geocentric colat given geographic colat as arg 
c  (n.b. fac=(1-f)**2,   argument in radians.
c     written by G. Masters (UCSD)
      real*8 arg
      real*8 pip2,fac
      data pip2,fac/1.570796326794895d0,0.993305621334896d0/
      geocen=pip2-datan(fac*dcos(arg)/dmax1(1.d-200,dsin(arg)))
      return
      end
c=======================================================================
      real*8 function geograph(arg)
c Returns geographic colat given geocentric colat as arg 
c  (n.b. fac=(1-f)**2,   argument in radians.
c     modified by K. Yoshizawa (RSES,ANU)
      real*8 arg
      real*8 pip2,fac
      data pip2,fac/1.570796326794895d0,0.993305621334896d0/
      geograph=pip2-datan(dcos(arg)/fac/dmax1(1.d-200,dsin(arg)))
      return
      end
c=======================================================================
      real*8 function geocen_deg(arg)
c Returns geocentric latitude given geographic latitude as arg (deg) 
c     K. Yoshizawa (RSES,ANU)
      real*8 arg,geocen,glat,clat
      include 'wkbj_common.inc'
      glat = (90d0-arg)*d2r
      clat = geocen(glat)
      geocen_deg = 90d0 - clat*r2d
      !write(*,*) arg,glat,clat,geocen_deg
      return
      end
c=======================================================================
      real*8 function geograph_deg(arg)
c Returns geographic latitude given geocentric latitude as arg (deg) 
c     K. Yoshizawa (RSES,ANU)
      real*8 arg,geograph,glat,clat
      include 'wkbj_common.inc'
      clat = (90d0-arg)*d2r
      glat = geograph(clat)
      geograph_deg = 90d0 - glat*r2d
      return
      end

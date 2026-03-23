
      integer, parameter :: maxdata=5000, numcomp=3, maxwindow=4
     &                     ,maxstas=150, nmx=20, lmx=560, jmx=510
     &                     ,maxdivs=500

      real(8), parameter :: prdmin=10d0, prdmax=500d0

      ! For normalization
      real(8) G,RHO,RA,EM,unit_seis
      data G,RHO,RA,EM/6.6723d-11,5.515d3,6.371d6,5.97d24/
      common /unitseis/unit_seis

      real(4) all_params(12,70000) !nsample_max*(nitr_max+1) 
      common /writeresult/all_params

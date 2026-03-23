c
c     band_pass : applies recursive (IIR) band-pass filters to
c                     a seismogram
c
c     input: x(time) ... an original seismogram in time domain
c
c     output: ywave(time,J) ... filtered seismograms and their envelopes
c             *NOTE J:   1 - J/2 : filtered seismograms
c                      J/2+1 - J : envelopes 
c     
c     K. Yoshizawa (RSES, ANU)
c     Jan. 2001.
c     rearranged by Higa 2024

      subroutine bandpass_fil(x,yy,iw,k)
      implicit real*8(a-h,o-x)
      implicit real*4(y)
      implicit complex*16(z)
      include 'swfi_param.inc'
      include 'swfi_common.inc'
      include 'wkbj_common.inc'
      real(4) x(maxdata), yy(maxdata)
      integer irek, npoles, npts, nft, iw, k
      real(4) ytd, yf1, yf2

      irek = 1
      npoles = 4
      ytd = real(tint)
      npts = n_data(k,iw)
      nft = n_fft(k,iw)

      yf1 = bpfrq_min(iw)
      yf2 = bpfrq_max(iw)

      call bwfilt(x,yy,ytd,npts,irek,npoles,yf1,yf2)
      return
      end

c-----------------------------------------------------------------bwfilt
c
        subroutine bwfilt(x,y,dt,n,irek,norder,f1,f2)
c
c   recursive filtering of data with butterworth filter
c   x:   input array
c   y:   output array
c   dt:  time increment
c   n:   number of data points
c   irek=0  : forward filtering only
c        1  : forward and backward filtering
c   norder  : order of butterworth filter
c         =0: only filtering, no determination of coefficients
c      .lt.0: no starplots of transfer function and impulse response
c   f1      : low cutoff frequency (hz)
c         =0: low pass filter
c   f2      : high cutoff frequency (hz)
c .gt.0.5/dt: high pass filter
c
        dimension x(1),y(1)
        dimension a(10),b1(10),b2(10)
CC        iunit=3
        if(norder.eq.0) goto 100
        npoles=iabs(norder)
c                            **determination of filter coefficients
        call bpcoeff(f1,f2,npoles,dt,a,b1,b2)
c                            ** filtering
100     if(n.ne.0)call rekurs(x,y,n,a,b1,b2,npoles,irek)

        return
        end
c-------------------------------------------------------------------rekurs
c
c   programs for recursive butterworth filters
c
        subroutine rekurs(x,y,ndat,a,b1,b2,npoles,iflag)
c
c  performs recursive filtering of data in array x of length ndat
c  filtered output in y
c  a, b1, b2 are the filter coefficients 
c                    (previously determined in bwcoef)
c  npoles is the number of poles
c  iflag =0: forward filtering only
c  iflag.ne.0: forward and backward filtering
c
        dimension x(1),y(1),a(1),b1(1),b2(1),
     &    z(10),z1(10),z2(10)
c
c  forward
c
        x1=0.
        x2=0.
        do 5 i=1,npoles
        z1(i)=0.
5       z2(i)=0.
c
        do 10 n=1,ndat
        z(1)=a(1)*(x(n)-x2)-b1(1)*z1(1)-b2(1)*z2(1)
        do 15 i=2,npoles
15      z(i)=a(i)*(z(i-1)-z2(i-1))-b1(i)*z1(i)-b2(i)*z2(i)
        x2=x1
        x1=x(n)
        do 17 i=1,npoles
        z2(i)=z1(i)
17      z1(i)=z(i)
        y(n)=z(npoles)
10      continue
c
        if(iflag.eq.0) return
c
c  backward
        x1=0.
        x2=0.
        do 24 i=1,npoles
        z1(i)=0.
24      z2(i)=0.
c
        do 20 n=ndat,1,-1
        z(1)=a(1)*(y(n)-x2)-b1(1)*z1(1)-b2(1)*z2(1)
        do 25 i=2,npoles
25      z(i)=a(i)*(z(i-1)-z2(i-1))-b1(i)*z1(i)-b2(i)*z2(i)
        x2=x1
        x1=y(n)
        do 27 i=1,npoles
        z2(i)=z1(i)
27      z1(i)=z(i)
        y(n)=z(npoles)
20      continue

        return
        end
c---------------------------------------------------------------------bpcoeff
c
        subroutine bpcoeff(f1,f2,npoles,dt,a,b1,b2)
c
c           determines filter coefficients 
c           for the recursive bandpass filter
c
        dimension a(10),b1(10),b2(10)
        complex s(12),t1,t2,p
        data pi/3.14159265/
c
        if(npoles.gt.10) stop  ' npoles greater 10: stop'
        d2=2./dt
        w1=d2*tan(2.*pi*f1/d2)
        w2=d2*tan(2.*pi*f2/d2)
        w0=0.5*(w2-w1)
c
        i=1
        npol2=npoles/2+1
        do 10 n=1,npoles
        p=cexp(cmplx(0.,float(2*n-1+npoles)*pi/float(2*npoles)))
        t1=p*w0
        t2=csqrt(t1*t1-w1*w2)
        s(i)=t1+t2
        s(i+1)=t1-t2
        i=i+2
10      continue
        do 20 n=1,npoles
        ssum=2.*real(s(n))
        sprod=real(s(n)*conjg(s(n)))
        fact1=d2*d2-d2*ssum+sprod
        fact2=2.*(sprod-d2*d2)
        fact3=d2*d2+d2*ssum+sprod
        a(n)=2.*d2*w0/fact1
        b1(n)=fact2/fact1
        b2(n)=fact3/fact1
20      continue

        return
        end


c-----------------------------------------------------------
c
c     gaus_filter : calculates gaussian filter
c
c     K. Yoshizawa (RSES, ANU)
c     Jan. 2001.
c
      real*8 function gaus_filter(frq,frq0,alpha,beta,npow)
      implicit real*8(a-h,o-x)

      range1 = (1d0 - beta) * frq0
      range2 = (1d0 + beta) * frq0

      if(frq.lt.range1.or.frq.gt.range2) then
         gaus_filter = 0d0
      else
         gaus_tmp = - alpha * dabs((frq-frq0)/frq0) ** npow
         gaus_filter = dexp(gaus_tmp)
      endif
c
      return
      end

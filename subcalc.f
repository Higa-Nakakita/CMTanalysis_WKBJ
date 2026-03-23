      ! === compute azimuth,back azimuth, epicentral distance
      subroutine azdel(theSc,phiSc,theRc,phiRc,az,baz,del)
      implicit none
      include 'wkbj_common.inc'
      real(8) theSc,phiSc,theRc,phiRc,az,baz,del
      real(8) Slat,Rlat,Dphi,Slat_rad,Rlat_rad,Dphi_rad,del_rad
      real(8) BBtmp,CCtmp,BB,CC
 
      ! latitude -> colatitude
      Slat = 90d0 - theSc
      Rlat = 90d0 - theRc
      Dphi = phiRc - phiSc

      ! evaluate Dphi
      if (Dphi > 180d0) then
              Dphi = Dphi - 360d0
      elseif (Dphi < -180) then
              Dphi = Dphi + 360d0
      endif

      ! deg -> rad
      Slat_rad = Slat * d2r
      Rlat_rad = Rlat * d2r
      Dphi_rad = Dphi * d2r

      !compute distance from src to rcv 
      del_rad = dacos(dcos(Slat_rad)*dcos(Rlat_rad)  
     &  +  dsin(Slat_rad)*dsin(Rlat_rad)*dcos(Dphi_rad))
      del = del_rad * r2d

      !compute Azimuth and Back Azimuth
      CCtmp = (dcos(Rlat_rad)-dcos(del_rad)*dcos(Slat_rad)) 
     &              /(dsin(del_rad)*dsin(Slat_rad))
      BBtmp = (dcos(Slat_rad)-dcos(del_rad)*dcos(Rlat_rad)) 
     &              /(dsin(del_rad)*dsin(Rlat_rad))
      CCtmp = dmax1(-1.0d0,CCtmp)
      CCtmp = dmin1(1.0d0,CCtmp)
      BBtmp = dmax1(-1.0d0,BBtmp)
      BBtmp = dmin1(1.0d0,BBtmp)
      CC = dacos(CCtmp)
      BB = dacos(BBtmp)
      az  = CC * r2d
      baz = BB * r2d

      ! evaluate Azimuth and Back Azimuth for larger values(more thean 180 deg)
      if (Dphi<0) then
              az = 360d0 - az
      endif
      if (Dphi>0) then
              baz = 360d0 - baz
      endif
      return 
      end


      real*8 function spline(n,x,y,xnn)
      implicit real*8(a-h,o-z)
      dimension x(n),y(n),wh(n),ws(n),bb(n-1),cc(n-1),dd(n-1)
      xn = xnn
c     special case for linear interpolation
      if(n.eq.2) then
         dx = x(2) - x(1)
         aaa = (y(2) - y(1)) / dx
         bbb = (y(1)*x(2) - y(2)*x(1)) / dx
         spline = aaa * xn + bbb
         return
      endif
      do 10 i=1,n-1
         wh(i)=x(i+1)-x(i)
 10      dd(i)=(y(i+1)-y(i))/wh(i)
      cc(2)=2d0*(wh(1)+wh(2))
      bb(2)=dd(2)-dd(1)
      do 20 i=3,n-1
         cc(i)=2d0*(wh(i-1)+wh(i))-wh(i-1)*wh(i-1)/cc(i-1)
 20      bb(i)=dd(i)-dd(i-1)-bb(i-1)*wh(i-1)/cc(i-1)
      ws(1)=0.0d0
      ws(n)=0.0d0
      ws(n-1)=bb(n-1)/cc(n-1)
      do i=n-2,2,-1
         ws(i)=(bb(i)-wh(i)*ws(i+1))/cc(i)
      end do
      do i=1,n-1
         bb(i)=dd(i)-wh(i)*(ws(i+1)+2d0*ws(i))
         cc(i)=3d0*ws(i)
         dd(i)=(ws(i+1)-ws(i))/wh(i)
      end do
      i=1
 30   if(i.ge.(n-1).or.(xn.le.x(i+1))) go to 40
      i=i+1
      go to 30
 40   xn=xn-x(i)
      spline=y(i)+(bb(i)+(cc(i)+dd(i)*xn)*xn)*xn
      return
      end


c------------------------------------------------------------------------------
c     
c     fast: fast Fourier transform
c       ind = -1: FFT
c              1: IFFT (x must be divided by N in a main program)
c
c     from Ohsaki (1994) (modified for double precision version)
c
      subroutine fast(n,x,nd,ind)
      complex*16  x(nd),temp
      integer ind,kmax,istep ,n,nd
      complex*16 zexp_i
      real*8 zzz,theta
c      zexp_i(zzz) = dcmplx(dcos(zzz),0d0) + dcmplx(0d0,dsin(zzz))
      zexp_i(zzz) = dcmplx(dcos(zzz),dsin(zzz))

      !write(*,*) 'n,nd',n,nd,ind
      !write(*,*) x
      pi = 4d0*datan(1d0)
      j = 1
      do 140 i=1,n
         if(i.ge.j) go to 110
         temp=x(j)
         x(j)=x(i)
         x(i)=temp
 110     m=n/2
 120     if(j.le.m) go to 130
         j=j-m
         m=m/2
         if(m.ge.2) go to 120
 130     j=j+m
 140  continue
      kmax=1
 150  if(kmax.ge.n) return
      istep=kmax*2
      do 170 k=1,kmax
c         theta=dcmplx(0d0,pi*dble(ind)*dble(k-1)/dble(kmax))
         theta = pi * dble(ind) * dble(k-1) / dble(kmax)
         do 160 i=k,n,istep
            j=i+kmax
            temp=x(j)*zexp_i(theta)
            x(j)=x(i)-temp
            x(i)=x(i)+temp
 160     continue
 170  continue
      kmax=istep
      go to 150
      end

      subroutine changed_fast(n,x,nd,ind)
      implicit none
      integer ind,kmax,istep ,n,nd,i,j,k,m
      complex(16)  x(nd),temp
      complex(16) zexp_i
      real(8) zzz,theta,pi
c      zexp_i(zzz) = dcmplx(dcos(zzz),0d0) + dcmplx(0d0,dsin(zzz))
      zexp_i(zzz) = dcmplx(dcos(zzz),dsin(zzz))

      pi = 4d0*datan(1d0)
      j = 1
      do 140 i=1,n
         if(i.ge.j) go to 110
         temp=x(j)
         x(j)=x(i)
         x(i)=temp
 110     m=n/2
 120     if(j.le.m) go to 130
         j=j-m
         m=m/2
         if(m.ge.2) go to 120
 130     j=j+m
 140  continue
      kmax=1
 150  if(kmax.ge.n) return
      istep=kmax*2
      do 170 k=1,kmax
c         theta=dcmplx(0d0,pi*dble(ind)*dble(k-1)/dble(kmax))
         theta = pi * dble(ind) * dble(k-1) / dble(kmax)
         do 160 i=k,n,istep
            j=i+kmax
            temp=x(j)*zexp_i(theta)
            x(j)=x(i)-temp
            x(i)=x(i)+temp
 160     continue
 170  continue
      kmax=istep
      go to 150
      end


      ! === word counter ====
      integer function lofw(word)
      character(*) word
      integer i, k, lw
      lw = len(word)
      k = 0
      do i =1, lw
         if(word(i:i).eq.' ') go to 99
         k = k + 1
      enddo
  99  lofw = k
      return
      end


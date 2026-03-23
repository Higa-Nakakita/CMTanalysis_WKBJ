
      subroutine srcterm_renew(rmodel,newdep,newazim,newdel,STmode,k)
      ! Renew the source term (See Dahlen&Tromp (1998) equation 16.170)
      ! Firstly, eigenfunction value at new source depth will be
      ! obatined using Akima-interporaion.
      ! Secondly, radiation pattern calculation will be performed 
      ! for the new moment tensor. half-duration correction will be
      ! applied after that.
      ! Then some adjustment related to unit will be applied.
      ! Return tha complex spectrum. 
      implicit none
      include 'swfi_param.inc'
      include 'swfi_common.inc'
      include 'wkbj_common.inc'
      include 'eigen_common.inc'
      include 'new_eigen_common.inc'

      real(8) rmodel(12),newdep,newazim,newdel
      integer STmode,k

      real(8) Mrr,Mtt,Mpp,Mrt,Mrp,Mtp,scalar,halfd

      real(8) rtmp, az_s, sin_D, geosp_tmp, sin_z, cos_z, sin_2z, cos_2z
      complex(16) zexp_posi, zexp_nega

      integer n, l, narray(0:nmx,2), ntmp, model_src
      real(8) dl(lmx), omgsrc ,omg_src(lmx)
      real(8) eif_srcU(lmx),eif_srcdU(lmx),eif_srcV(lmx),eif_srcdV(lmx)
      real(8) eif_srcW(lmx), eif_srcdW(lmx)

      integer j, jjtmp
      real(8) frqtmp, omgtmp
      real(8) j_skip1_src(0:nmx,2),j_skip2_src(0:nmx,2)
      real(8) spline, dl_src, dks
      real(8) Us, dUs, Vs, dVs, Ws, dWs
      real(8) Atmp1,Atmp2_1,Atmp2_2,Atmp3_1,Atmp3_2
      complex(16) zAexp_s,zAtmp,zvel,zacc

      ! eigen func of source
      real(8) f, df, ddf
      real(4) Us_z(1:numlay),dUs_z(1:numlay)
      real(4) Vs_z(1:numlay),dVs_z(1:numlay)
      real(4) Ws_z(1:numlay),dWs_z(1:numlay)
      real(8) Usrc(0:nmx,lmx),dUsrc(0:nmx,lmx)
      real(8) Vsrc(0:nmx,lmx),dVsrc(0:nmx,lmx)
      real(8) Wsrc(0:nmx,lmx),dWsrc(0:nmx,lmx)

      ! halfduration correction
      real(8) con, arg, sinc2, freal, fimag
      complex(16) dexp_MT 

      ! For calc radiation pattern
      integer iaz, maxradas
      real(4) radia_sum(2,360)

      model_src = 2

      Mrr   = rmodel(5)
      Mtt   = rmodel(6)
      Mpp   = rmodel(7)
      Mrt   = rmodel(8)
      Mrp   = rmodel(9)
      Mtp   = rmodel(10)
      scalar= rmodel(11)
      halfd = rmodel(12) 

      !write(*,*) "Mrr,Mtt,Mpp,Mrt,Mrp,Mtp,scalar,halfd"
      !write(*,*)  Mrr,Mtt,Mpp,Mrt,Mrp,Mtp,scalar,halfd

      rtmp = 1d0 - newdep / 6371d0

      az_s  = pi - newazim
      sin_D = dsin(newdel)
      geosp_tmp = 8d0*pi*dabs(sin_D)
      sin_z  = dsin(az_s)
      cos_z  = dcos(az_s)
      sin_2z = dsin(2d0*az_s)
      cos_2z = dcos(2d0*az_s)
      zexp_posi = dcmplx(dcos(pi/4d0),dsin(pi/4d0))
      zexp_nega = dcmplx(dcos(pi/4d0),-dsin(pi/4d0))


      !write(*,*) 'newdel,newazim, az_s,sin_z,cos_z,sin_2z,cos_2z'
      !write(*,*) newdel,newazim,az_s,sin_z,cos_z,sin_2z,cos_2z

      !open(1,file='ck_sinc.txt')
      do n = nn1, nn2
        narray(n,STmode) = 0

        do l = lbeg_path(n,STmode,k), lend_path(n,STmode,k), mode_skip
          narray(n,STmode) = narray(n,STmode) + 1
          ntmp = narray(n,STmode)
          dl(ntmp) = dble(l)
          omg_src(ntmp) = dble(eig_para(model_src,STmode,n,l,1,k))
          omgsrc = omg_src(ntmp)

          ! store eigen function of source corresponds to newdep
          Us_z(1:numlay) = src_eifS_z(n,l,1:numlay,1)
          dUs_z(1:numlay)= src_eifS_z(n,l,1:numlay,2)
          Vs_z(1:numlay) = src_eifS_z(n,l,1:numlay,3)
          dVs_z(1:numlay)= src_eifS_z(n,l,1:numlay,4)
          Ws_z(1:numlay) = src_eifT_z(n,l,1:numlay,1)
          dWs_z(1:numlay)= src_eifT_z(n,l,1:numlay,2)

          ! Akima interpolation for depth of eigen func
          call makima1p_modi(f,df,ddf,
     &        dble(newdep),numlay,dble(Us_z),dble(src_depth(1:numlay)))
          Usrc(n,l) = f

          call makima1p_modi(f,df,ddf,
     &        dble(newdep),numlay,dble(dUs_z),dble(src_depth(1:numlay)))
          dUsrc(n,l) = f
          
          call makima1p_modi(f,df,ddf,
     &        dble(newdep),numlay,dble(Vs_z),dble(src_depth(1:numlay)))
          Vsrc(n,l) = f

          call makima1p_modi(f,df,ddf,
     &        dble(newdep),numlay,dble(dVs_z),dble(src_depth(1:numlay)))
          dVsrc(n,l) = f

          call makima1p_modi(f,df,ddf,
     &        dble(newdep),numlay,dble(Ws_z),dble(src_depth(1:numlay)))
          Wsrc(n,l) = f

          call makima1p_modi(f,df,ddf,
     &        dble(newdep),numlay,dble(dWs_z),dble(src_depth(1:numlay)))
          dWsrc(n,l) = f
          
          ! Interpolation will be performed for the Eigen function
          ! above the free surface when dep < 0
          ! Its value does not make sense
          ! Thus Let us avoid it.
          if (newdep < 0d0) then  
             Usrc(n,l)  = 0d0
             dUsrc(n,l) = 0d0
             Vsrc(n,l)  = 0d0
             dVsrc(n,l) = 0d0
             Wsrc(n,l)  = 0d0
             dWsrc(n,l) = 0d0
          endif

          if (STmode==1) then ! Spheroidal mode
             eif_srcU(ntmp)  = dble(Usrc(n,l))  * omgsrc
             eif_srcdU(ntmp) = dble(dUsrc(n,l)) * omgsrc
             eif_srcV(ntmp)  = dble(Vsrc(n,l))  * omgsrc
             eif_srcdV(ntmp) = dble(dVsrc(n,l)) * omgsrc

          elseif(STmode==2) then ! Toroidal mode
             eif_srcW(ntmp)  = dble(Wsrc(n,l))  * omgsrc
             eif_srcdW(ntmp) = dble(dWsrc(n,l)) * omgsrc

          endif

          !if (k==1) then
          !   open(3,file='ck_eifsrc_nl.txt',position='append')
          !   write(3,*) n,l,real(Usrc(n,l)),real(dUsrc(n,l)),
          !&               real(Vsrc(n,l)),real(dVsrc(n,l)),
          !&               real(Wsrc(n,l)),real(dWsrc(n,l))
          !   close(3)
          !endif

         ! if (STmode==1) then ! Spheroidal mode
         !    eif_srcU(ntmp)  = dble(eif_srcS(n,l,1,k)) * omgsrc
         !    eif_srcdU(ntmp) = dble(eif_srcS(n,l,2,k)) * omgsrc
         !    eif_srcV(ntmp)  = dble(eif_srcS(n,l,3,k)) * omgsrc
         !    eif_srcdV(ntmp) = dble(eif_srcS(n,l,4,k)) * omgsrc

         ! elseif(STmode==2) then ! Toroidal mode
         !    eif_srcW(ntmp)  = dble(eif_srcT(n,l,1,k)) * omgsrc
         !    eif_srcdW(ntmp) = dble(eif_srcT(n,l,2,k)) * omgsrc
         ! endif

        enddo ! enddo of l

        jjtmp = 0
        do j = 1, nyquist(k)
          frqtmp = frqprt(j)
          omgtmp = omgprt(j)
          !write(*,*) j, frqtmp, omgtmp

          if (frqtmp < freq_min(n,STmode,k)) then
             j_skip1_src(n,STmode) = j
             zAexp_src(n,j,STmode,k) = dcmplx(0d0,0d0)
          elseif (frqtmp > freq_max(n,STmode,k)) then
             zAexp_src(n,j,STmode,k) = dcmplx(0d0,0d0)
          else

             ! ==== Correction for half-duration
             ! See textbook zishingaku (2015) equation 13.28 (p182)
             con = 0.5d0 * (halfd*2d0) * omgtmp 
             arg = 0.5d0 * con
             sinc2 = (dsin(arg)/arg) **2
             freal = sinc2 * dcos(con)
             fimag = -sinc2 * dsin(con)
             dexp_MT = dcmplx(freal,fimag)
             !write(1,*) n,omgtmp, 
             !&       real(real(dexp_MT)), real(aimag(dexp_MT))
             ! ==== halfd calc done

             j_skip2_src(n,STmode) = j
             jjtmp = jjtmp + 1
             dl_src = spline(ntmp,omg_src,dl,omgtmp)
             dks = dsqrt(dl_src*(dl_src + 1d0))
             !write(*,*) dl_src,jjtmp,dks

             if (STmode==1) then ! Spheroidal mode
                Us  = spline(ntmp,omg_src,eif_srcU,omgtmp)
                dUs = spline(ntmp,omg_src,eif_srcdU,omgtmp)
                Vs  = spline(ntmp,omg_src,eif_srcV,omgtmp)
                dVs = spline(ntmp,omg_src,eif_srcdV,omgtmp)
                Atmp1   = Mrr*dUs + (Mtt+Mpp)*(Us-0.5d0*dks*Vs)/rtmp
                Atmp2_1 = dVs - Vs/rtmp + dks * Us / rtmp
                Atmp2_2 = (Mrp * sin_z + Mrt * cos_z)
                Atmp3_1 = dks * Vs/rtmp
                Atmp3_2 = (Mtp * sin_2z + 0.5d0*(Mtt-Mpp)*cos_2z)
                zAtmp = (Atmp1-Atmp3_1*Atmp3_2) * zexp_posi
     &                 - Atmp2_1 * Atmp2_2 * zexp_nega
               !write(*,*) Atmp1,Atmp2_1,Atmp2_2,Atmp3_1,Atmp3_2,zAtmp

             elseif (STmode==2) then ! Toroidal mode
                Ws  = spline(ntmp,omg_src,eif_srcW,omgtmp)
                dWs = spline(ntmp,omg_src,eif_srcdW,omgtmp)
                Atmp1 = (dWs-Ws/rtmp)*(Mrt*sin_z - Mrp*cos_z)
                Atmp2_1 = dks * Ws / rtmp
                Atmp2_2 = 0.5d0*(Mtt-Mpp)*sin_2z - Mtp*cos_2z
                zAtmp = -Atmp1*zexp_nega -Atmp2_1*Atmp2_2*zexp_posi
             endif

             ! make halfd coorection
             zAtmp = zAtmp * dexp_MT

             if (iobstype==1) then
                 zAexp_s = - zAtmp / omgtmp
             elseif (iobstype==2) then
                 zvel = dcmplx(0d0,omgtmp)
                 zAexp_s = - zvel * zAtmp / omgtmp
             elseif (iobstype==3) then
                 zacc = dcmplx(0d0,omgtmp) * dcmplx(0d0,omgtmp)
                 zAexp_s = - zacc * zAtmp / omgtmp
             endif
             zAexp_src(n,j,STmode,k) = zAexp_s

             !if (k == 1) then ! Only execue when first rcv
             !call calc_radpattern(STmode,n,j,dks,rtmp,
             !&          rmodel,Us,dUs,Vs,dVs,Ws,dWs,radia_sum)
             !endif

          endif
        enddo 
      enddo ! enddo of n
      !close(1)

      !if (k==1) then
      !   maxradas = 360
      !   open(3,file='ck_radpattern.txt',position='append')
      !   do iaz = 1, maxradas
      !      write(3,*) STmode, iaz, radia_sum(STmode,iaz)
      !      !write(*,*) STmode, iaz, radia_sum(STmode,iaz)
      !   enddo
      !   close(3)
      !endif


      !open(3,file='ck_zAexpsrc.txt',position='append')
      !do n = nn1, nn2
      !   do j = int(j_skip1_src(n,STmode)), int(j_skip2_src(n,STmode))
      !      write(3,*) n,j,real(frqprt(j)),STmode,k,
      !&    real(real(zAexp_src(n,j,STmode,k))),
      !&    real(aimag(zAexp_src(n,j,STmode,k)))
      !   enddo
      !enddo
      !close(3)

      endsubroutine srcterm_renew

      subroutine reflect_M0(scalar,de_nrm,Mrr,Mtt,Mpp,Mrt,Mrp,Mtp)
      ! Compute the scaling factor that reflects the seismic moment 
      !in the waveform amplitude.
      !In this calculation, only the exponent of the MT is used,
      ! and the mantissa is excluded. (251117)
      implicit none
      include "swfi_param.inc"
      real(8) scalar, de_nrm
      real(8) M0, dtmp, Mrr,Mtt,Mpp,Mrt,Mrp,Mtp

      M0 = 1/sqrt(2d0) * 10d0 ** scalar 
      ! You must NOT include mantissa of tensor in this calculation
      ! because mantissa is already reflected in source term.
      ! (as a weight of eigen function)
      !M0 =M0 * sqrt(Mrr**2d0+Mtt**2d0+Mpp**2d0+
      !&  2d0 * (Mrt**2d0+ Mrp**2d0+ Mtp**2d0)) 
      dtmp = G * RHO * EM *RA * (3d0/4d0)
      de_nrm = M0 * 1d-7 * unit_seis / dtmp
      
      endsubroutine reflect_M0


      subroutine calc_radpattern(STmode,n,j,dks,rtmp,
     &            rmodel,Us,dUs,Vs,dVs,Ws,dWs,radia_sum)
      implicit none
      include 'swfi_param.inc'
      include 'swfi_common.inc'
      include 'wkbj_common.inc'
      include 'eigen_common.inc'
      include 'new_eigen_common.inc'

      integer STmode, n, j
      integer iaz, maxradaz
      real(8) omgtmp
      real(8) rmodel(12),Mrr,Mtt,Mpp,Mrt,Mrp,Mtp
      real(8) dks, rtmp, Us,dUs,Vs,dVs,Ws,dWs
      real(8) aztmp, az_s, sin_z, cos_z, sin_2z, cos_2z 
      real(8) Atmp1,Atmp2_1,Atmp2_2,Atmp3_1,Atmp3_2
      complex(16) zAexp_s,zAtmp,zvel,zacc
      complex(16) zexp_posi, zexp_nega
      real(8) radia
      real(4) radia_sum(2,360)
      !real(8) radia_array(2,0:nmx,maxdata,360)

      zexp_posi = dcmplx(dcos(pi/4d0),dsin(pi/4d0))
      zexp_nega = dcmplx(dcos(pi/4d0),-dsin(pi/4d0))
      Mrr   = rmodel(5)
      Mtt   = rmodel(6)
      Mpp   = rmodel(7)
      Mrt   = rmodel(8)
      Mrp   = rmodel(9)
      Mtp   = rmodel(10)

      omgtmp = omgprt(j)
      maxradaz = 360
      do iaz = 1, maxradaz
         aztmp = dble(iaz-1) * d2r
         az_s = pi - aztmp
         sin_z = dsin(az_s)
         cos_z = dcos(az_s)
         sin_2z= dsin(2d0*az_s)
         cos_2z= dcos(2d0*az_s)

         if (STmode==1) then ! Spheroidal mode
             Atmp1   = Mrr*dUs + (Mtt+Mpp)*(Us-0.5d0*dks*Vs)/rtmp
             Atmp2_1 = dVs - Vs/rtmp + dks * Us / rtmp
             Atmp2_2 = (Mrp * sin_z + Mrt * cos_z)
             Atmp3_1 = dks * Vs/rtmp
             Atmp3_2 = (Mtp * sin_2z + 0.5d0*(Mtt-Mpp)*cos_2z)
             zAtmp = (Atmp1-Atmp3_1*Atmp3_2) * zexp_posi
     &                 - Atmp2_1 * Atmp2_2 * zexp_nega

          elseif (STmode==2) then ! Toroidal mode
             Atmp1 = (dWs-Ws/rtmp)*(Mrt*sin_z - Mrp*cos_z)
             Atmp2_1 = dks * Ws / rtmp
             Atmp2_2 = 0.5d0*(Mtt-Mpp)*sin_2z - Mtp*cos_2z
             zAtmp = -Atmp1*zexp_nega -Atmp2_1*Atmp2_2*zexp_posi
          endif

          if (iobstype==1) then
              zAexp_s = - zAtmp / omgtmp
          elseif (iobstype==2) then
              zvel = dcmplx(0d0,omgtmp)
              zAexp_s = - zvel * zAtmp / omgtmp
          elseif (iobstype==3) then
              zacc = dcmplx(0d0,omgtmp) * dcmplx(0d0,omgtmp)
              zAexp_s = - zacc * zAtmp / omgtmp
          endif

          radia = dsqrt(dble(real(zAexp_s)**2d0) 
     &            + dble(aimag(zAexp_s)**2d0))
          !radia_array(STmode,n,j,iaz) = radia

          radia_sum(STmode,iaz) = radia_sum(STmode,iaz) + real(radia)
          !write(*,*) 'r',aztmp,radia_sum(STmode,iaz),radia

      enddo

      
      endsubroutine calc_radpattern
     

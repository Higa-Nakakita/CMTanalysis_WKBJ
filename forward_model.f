
      subroutine forward_modelling(models,nd,nitr)
      ! Perform the forward calculations for each station and each component.
      implicit none
      integer k, icomp, nitr, nd
      real(4) models(nd)
      include "swfi_param.inc"
      include "wkbj_common.inc"
      include "swfi_common.inc"

      write(*,*) 
      write(*,*) 'forward_modelling nitr=',nitr
      do k = 1, rcvnum
        do icomp = 1, 3
          call wkbj_pertb(models,nd,icomp,k,nitr)
        enddo
      enddo
      !write(*,*) 'stop @forward_modelling'
      !stop

      endsubroutine forward_modelling


      subroutine wkbj_pertb(models,nd,icomp,k,nitr)
      ! Renew the theoretical seismograms,
      ! by update source terms and perturbation of phase and attenuation.
      implicit none
      include "swfi_param.inc"
      include "wkbj_common.inc"
      include "swfi_common.inc"
      include "eigen_common.inc"
      include "new_eigen_common.inc"

      integer nd, icomp, nitr, k, STmode
      real(4) models(nd)

      real(8) d_t,d_lat,d_lon,d_dep,t_shift,M0
      real(8) Mrr,Mtt,Mpp,Mrt,Mrp,Mtp,scalar,halfd
      real(8) rmodel(12)

      real(8) new_lat, new_lon, new_dep, theSc_geoc, geocen_deg
      real(8) d_del, new_azim, new_del

      integer n, j, jjtmp, jinv, nt
      real(8) dktmp, gvptmp,omgtmp,Qptmp,Atte
      complex(16) src_term,phase_pertb,atten_pertb,pertb_all
      complex(16) zwave_frq(maxdata)

      real(8) de_nrm, correct_ift
      real(4) ywave_time(maxdata)
      integer iw, nstart, nend, npts, nnr, nnl
      real(4) ytb, ytd
      real(4) syntmp(maxdata),syntmp_tap(maxdata),syntmp_fil(maxdata)

      integer nt_time
      integer jj

      real(8) calc_halfd


      if (icomp==1 .or. icomp==2) then ! Vertical or Radial
          STmode = 1  ! Speroidal mode
      elseif (icomp==3) then  ! Transverse
          Stmode = 2  ! Toroidal mode
      endif

      if (MT_CMT==0) then !MT inversion
        d_t   = 0d0
        d_lat = 0d0
        d_lon = 0d0
        d_dep = 0d0
        if (trM==0) then
          Mrr = dble(models(1))
          Mtt = dble(models(2))
          Mpp = -(Mrr+Mtt)
          Mrt = dble(models(3))
          Mrp = dble(models(4))
          Mtp = dble(models(5))
          scalar = dble(models(6))
          halfd = calc_halfd(scalar,Mrr,Mtt,Mpp,Mrt,Mrp,Mtp)
        elseif (trM==1) then
          Mrr = dble(models(1))
          Mtt = dble(models(2))
          Mpp = dble(models(3))
          Mrt = dble(models(4))
          Mrp = dble(models(5))
          Mtp = dble(models(6))
          scalar = dble(models(7))
          halfd = calc_halfd(scalar,Mrr,Mtt,Mpp,Mrt,Mrp,Mtp)
        endif
      elseif (MT_CMT==1) then !CMT inversion
        if (origT==1.or.origT==2) then !origtime: free
            d_t   = dble(models(1))
            d_lat = dble(models(2))
            d_lon = dble(models(3))
            d_dep = dble(models(4))
            Mrr   = dble(models(5))
            Mtt   = dble(models(6))
            if (trM==1) then !trM free
               Mpp   = dble(models(7))
               Mrt   = dble(models(8))
               Mrp   = dble(models(9))
               Mtp   = dble(models(10))
               scalar= dble(models(11))
               if (origT==1) then ! only dt: search
                    halfd = calc_halfd(scalar,Mrr,Mtt,Mpp,Mrt,Mrp,Mtp) 
               elseif (origT==2) then ! dt & halfd: search
                    halfd = dble(models(12))
               endif
            elseif (trM==0) then !trM=0
                Mpp   = - (Mrr+Mtt)
                Mrt   = dble(models(7))
                Mrp   = dble(models(8))
                Mtp   = dble(models(9))
                scalar= dble(models(10))
               if (origT==1) then
                    halfd = calc_halfd(scalar,Mrr,Mtt,Mpp,Mrt,Mrp,Mtp) 
               elseif (origT==2) then
                    halfd = dble(models(11))
               endif
            endif
        elseif (origT==0) then !origtime=0
            d_t   = 0d0
            d_lat = dble(models(1))
            d_lon = dble(models(2))
            d_dep = dble(models(3))
            Mrr   = dble(models(4))
            Mtt   = dble(models(5))
            if (trM==1) then !trM free
                Mpp   = dble(models(6))
                Mrt   = dble(models(7))
                Mrp   = dble(models(8))
                Mtp   = dble(models(9))
                scalar= dble(models(10))
                halfd = dble(models(11))
            elseif (trM==0) then !trM=0
                Mpp   = - (Mrr+Mtt)
                Mrt   = dble(models(6))
                Mrp   = dble(models(7))
                Mtp   = dble(models(8))
                scalar= dble(models(9))
                halfd = dble(models(10))
            endif
        else
            write(*,*) 'Wrong in origT? @forward'
            stop
        endif
      endif !MT or CMT

      d_t   = 4.0d0 !Centroid Time
      d_lat = 0d0
      d_lon = 0d0
      d_dep = 0d0
      Mrr   = 1.35d0 !1d0
      Mtt   = 0.021d0
      Mpp   = -1.370 !-1d0
      Mrt   = 0.543d0
      Mrp   = 1.300d0
      Mtp   = -0.372d0
      scalar= 25d0
      halfd = 2.8d0



      ! just check ===
      !write(*,*) calc_halfd(scalar,Mrr,Mtt,Mpp,Mrt,Mrp,Mtp) 
      !M0 = 1/sqrt(2d0) * 10d0 ** scalar
      !M0 =M0 * sqrt(Mrr**2d0+Mtt**2d0+Mpp**2d0+
      !&  2d0 * (Mrt**2d0+ Mrp**2d0+ Mtp**2d0))
      !write(*,*) 'Mo=', M0
      ! =======
     
      ! create model vector 
      rmodel = (/d_t,d_lat,d_lon,d_dep,Mrr,Mtt,Mpp,Mrt,Mrp,Mtp,
     &          scalar,halfd/)

      !write(*,*) 'd_t,d_lat,d_lon,d_dep,Mrr,Mtt,Mpp,Mrt,Mrp,Mtp',
      !&          'scalar,halfd'
      !write(*,*) rmodel
      !write(*,*) 'stop @forward_modelling'
      !stop

   
      ! Obtain perturbating disatnce and azimuth 
      ! due to the renew of the point-source location. 
      if (icomp==1) then ! Avoid overlap calc for the efficiency
         call calc_newloca(k,rmodel,d_del,new_dep,new_azim,new_del)
      endif

      ! Renew the source term.
      ! For the efficency, avoid overlap calculation.
      ! source term will be calculated based on S/T 2modes not 3comps.
      if (icomp==1 .or. icomp==3) then ! Z, T -> Spheroidal, Toroidal
         call srcterm_renew(rmodel,new_dep,new_azim,new_del,STmode,k)
      endif


      do n = nn1, nn2
         !do j = 1, j_skip1(n,STmode,k)
         do j = 1, j_skip1_adjst(n,STmode)
            zAexp_frq(n,j,STmode,icomp,k) = dcmplx(0d0,0d0)
         enddo

         !do j = j_skip1(n,STmode,k)+1, j_skip2(n,STmode,k)
         do j =  j_skip1_adjst(n,STmode)+1,j_skip2_adjst(n,STmode)
            jjtmp = j_frq(n,j,STmode)

            !! SOURCE term  
            src_term = zAexp_src(n,j,STmode,k)

            !! PATH term correction (petrurbating phase and attenuation)
            ! phase term
            dktmp = dkp(n,jjtmp,STmode)
            phase_pertb = zexp(dcmplx(0d0,-dktmp*d_del))
            ! atten term
            gvptmp = gv_array(n,jjtmp,STmode) !* dkm2d * d2r
            omgtmp = omgprt(j)
            Qptmp  = Qp_array(n,jjtmp,STmode)
            Atte   = d_del * omgtmp / (2d0*Qptmp*gvptmp)
            atten_pertb = zexp(dcmplx(-Atte,0d0))

            ! Sum all perturbation effect
            pertb_all = phase_pertb * atten_pertb * src_term

            ! Perturbation
            zAexp_frq(n,j,STmode,icomp,k) 
     &      =       zAexp_fix(n,j,STmode,icomp,k) * pertb_all
         enddo

         !do j = j_skip2(n,STmode,k)+1, nyquist(k)
         do j = j_skip2_adjst(n,STmode)+1, nyquist(k)
           zAexp_frq(n,j,STmode,icomp,k) = dcmplx(0d0,0d0)
         enddo
      enddo

      !open(4,file='ck_zAexpfrq.txt',position='append')
      !do n = nn1, nn2
      !   do j = j_skip1(n,STmode,k)+1, j_skip2(n,STmode,k)
      !      write(4,*) n,j,real(frqprt(j)),STmode,icomp,k,
      !&   real(real(zAexp_frq(n,j,STmode,icomp,k))),
      !&   real(aimag(zAexp_frq(n,j,STmode,icomp,k)))
      !    enddo
      !enddo
      !close(4)


      do j = 1, nyquist(k)
         zwave_frq(j) = dcmplx(0d0,0d0)
      enddo

      ! Create the function of freq (sum over modes)
      do j = 1, nyquist(k)
         do n = nn1, nn2
            if (frqprt(j)<freq_min(n,STmode,k) .or.
     &          frqprt(j)>freq_max(n,STmode,k)) goto 99
            zwave_frq(j) = zwave_frq(j) + zAexp_frq(n,j,STmode,icomp,k)
         enddo
99    enddo

      !open(10,file='ck_zwavefrq.txt',position='append')
      !do j = 1, nyquist
      !   write(10,*) j,real(frqprt(j)),STmode,icomp,k,
      !&   real(real(zwave_frq(j))), real(aimag(zwave_frq(j)))
      !enddo
      !close(10)

      !open(101,file='ck_synconjg.txt',position='append')
      ! set complex conlugate parts
      do j = 2, nyquist(k)
         jinv = ntmax(k) - j + 2
         zwave_frq(jinv) = conjg(zwave_frq(j))
      !   write(101,*) j,jinv,zwave_frq(jinv)
      enddo
      !close(101)

      !open(11,file='ck_zwavefrq_conjg.txt',position='append')
      !do nt = 1, ntmax
      !   write(11,*) nt,STmode,icomp,k,
      !&   real(real(zwave_frq(nt))), real(aimag(zwave_frq(nt)))
      !enddo
      !close(11)


      ! IFFT
      call changed_fast(ntmax(k),zwave_frq,maxdata,1)
      !call fast(ntmax,dble(zwave_frq),maxdata,1)
      
      ! reflect scalar moment     
      call reflect_M0(scalar,de_nrm,Mrr,Mtt,Mpp,Mrt,Mrp,Mtp)
      correct_ift = de_nrm / (tint * dble(ntmax(k)))
      !write(*,*) correct_ift
       
      ! zero padding
      do nt = 1, maxdata
         syntmp(nt) = 0d0
      enddo


      ! t_shift = centroidtime - halfd
      t_shift = d_t - halfd

      ! time correction for full-waveforms (20250205)
      do nt = 1, ntmax(k)
         !syntmp(nt) = real(dble(zwave_frq(nt)))
         !syntmp(nt) = correct_ift * syntmp(nt) 
         !nt_time = nt+int(real(d_t)+t0_buffer)
         nt_time = nt+int(real(t_shift)+t0_buffer)
         syntmp(nt_time) = real(dble(zwave_frq(nt)))
         syntmp(nt_time) = correct_ift * syntmp(nt_time) 
      enddo

      !open(12,file='ck_syntmp.txt',position='append')
      !do nt = 1, ntmax
      !   write(12,*) nt,STmode,icomp,k,
      !&   real(real(zwave_frq(nt))), !, real(aimag(zwave_frq(nt))),
      !&   real(syntmp(nt))
      !enddo
      !close(12)


      ! Convert Transverse polarity from Right(Dahlen) -> Left(Obspy)
      if (icomp==3) then
         do nt = 1, maxdata
            syntmp(nt) = syntmp(nt) * (-1d0)
         enddo
      endif
    
      ! Just for checking full synthetic waveform
      !do iw = 1, nwave
      !   call w_seis_full(syntmp,iw,1,k,icomp)
      !enddo
      !call w_seis_full(syntmp,1,1,k,icomp)

      ! partitioning of time-windows
      ! wave procesiing must be same as User_init
      do iw = 1, nwave
         ytb = t_begin(k,iw)
         nstart = nint(ytb/tint) + 1
         nend = nstart + n_data(k,iw) - 1
         npts = n_data(k,iw)
         nnl = int(npts * 0.05)
         nnr = nnl
         !write(*,*) 'nstart,nend,nnl,nnr'
         !write(*,*) nstart,nend,nnl,nnr
         call cos_tap(syntmp,syntmp_tap,nstart,nend,nnl,nnr)
         !write(*,*) 'syntmp',syntmp(1:npts)
         !write(*,*) 'syntmp_tap',syntmp_tap(1:npts)

         call bandpass_fil(syntmp_tap,syntmp_fil,iw,k)
         !write(*,*) 'syntmp_fil',syntmp_fil(1:npts)

         synthetic_data(1:npts,icomp,iw,k) = syntmp_fil(1:npts)

         call w_seis(syntmp_fil,iw,3,k,icomp)

         !if (k==7 .and. icomp==3) then stop
         !endif
      enddo

      ! store all parameters (dependent variables are included)
      jj = nitr
      all_params(1,jj) = real(d_t)
      all_params(2,jj) = real(d_lat)
      all_params(3,jj) = real(d_lon)
      all_params(4,jj) = real(d_dep)
      all_params(5,jj) = real(Mrr)
      all_params(6,jj) = real(Mtt)
      all_params(7,jj) = real(Mpp)
      all_params(8,jj) = real(Mrt)
      all_params(9,jj) = real(Mrp)
      all_params(10,jj)= real(Mtp)
      all_params(11,jj)= real(scalar)
      all_params(12,jj)= real(halfd)      

      endsubroutine wkbj_pertb


      subroutine calc_newloca(k,rmodel,d_del,new_dep,new_azim,new_del)
      ! just obtaion perurbation of delta and azimuth. 
      ! (Dahlen&Tromp (1998) eq.16.166) 
      implicit none
      include "swfi_param.inc"
      include "wkbj_common.inc"
      include "swfi_common.inc"
      integer k
      real(8) rmodel(12)
      real(8) d_del, new_dep, new_azim

      real(8) d_lat, d_lon, d_dep
      real(8) new_lat, new_lon
      real(8) geocen_deg, theSc_geoc

      real(8) ref_del,rcv_lat,rcv_lon,az,baz,new_del

      ! obtain new source location
      d_lat = rmodel(2)
      d_lon = rmodel(3)
      d_dep = rmodel(4)
      new_lat = theSc + d_lat
      new_lon = phiSc + d_lon
      new_dep = depth + d_dep
      theSc_geoc = geocen_deg(new_lat)

      ! obtain receiver location
      rcv_lat = geocen_deg(rcv_latlon(k,1))
      rcv_lon = rcv_latlon(k,2)
      ref_del = del_array(k) * d2r

      ! calc azimuth, new_del for each pair
      call azdel(theSc_geoc,new_lon,rcv_lat,rcv_lon,az,baz,new_del)
      new_del =  new_del * d2r
      d_del = new_del - ref_del
      new_azim = az*d2r

!      write(*,*) 'ref'
!      write(*,*) theSc, phiSc, depth
!      write(*,*) 'new'
!      write(*,*) new_lat, new_lon, new_dep
      endsubroutine calc_newloca

      real(8) function calc_halfd(scalar,Mrr,Mtt,Mpp,Mrt,Mrp,Mtp) 
      ! Half-duration calculation from the Seismic Moment
      ! according to emprical relationship btw Mo and halfd
      ! (see Ekstrom etal. 2012 or Dahlen&Tromp chap5.4.5 )
      implicit none
      real(8) M0, scalar, Mrr,Mtt,Mpp,Mrt,Mrp,Mtp
      real(8) coeff

      coeff = 1.05d0 * 10d0 ** (-8d0)

      M0 = 1/sqrt(2d0) * 10d0 ** scalar
      M0 =M0 * sqrt(Mrr**2d0+Mtt**2d0+Mpp**2d0+
     &  2d0 * (Mrt**2d0+ Mrp**2d0+ Mtp**2d0))

      coeff = 1.05d0 * 10d0 ** (-8d0)
      calc_halfd = coeff * M0 ** (1d0/3d0)

      endfunction calc_halfd

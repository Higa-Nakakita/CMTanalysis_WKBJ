      subroutine read_param(ranges,scales,nd)
      !Read model_parametr according to the settings specified in swift.in.
      implicit none
      real , intent(out) :: ranges(2,*)
      real , intent(out) :: scales(*)
      integer , intent(out) :: nd
      integer i , nmecparam
      real(4) tmp_ranges(2,20),tmp_scales(20)
      real(4) tmp_ranges2(2,20),tmp_scales2(20)
      integer tmp_numdim, numdim, ndim
      include "swfi_param.inc"
      include "wkbj_common.inc"
      include "swfi_common.inc"

      !open(1,file='model_parameter',status='old')
      !read(1,*) ndim
      !read(1,*) dt_min,    dt_max,    dt_scale
      !read(1,*) dlat_min,  dlat_max,  dlat_scale
      !read(1,*) dlon_min,  dlon_max,  dlon_scale
      !read(1,*) ddep_min,  ddep_max,  ddep_scale
      !read(1,*) Mrr_min,   Mrr_max,   Mrr_scale
      !read(1,*) Mtt_min,   Mtt_max,   Mtt_scale
      !read(1,*) Mpp_min,   Mpp_max,   Mpp_scale
      !read(1,*) Mrt_min,   Mrt_max,   Mrt_scale
      !read(1,*) Mrp_min,   Mrp_max,   Mrp_scale
      !read(1,*) Mtp_min,   Mtp_max,   Mtp_scale
      !read(1,*) scalar_min,scalar_max,scalar_scale
      !read(1,*) halfd_min, halfd_max, halfd_scale
      !close(1)

      !open(1,file='../../model_parameter',status='old')
      open(1,file='model_parameter',status='old')
      read(1,*) tmp_numdim
      do ndim = 1, tmp_numdim
         read(1,*) tmp_ranges(1,ndim),tmp_ranges(2,ndim),
     &             tmp_scales(ndim+1)
      enddo
      close(1)

      do ndim = 1, tmp_numdim
         write(*,*)  tmp_ranges(1,ndim),tmp_ranges(2,ndim),
     &             tmp_scales(ndim+1)
      enddo


      if (MT_CMT == 0) then 
        write(*,*) 'MTinv'
        if (trM==0) then
          numdim = tmp_numdim - 5
          do ndim = 1, 2
             ranges(1,ndim) = tmp_ranges(1,ndim+4)
             ranges(2,ndim) = tmp_ranges(2,ndim+4)
             scales(ndim+1) = tmp_scales(ndim+5) 
             !write(*,*) ndim, ranges(1,ndim), ranges(2,ndim)
          enddo
          do ndim = 4, numdim 
             ranges(1,ndim-1) = tmp_ranges(1,ndim+4)
             ranges(2,ndim-1) = tmp_ranges(2,ndim+4)
             scales(ndim) = tmp_scales(ndim+5) 
             !write(*,*) ndim-1, ranges(1,ndim-1), ranges(2,ndim-1)
          enddo
          numdim = numdim - 1
        elseif (trM==1) then
          numdim = tmp_numdim - 5
          do ndim = 1, numdim
             ranges(1,ndim) = tmp_ranges(1,ndim+4)
             ranges(2,ndim) = tmp_ranges(2,ndim+4)
             scales(ndim+1) = tmp_scales(ndim+5)   
             !write(*,*) ndim, ranges(1,ndim), ranges(2,ndim)
          enddo          
        endif
      elseif (MT_CMT == 1) then  ! CMT inversion
        !! TraceM = Free 
        if (trM==1) then
           if (origT==2) then ! dt & halfd: search
              numdim = tmp_numdim
              do ndim = 1, numdim
                 ranges(1,ndim) = tmp_ranges(1,ndim)   
                 ranges(2,ndim) = tmp_ranges(2,ndim)   
                 scales(ndim+1) = tmp_scales(ndim+1)   
              enddo
           elseif (origT==1) then ! dt : search
              numdim = tmp_numdim - 1
              do ndim = 1, numdim
                 ranges(1,ndim) = tmp_ranges(1,ndim)   
                 ranges(2,ndim) = tmp_ranges(2,ndim)   
                 scales(ndim+1) = tmp_scales(ndim+1)   
              enddo
           elseif (origT==0) then ! halfd : search
              numdim = tmp_numdim - 1
              do ndim = 1, numdim
                 ranges(1,ndim) = tmp_ranges(1,ndim+1)
                 ranges(2,ndim) = tmp_ranges(2,ndim+1)
                 scales(ndim+1) = tmp_scales(ndim+2)
              enddo
           else
              write(*,*) 'Something wrong in Torig?'
              stop
           endif

        elseif (trM==0) then ! traceM = 0; Mpp=-(Mrr+Mtt)
           do ndim = 1, 6 
              tmp_ranges2(1,ndim) = tmp_ranges(1,ndim)   
              tmp_ranges2(2,ndim) = tmp_ranges(2,ndim)   
              tmp_scales2(ndim+1) = tmp_scales(ndim+1)   
           enddo
           do ndim = 8,tmp_numdim ! skip Mpp
              tmp_ranges2(1,ndim-1) = tmp_ranges(1,ndim)   
              tmp_ranges2(2,ndim-1) = tmp_ranges(2,ndim)   
              tmp_scales2(ndim) = tmp_scales(ndim+1)   
           enddo
           if (origT==2) then ! dt & halfd: search
              numdim = tmp_numdim - 1
              do ndim = 1, numdim
                 ranges(1,ndim) = tmp_ranges2(1,ndim)
                 ranges(2,ndim) = tmp_ranges2(2,ndim)
                 scales(ndim+1) = tmp_scales2(ndim+1)
              enddo
           elseif (origT==1) then ! dt : search
              numdim = tmp_numdim - 2
              do ndim = 1, numdim
                 write(*,*) tmp_ranges2(1,ndim)
                 write(*,*) tmp_ranges2(2,ndim)
                 write(*,*) tmp_scales2(ndim+1)
                 write(*,*) 
                 ranges(1,ndim) = tmp_ranges2(1,ndim)
                 ranges(2,ndim) = tmp_ranges2(2,ndim)
                 scales(ndim+1) = tmp_scales2(ndim+1)
              enddo
           elseif (origT==0) then ! halfd : search
              numdim = tmp_numdim - 2
              do ndim = 1, numdim
                 ranges(1,ndim) = tmp_ranges2(1,ndim+1)
                 ranges(2,ndim) = tmp_ranges2(2,ndim+1)
                 scales(ndim+1) = tmp_scales2(ndim+2)
              enddo
           else
              write(*,*) 'Something wrong in Torig?'
              stop
           endif

        else 
           write(*,*) 'Something wrong in trM?'
           stop
        endif
      else
        write(*,*) 'Wrong in MT_CMT',MT_CMT
        stop
      endif
      
      nd = numdim      
      scales(1) = 1

      !! write out for confirmation
      open(1,file='userinit_conf.txt',position='append')
      write(1,*) 
      write(1,*) 'Inversion condition'
      if (MT_CMT == 0) then
         write(1,*) 'MT inversion (1st setp)'
      elseif (MT_CMT == 1) then
         write(1,*) 'CMT inversion (2nd setp)'
      endif
      if (trM==0) then
         write(1,*) 'trM=0 .i.e. Mpp=-(Mrr+Mtt)'
      elseif (trM==1) then
         write(1,*) 'trM->Free'
      endif
      if (MT_CMT == 1) then
         if (origT==0) then
            write(1,*) 'dt=0, halfd:Free'
         elseif (origT==1) then
            write(1,*) 'dt:Free, halfd: Dependent variable of Mo'
         else 
            write(1,*) 'dt:Free, halfd:Free'
         endif
      endif
      write(1,*) 'The number of dimension of parameter space:',nd
      write(1,*) '===== Parameters Search Range ==='
      do i = 1, nd
          write(1,*) ranges(1,i),ranges(2,i),scales(i+1)
      enddo 
      close(1)      
      !stop

      
      endsubroutine read_param


      subroutine user_init_cmtinv
      ! Read various input files and eigen-tables created by minos, and
      ! compute fix-term in this inversion.
      ! 1st step: Obtain the coordinates of stations and the initial source, 
      ! and  the paths of the observed waveforms.
      ! 2nd step: Read swfi.in.
      ! 3rd step: Read and preprocess the observed waveforms.
      ! 4th step: Read the eigenfunctions, eigenparameters files and
      ! Compute the fix-term of WKBJ formula (i.e., path-term,
      ! receiver-term). 

      implicit none
      integer i, n, k, lw, lofw
      integer U_num1, iicomp, pos, ios

!      real(4) pi,pi2,r2d,d2r,dkm2d,d2km,dkm2r <- already declared at *.inc
      include "swfi_param.inc"
      include "wkbj_common.inc"
      include "swfi_common.inc"
      include "eigen_common.inc"
      include "new_eigen_common.inc"

      integer nlen, nerr
      character(len=200) obswv_path(maxstas,numcomp)
      character(len=300) obsdata_path,rcv_coordinates,cha
      character(len=300) src_coordinates
      real(8) theSc_geoc, geocen_deg
      real(8) rlat,rlon,rlat_geoc,az,baz,del

      ! 2nd STEP
      integer U_num2, iunit
      real(8) gvmin(maxwindow),gvmax(maxwindow)
      character(len=20) swfi_in
      character(len=200) gv_table
      ! For Normalisation
!      real(8) G,RHO,RA,EM,unit_seis
!      data G,RHO,RA,EM/6.6723d-11,5.515d3,6.371d6,5.97d24/
      real(8) omgnrm,velnrm,dtmp


      ! 3rd STEP
      integer iw
      real(8) gvmin_temp, gvmax_temp, ytmin_temp, ytmax_temp
      ! read obsdata
      real(4) obstmp(maxdata),obstmp_tap(maxdata),obstmp_fil(maxdata)
      real(4) ytb_org,ytd_org,yte_org
      integer npts_org,npts_max,nntmp
      real(8) ytime_min,ytime_max
      integer npts_new_fit
      integer nstart,nend,nnl,nnr

      ! freq adjustment
      integer tmp1, tmp2, tmp1_old, tmp2_old, STmode,jtmp_sup

      ! === time_buffer ====
      t0_buffer = 300d0
      !t0_buffer = 120d0
      !t0_buffer = 60d0
      !t0_buffer = 0d0

      ! === set file names =========
      obsdata_path = "obsdata_path"
      rcv_coordinates="rcv_coordinates"
      src_coordinates="src_coordinates"
      swfi_in  = "swfi.in"
      !swfi_in  = "../swfi.in"
      !swfi_in = "../../swfi.in"

      ! === delete allparameters fils ==
      open(100,file='all_params.dat')
      close(100,status='delete')

      ! === delete some check files ==
      open(100,file='ck_eifrcv_nl.txt')
      close(100,status='delete')
      open(100,file='ck_eifsrc_nl.txt')
      close(100,status='delete')
      open(100,file='ck_radpattern.txt')
      close(100,status='delete')
      open(100,file='ck_zAexpfix.txt')
      close(100,status='delete')
      open(100,file='ck_zAexpfrq.txt')
      close(100,status='delete')

      ! === set unit number for files
      U_num1 = 10
      U_num2 = 11

      ! ==== compute constant number (wkbj_common.inc) ======
      pi  = 4d0 * datan(1d0)
      pi2 = 2d0 * pi
      r2d = 180d0 / pi
      d2r = pi / 180d0
      dkm2d = r2d / 6371d0 ! convert km  -> deg on the earth 
      d2km  = d2r * 6371d0 ! convert deg -> km  on the earth
      dkm2r = 1d0 / 6371d0 ! convert km -> rad  on the earth

      ! ===== 1st Step; Obtain information of path for dir & data ==========  
      ! ==== obtain coordinate of stations
      k = 0
      ios = 0
      open(U_num1,file=rcv_coordinates,status='old')
      do i = 1, 100000
         read(U_num1,*,iostat=ios)
         if (ios/=0) then
             exit
         endif
         k = k + 1
      enddo
      close(U_num1)
      rcvnum = k
      write(*,*) 'rcvnum',rcvnum

      open(U_num1,file=rcv_coordinates,status="old")
      do k = 1, rcvnum
         read(U_num1,*) rcvnm(k), rcv_latlon(k,1),rcv_latlon(k,2)
         write(*,*) rcvnm(k), rcv_latlon(k,1),rcv_latlon(k,2)
      enddo
      close(U_num1)



      ! ===== obtain path of obs_waveforms data ====
      !write(*,*) rcvnum
      open(U_num1,file=obsdata_path,status="old")
      do k = 1, rcvnum
         do iicomp = 1, 3 !1:Z, 2:R, 3:T
           read(U_num1,'(a)') cha

           ! check whether input comp is correct or not
           pos = 1 ! dummy
           if (iicomp == 1) then
                   pos = index(cha,"LHZ")
                   if (pos == 0) then
                        write(*,*) "WRONG in comp Z at",k
                        stop
                   endif
           elseif (iicomp == 2) then
                   pos = index(cha,"LHR")
                   if (pos == 0) then
                        write(*,*) "WRONG in comp R at",k
                        stop
                   endif
           elseif (iicomp == 3) then
                   pos = index(cha,"LHT")
                   if (pos == 0) then
                        write(*,*) "WRONG in comp T at",k
                        stop
                   endif
           endif

           lw = lofw(cha)
           obswv_path(k,iicomp) = cha(1:lw)
           !write(*,*) k,iicomp,lw,obswv_path(k,iicomp)
         enddo
      enddo
      close(U_num1)
  

      ! ==== obtain coordinates of source (initial)
      open(U_num1,file=src_coordinates,status="old")
      read(U_num1,*) theSc,phiSc,depth
      close(U_num1)
      ! --- convert geographic coordinates to geocentric
      theSc_geoc = geocen_deg(theSc)

      ! === compute azimuth and distance for each path
      do k = 1, rcvnum
         rlon = rcv_latlon(k,2)
         rlat = rcv_latlon(k,1)
         rlat_geoc = geocen_deg(rlat)
         call azdel(theSc_geoc,phiSc,rlat_geoc,rlon,az,baz,del)
         az_array(k) = az
         baz_array(k)= baz
         del_array(k)= del
         distance(k) = del * d2km
      enddo
      ! ============ 1st STEP: OVER  ========================



      ! =========== 2nd STEP: READ IN "swfi.in" ===========
      open(U_num2,file=trim(swfi_in),status="old") 
      read(U_num2,*) 
      read(U_num2,*)
      read(U_num2,*) nn1,nn2  ! min,max mode branches for calculation
      read(U_num2,*) iobstype ! input wave type (1:Disp, 2:Vel, 3:Acc)
      read(U_num2,*) iunit    ! unit of input data (1:m, 2:mm,3:um, 4:nm) 

      ! --- Time window based on group velocity 
      read(U_num2,*) nwave  ! the number of time windows
      do i = 1, nwave
         read(U_num2,*) gvmin(i), gvmax(i)  ! min, max group velocity (km/s)
         read(U_num2,*) weight_data(i)      ! weight factor for this time window
         read(U_num2,*) bp1(i), bp2(i)      ! frequency band (mHz)
         bpfrq_min(i) = bp1(i) / 1000d0   ! (Hz)
         bpfrq_max(i) = bp2(i) / 1000d0
      enddo
      ! skip rows when nwave < maxwidnow=4
      if (nwave < maxwindow) then
         do i = 1, maxwindow - nwave
            read(U_num2,*)
            read(U_num2,*)
            read(U_num2,*)
         enddo
      endif

      ! --- for some info for inversion
      read(U_num2,*) Z_wei,R_wei,T_wei !weighting factor for each comp
      read(U_num2,*) ZRT_ratio  !weighting factor based on area ratio of each comp
      read(U_num2,*) geomet_wei !geometrical attenuation correction of misfit (1:ON,0:OFF)
      read(U_num2,*) rnrm       ! norm index i (L_i norm)
      read(U_num2,*) MT_CMT     !Flag of inversion flow(0:MT inv, 1:CMT inv)
      read(U_num2,*) trM        !(trM=0; 0:ON, 1:OFF)
      read(U_num2,*) origT      !(0: dt=0, 1:dt=freeparam)
      read(U_num2,*) conf_rcv   !(0:OFF, 1:ON)
      read(U_num2,*) conf_path
      read(U_num2,*) conf_fix
      read(U_num2,*) conf_src
      read(U_num2,*) iso_ani !iso=0, ani=1
      read(U_num2,*) which_struc !(0:PREM, 1:Yoshizawa2010, 2:Yoshizawa2010+jma2001)
      read(U_num2,*) src_iso_ani !prem iso=0, ani=1
      read(U_num2,*) cut_cost1, prd_inf   !0:OFF 1:ON
      read(U_num2,*) mode_skip   
      close(U_num2)
      !stop
      ! ==== 2nd STEP: OVER =======================

      ! Normalisation factors
      if(iunit==2) then
              unit_seis = 1d3
      elseif(iunit==3) then
              unit_seis = 1d6
      elseif(iunit==4) then
              unit_seis = 1d9
      else
              unit_seis = 1d0
      endif
      omgnrm = dsqrt(pi*G*RHO)
      velnrm = RA * omgnrm
      dtmp = G * RHO * EM * RA * (3d0/4d0)


      ! ~~~~~ Local modes 0409 ~~~~~~~~~~~~~~~~~~~~~~
      call user_init_localmodes
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      

      ! ==== 3rd STEP =======================
      ! READ obsdatas in SAC format, SET time-windows, taper, filter
      ! do loop for each rcv, comp

      ! set TEMPORAL group velocity (gv_temp)
      !gvmin_temp = min(gvmin(1:nwave)
      !gvmax_temp = max(gvmax(1:nwave)
      gvmin_temp = gvmin(1)
      gvmax_temp = gvmax(1)
      if (nwave > 1) then
         do iw = 1, nwave
            if(gvmin_temp < gvmin(iw)) gvmin_temp = gvmin(iw)
            if(gvmax_temp > gvmax(iw)) gvmax_temp = gvmax(iw)
         enddo
      endif


      do k = 1, rcvnum
          write(*,*)
          write(*,*) '====== rcv number:',k,'======='
          ! calc temporal timewindow based on gvmin, gvmax
          ytmin_temp = distance(k)/gvmax_temp
          ytmax_temp = distance(k)/gvmin_temp

          do iicomp = 1, 3
             write(*,*)
             if (iicomp.eq.1) then
                write(*,*) '---- comp: Vertical ---- '
             elseif (iicomp.eq.2) then
                write(*,*) '---- comp: Radius   ----'
             elseif (iicomp.eq.3) then
                write(*,*) '---- comp: Transverse ---- '
             else
                write(*,*) 'Error: Something wrong in iicomp',iicomp
                stop
             endif

             ! ---- read obsdata for each component 
             cha = obswv_path(k,iicomp)
             lw = lofw(cha)
          !        call readdata(obstmp,ytb_org,ytd_org,ydpmax,ydpmin,
          !&       npts_org,cha(1:lw),lw)
             call rsac1(trim(cha(1:lw)),obstmp,npts_org,ytb_org,
     &       ytd_org,maxdata,nerr) !sac func not object file (2025, 4/2)
             yte_org = ytb_org + ytd_org * real(npts_org)


             if(ytmin_temp<ytb_org .or. ytmax_temp>yte_org) then
                     write(*,*) 'ERROR: timewindow of obs is too short!'
                     write(*,*) "ytmin_temp,ytb_org,ytmax_temp,yte_org"
                     write(*,*) ytmin_temp,ytb_org,ytmax_temp,yte_org
                     stop
             endif

             ! (TEMP) Calc nyquist freq, freq interval for FFT
             ! I think (ntmax,nyquist,dfrq,domega) donot depend on rcv 
             ! as long as (npts_max < 2028) is satisfied
             npts_max = nint(ytmax_temp/ytd_org) + 1
             do n = 11, 15
                nntmp = 2 ** n
                if(nntmp>npts_max) then
                        ntmax(k) = nntmp
                        exit
                endif
             enddo
             nyquist(k) = ntmax(k) / 2 + 1
             dfrq(k) = 1d0 / (dble(ntmax(k))*ytd_org)
             domega(k) = pi2 * dfrq(k)
             write(*,*) "ntmax,nyquist,dfrq,domega"
             write(*,*) ntmax(k),nyquist(k),dfrq(k),domega(k)

             !open(100,file='ck_obstmp.txt',position='append')
             !write(100,*) 'k,iicomp',k,iicomp
             !do i = 1, ntmax
             !   write(100,*) i, obstmp(i)
             !enddo
             !close(100)

             ! set time-window for each path correspond to gv
             do iw = 1, nwave
                ytime_min = distance(k) / gvmax(iw) + t0_buffer
                ytime_max = distance(k) / gvmin(iw) + t0_buffer
                npts_new_fit = nint((ytime_max-ytime_min)/ytd_org) + 1
                if (npts_new_fit > npts_org) then
                    write(*,*) 'ERROR: input timewindow is too short'
                    write(*,*) 'Must be: npts_new_fit < npts_org'
                    write(*,*) npts_new_fit < npts_org
                    stop
                else
                    n_data(k,iw) = npts_new_fit
                endif
                nstart = nint((ytime_min-ytb_org)/ytd_org) + 1
                nend   = nstart + n_data(k,iw) - 1
                nnl    = int(n_data(k,iw) * 0.05) !Left side teper portion(5%) 
                nnr    = nnl
                call cos_tap(obstmp,obstmp_tap,nstart,nend,nnl,nnr)

                t_begin(k,iw) = ytb_org + (nstart-1) * ytd_org
                t_end(k,iw)   = t_begin(k,iw) + ytd_org * (n_data(k,iw))
                tint = ytd_org

                ! the number of points for FFT
                do n = 6, 20
                   nntmp = 2 ** n
                   if(nntmp>n_data(k,iw)) then
                           n_fft(k,iw) = nntmp
                           exit
                   endif
                enddo
                call bandpass_fil(obstmp_tap,obstmp_fil,iw,k)

                observed_data(1:n_data(k,iw),iicomp,iw,k) = 
     &          obstmp_fil(1:n_data(k,iw))

                call w_seis(obstmp_fil,iw,0,k,iicomp)
                enddo  ! enddo nwave
           enddo !enddo iicomp
      enddo  ! enddo rcvnum
      !stop
      ! ==== 3rd STEP: OVER  =======================


      ! Just write out summary info of User_init 
      write(*,*) 'Attention'
      write(*,*) 'This is based on Right-hand rule (Dahlen&Tromp p365)'
      write(*,*) 'conversion Right -> Left will be applied (at syntmp)'
      open(1,file="userinit_conf.txt")
      write(1,*) 'based on Right-hand rule (Dahlen&Tromp p365)'
      write(1,*) 'conversion Right -> Left will be applied (at syntmp)'
      write(1,*)
      write(1,*) 'initial source: ',real(theSc),real(phiSc),real(depth)
      do k = 1, rcvnum
         write(1,*) '== rcv info (rcv:',trim(rcvnm(k)),' num:',k,') ==='
         write(1,*) 'rcv: ',real(rcv_latlon(k,1)),real(rcv_latlon(k,2))
         write(1,*) 'epicental distance(km): ',real(distance(k))
         write(1,*) 'azimuth(deg):',real(az_array(k)),'baz(deg):',
     &    real(baz_array(k))
         do iw = 1, nwave
           write(1,*) '--- time window No:',iw,' ---'
           write(1,*) 'gvmin gvmax (km/s):',real(gvmin(iw)),
     &    real(gvmax(iw))
           write(1,*) 'weight:',real(weight_data(iw))
           write(1,*) 'band pass (mHz):',real(bp1(iw)),real(bp2(iw))
           write(1,*) 'ytd=',real(ytd_org),' num of points',
     &    real(n_data(k,iw))
           write(1,*) ' t_begin ~ t_end (sec)'
           write(1,*) real(t_begin(k,iw)),'--',real(t_end(k,iw))
           write(1,*)
           write(1,*) 'ntmax:',ntmax(k),' nyquist:',nyquist(k)
           write(1,*) 'dfrq:',real(dfrq(k)),' domega:',real(domega(k))
         enddo
       enddo
       write(1,*) 't0_buffer:',t0_buffer
       write(1,*)
       close(1)

       ! ======== 4th Step ========== 
       ! Read eiegnparameters and eigen functions
       do k = 1, rcvnum
          do iicomp = 1, numcomp
             call rd_eigen_minos(iicomp,k)
          enddo
       enddo
       !stop

       ! Calculate fixterm of WKBJ formula.
       ! i.e., receiver term and path term.
       do k = 1, rcvnum
          do iicomp = 1, numcomp
             call fixterm_calc(iicomp,k)
          enddo
       enddo
       !stop

       ! READ Eigen function of the source (PREM) 
       call rd_eif_src
       ! ======== 4th Step: OVER ==========  


       ! Just performing some adjustment for freq range.
       ! During spectral-domain synthesis, 
       ! it is necessary to match the frequency range.
       ! 2025 0717
       do STmode = 1, 2
         do n = nn1, nn2
            tmp1_old = 0
            tmp2_old = 500
            do k = 1, rcvnum
              tmp1 = j_skip1(n,STmode,k)
              tmp2 = j_skip2(n,STmode,k)
              if (tmp1>tmp1_old) then
                  j_skip1_adjst(n,STmode) = tmp1
                  tmp1_old = tmp1
              endif
              if (tmp2<tmp2_old) then
                  j_skip2_adjst(n,STmode) = tmp2
                  tmp2_old = tmp2
              endif
            enddo
         enddo
       enddo


       ! j_skip2_adjst stores the idx of freq that considering
       ! for local modes at each station.
       ! j_frq_path_max stores the max of idx of freq of the mode 
       ! for perturbing the source location.
       ! When the waveform calculation at each forward modelling,
       ! it is required that j_skip2_adjst = < j_frq_path_max, 
       ! but this condition is sometimes not satisfied.
       ! below code will coorect j_skip2_adjst in that case.
       !2025 1129
       open(1,file='correct_jjtmp_record.txt')
       do STmode = 1,2
         do n = nn1, nn2
            jtmp_sup = j_skip2_adjst(n,STmode)
            if (j_frq(n,jtmp_sup,STmode)>j_frq_path_max(n,STmode)) then
                write(1,*) n,STmode,j_frq(n,jtmp_sup,STmode),
     &                     j_frq_path_max(n,STmode)
                j_skip2_adjst(n,STmode) = j_frq_path_max(n,STmode)
            endif
         enddo
       enddo
       close(1)



       count_ns = 0

       !write(*,*) 'User init stop'
       !stop
           
      endsubroutine user_init_cmtinv



      subroutine user_init_localmodes
      ! Obtain path to eigenparameters file for local mode of each grid
      ! and ray-path.(make array "cha_div_deig" from "divdel_info").
      ! Create the array of the epicentral distances "gcp_dels" for each
      ! ray-path (used for ray-path integration). 
      ! If there is no eigenparameters file for the certain grid,
      ! that local mode will be for PREM.
      implicit none
      include "swfi_param.inc"
      include "wkbj_common.inc"
      include "swfi_common.inc"
      include "new_eigen_common.inc"

      integer i, n, k, ios, alldivnum
      integer U_num1, U_num2, U_num3
      real(8) lat, lon
      character(len=20) rcv_coordinates,src_coordinate,div_coordinates
      character(len=200) cha_tmp1,cha_tmp1_1,cha_tmp1_2,cha_tmp1_3
      character(len=200) cha_tmp2
      character(len=200) cha_tmp3,cha_tmp4
      character(len=4) cha_lat
      character(len=5) cha_lon

      integer ii, idx
      real(8) clat1(10000),clon1(10000),clat2(10000),clon2(10000)
      real(8) ref_lat(10000), ref_lon(10000)
      real(8) reflat, reflon, slat, slon, rlat, rlon
      real(8) d_dels(10000),d_dists(10000),sum_dist(10000)
      real(8) d_del_array(maxdivs,maxstas),d_dist_array(maxdivs,maxstas)
      real(8) sum_dist_array(maxdivs,maxstas)
      real(8) divs_array(maxdivs,maxstas), sumdel
      integer del_idx_array(10000)
      integer rcv_div_idxs(maxstas)
      integer Flag_src, Flag_rcv, idx_beg, idx_end
      integer div_idx_array(2,maxdivs)
      integer access, fil_ck1, fil_ck2
      character(len=200) ani_prem_path1, ani_prem_path2, ani_prem_path
      character(len=200) iso_prem_path1, iso_prem_path2, iso_prem_path


      ! === delete allparameters fils ==
      open(100,file='ref_latlon_out.txt')
      close(100,status='delete')
      open(100,file='ref_eifpath.txt')
      close(100,status='delete')
      open(100,file='modes_range_info.txt')
      close(100,status='delete')
      open(100,file='ck_local_eigp.txt')
      close(100,status='delete')
      open(100,file='ck_local_integ.txt')
      close(100,status='delete')
      open(100,file='remedy_eigps.txt')
      close(100,status='delete')


      rcv_coordinates = 'rcv_coordinates'
      src_coordinate  = 'src_coordinate'
      div_coordinates = 'div_coordinates'


      !====================================================
      !======= Path to eigen-relation Setting =============
      ani_prem_path1 = '/work94A/higa/minos_Eigenfunc'
      ani_prem_path2 = '/noocean_prem_1016/calc_Eigen_1016/'//
     &   'execute_modi_rcvnl_251216'
      ani_prem_path = trim(ani_prem_path1)//trim(ani_prem_path2)

      iso_prem_path1 = '/work94A/higa/minos_Eigenfunc'
      iso_prem_path2 = '/noocean_prem_20250207/'//
     &   'execute_iso_modi_rcvnl_251216'
      iso_prem_path = trim(iso_prem_path1)//trim(iso_prem_path2)

      ! Must set path to local mode path,
      ! to do below procedure regardless of assumed model
      cha_tmp1_1 = '/work/higa/calc_eig_3D/Yoshizawa2010/'
      if (iso_ani==0) then !iso
           prem_path = trim(iso_prem_path)
           cha_tmp1_2 = '0422_iso/outdir/iso/'
      elseif (iso_ani==1) then !ani
           prem_path = trim(ani_prem_path)
           cha_tmp1_2 = '0422_ani/outdir/aniso/'
      else
           write(*,*) 'wrong in iso_ani',iso_ani
           stop
      endif

      if (which_struc==0) then
            ! it will be regarded as out-off grid for all grid,
            ! then PREM path will be inserted in cha_div_eig(idx,k)
            write(*,*) 'path: PREM 1D'
      elseif (which_struc==1) then
            write(*,*) 'path: Yoshizawa2010 3D model'
      elseif (which_struc==2) then
            write(*,*) 'path: Yoshizawa2010 3D & JMA2001 model'

            ! path to local dir must be over write for jma2001&kazu case
            cha_tmp1_1 = '/work/higa/calc_eig_3D/jma_Yoshizawa2010/'
            if (iso_ani==0) then !iso
                 cha_tmp1_2 = '1116_iso/outdir/iso/'
            elseif (iso_ani==1) then !ani
                 cha_tmp1_2 = '1116_ani/outdir/aniso/'
            endif
      else
            write(*,*) 'Wrong in which_struc',which_struc
            stop
      endif
      !====================================================
      !======Over: Path to eigen-relation Setting =============

      
      U_num1 = 20
      U_num2 = 21
      U_num3 = 22

      ! Read grids information to obtain reference coordinates and
      ! and designate the the path to directory. 
      ! Furtermore obtain the d_dels info.
      k = 0
      open(U_num3,file='div_del_info',
     &            status='old',action='read')
      do i = 1, 10000
         read(U_num3,*,iostat=ios)
         if (ios/=0) then
             exit
         endif
         k = k + 1
      enddo
      close(U_num3)
      alldivnum = k
      open(U_num3,file='div_del_info',
     &            status='old',action='read')
      do i = 1, alldivnum
         read(U_num3,*) clat1(i),clon1(i),clat2(i),clon2(i),ref_lat(i),
     &   ref_lon(i),d_dels(i),d_dists(i),sum_dist(i),del_idx_array(i)
      enddo
      close(U_num3)

      ! Find eigen-relation's directory path
      ! receiver
      do k = 1, rcvnum
         lat = rcv_latlon(k,1)
         lon = rcv_latlon(k,2)
         do ii = 1, alldivnum
         !write(*,*) lat,clat2(ii),lon,clon2(ii)
            if(lat==clat2(ii) .and. lon==clon2(ii)) then
                reflat = ref_lat(ii)
                reflon = ref_lon(ii)
                write(cha_lat,'(f4.1)') reflat
                write(cha_lon,'(f5.1)') reflon
                !cha_tmp1 = 'rcvs/'//'N'//cha_lat//'/'
                !cha_tmp2 = cha_lat//'.'//cha_lon
                !cha_rcv_eif(k) = trim(cha_tmp1)//trim(cha_tmp2)
                !cha_tmp1 = '/work94A/higa/minos_Eigenfunc/'
                !cha_tmp2= 'noocean_prem_1016/calc_Eigen_1016/execute'
                !cha_rcv_eif(k) = trim(cha_tmp1)//trim(cha_tmp2)
               ! write(*,*) cha_rcv_eif(k)
                !rcv_div_idxs(k) = del_idx_array(ii)
        !           write(*,*) rcv_div_idxs(k), rcvnm(k),
        !&                     lat,reflat, lon,reflon
                exit
            endif
         enddo
      enddo

      ! source
      slat = theSc
      slon = phiSc
      do ii = 1, alldivnum
         if (slat==clat1(ii) .and. slon==clon1(ii)) then
              reflat = ref_lat(ii)
              reflon = ref_lon(ii)
              !write(*,*) ii,  lat,reflat, lon,reflon
              write(cha_lat,'(f4.1)') reflat
              write(cha_lon,'(f5.1)') reflon
              !cha_tmp1 = 'src/'//'N'//cha_lat//'/'
              !cha_tmp2 = cha_lat//'.'//cha_lon
              !cha_src_eif = trim(cha_tmp1)//trim(cha_tmp2)
              !cha_tmp1 = '/work94A/higa/minos_Eigenfunc/'
              !cha_tmp2= 'noocean_prem_1016/calc_Eigen_1016/execute'
              !cha_src_eif = trim(cha_tmp1)//trim(cha_tmp2)
              !write(*,*) cha_rcv_eif(k)
              !write(*,*) 'src'
              !write(*,*) cha_src_eif
         endif
      enddo

      ! Make array of div
      !k = 1
      !do ii = 1, alldivnum
      !  Flag_rcv = 0
      !  if (slat==clat1(ii) .and. slon==clon1(ii)) then
      !       Flag_src = 1
      !       idx_beg = ii
      !  endif
      !  if (Flag_src==1 .and. Flag_rcv==0) then
      !       rlat = rcv_latlon(k,1)
      !       rlon = rcv_latlon(k,2)
      !       if(rlat==clat2(ii) .and. rlon==clon2(ii)) then
      !           idx_end = ii
      !          ! write(*,*) idx_beg,"~", idx_end
      !           div_idx_array(1,k) = idx_beg
      !           div_idx_array(2,k) = idx_end
      !           ! initiallize
      !           Flag_rcv = 1
      !           Flag_src = 0
      !           k = k + 1
      !           idx_beg = 0
      !           cycle
      !       endif
      !   endif
      !enddo


      ! Mae array of div (20251208)
      ! can handle cases when rcv_coordinates is not in alphabetical order.
      do k = 1, rcvnum
        Flag_src = 0
        Flag_rcv = 0
        do ii = 1, alldivnum
          if (slat==clat1(ii) .and. slon==clon1(ii)) then
            idx_beg = ii
            Flag_src = 1
          endif
          if (Flag_src==1 .and. rcv_latlon(k,1)==clat2(ii)
     &                    .and. rcv_latlon(k,2)==clon2(ii)) then
            idx_end = ii
            div_idx_array(1,k) = idx_beg
            div_idx_array(2,k) = idx_end
            exit
          endif
        enddo
      enddo




      do k = 1, rcvnum
         write(*,*) rcvnm(k),div_idx_array(1,k),'~', div_idx_array(2,k)
         if (div_idx_array(1,k)<=0 .or.  div_idx_array(2,k)<=0) then
                 write(*,*) "stop @user_init_localmodes"
                 write(*,*) "idx<=0 is detected, check div_del_info"
                 stop
         endif
         idx_beg = div_idx_array(1,k)
         idx_end = div_idx_array(2,k)
         idx = 1
         sumdel = 0d0
         do ii = idx_beg, idx_end
      !       write(*,*) ref_lat(ii),ref_lon(ii),
      !&                 d_dels(ii),d_dists(ii),sum_dist(ii)
            write(cha_lat,'(f4.1)') ref_lat(ii)
            write(cha_lon,'(f5.1)') ref_lon(ii)
            !cha_tmp1_1 = '/home/higa/Experimatal_Site/Apr_2025/'
            !cha_tmp1_2 = 'relfect_3Deffect/1_Stage/0419_2/outdir/iso/' 
            !cha_tmp1_1 = '/work/higa/calc_eig_3D/Yoshizawa2010/'
            !cha_tmp1_2 = '0422_iso/outdir/iso/'
            cha_tmp1_3 = 'N'//cha_lat//'/' 
         cha_tmp1 = trim(cha_tmp1_1)//trim(cha_tmp1_2)//trim(cha_tmp1_3)
            cha_tmp2 = cha_lat//'_'//cha_lon

            cha_tmp3 = trim(cha_tmp1)//trim(cha_tmp2)//'/eigen_pathav_S'
            cha_tmp4 = trim(cha_tmp1)//trim(cha_tmp2)//'/eigen_pathav_T'
            fil_ck1 = access(cha_tmp3,' ')
            fil_ck2 = access(cha_tmp4,' ')
            !write(*,*) cha_tmp3 
            !write(*,*) cha_tmp4

            ! In case of out of all-grids or PREM ++++++++++++++++++
            if (fil_ck1 /=0 .or. fil_ck2 /=0 .or. which_struc==0) then
                !write(*,*) 'Not found',cha_tmp4,' or ',cha_tmp3
                cha_div_eig(idx,k) = trim(prem_path)
            else
                cha_div_eig(idx,k) = trim(cha_tmp1)//trim(cha_tmp2)
            endif
            ! +++++++++++++++++++++++++++++++++++++++++++++++++++++

            write(*,*) cha_div_eig(idx,k)
            !cha_tmp1 = '/work94A/higa/minos_Eigenfunc/'
            !cha_tmp2= 'noocean_prem_1016/calc_Eigen_1016/execute'
            !cha_div_eig(idx,k) = trim(cha_tmp1)//trim(cha_tmp2)
            !write(*,*) k,idx,cha_div_eig(idx,k)
            d_del_array(idx,k)     = d_dels(ii)
            d_dist_array(idx,k)    = d_dists(ii)
            sum_dist_array(idx,k)  = sum_dist(ii)
            sumdel = sumdel + d_dels(ii) 
            gcp_dels(idx,k) = sumdel
      !      write(*,*) k,idx,trim(cha_div_eig(idx,k)),
       !&                 sum_dist_array(idx,k)
            idx = idx + 1
         enddo
         div_num_array(k) = idx - 1 ! remove last sum
      enddo
      !stop

      ! Just write out record 
      open(1,file='ref_latlon_out.txt')
      do k = 1, rcvnum
         write(1,*) k, rcvnm(k)
         idx = 1
         do ii = div_idx_array(1,k), div_idx_array(2,k)
            write(1,*) idx,ref_lon(ii), ref_lat(ii), d_del_array(idx,k),
     &               gcp_dels(idx,k),sum_dist_array(idx,k)
            idx = idx + 1
         enddo
       enddo
       close(1) 

       open(1,file='ref_eifpath.txt')
       do k = 1, rcvnum
          write(1,*) k,rcvnm(k)
          do idx = 1, div_num_array(k)
            write(1,*) idx,cha_div_eig(idx,k)
          enddo
       enddo
       close(1)

      ! stop
            


      endsubroutine

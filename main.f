
      Program main
      ! Just calling Neighborhood Algorithm
      ! NA consits of 3 programs
      ! user_init, forward, writemodels
      implicit none
      call na
      stop
      endprogram main

      subroutine user_init(nd,ranges,scales)
      implicit none
      integer i, n, nd
      real(4) ranges(2,*), scales(*)
 
      call user_init_cmtinv
      call read_param(ranges,scales,nd)
      !stop
      endsubroutine user_init

      subroutine forward(nd,nitr,models,misfitval,misfitvalWAV)
      implicit none
      integer i, n, nd, nitr
      real(4) models(nd), misfitval, misfitvalWAV
      real(4) t1,t2,t3

      !call cpu_time(t1)
      call forward_modelling(models,nd,nitr)
      call misfit_func(nd,misfitval,misfitvalWAV)
      !call cpu_time(t3)
      !write(*,*) t3,t1,'CPU time for forward=',t3-t1
      !write(*,*) t3,t1,'CPU time for fixterm =',t1
      stop
      endsubroutine forward

      subroutine writemodels(nd,ntot,models,misfit,misfitWAV,
     &           ns1,ns2,itrmax,nh_max,nh,header)
      implicit none
      include "swfi_param.inc"
      include "wkbj_common.inc"
      include "swfi_common.inc"
      include "eigen_common.inc"

      ! Necessary for write out SAC
      integer k, icomp, iw, npts
      real(4) syntmp_fil(maxdata), restmp(maxdata)

      ! Necessary for NA
      integer nd, ntot, itrmax, ns1, ns2, nh, nh_max
      character(*) header
      real(4) models(nd,*), misfit(ntot), mfitmin, mfitminc, mfitmean
      real(4) misfitWAV(ntot), mfitminWAV, mfitmincWAV, mfitmeanWAV
      integer lu_mod, mopt, np, ns, it, i, jj

      ! Variance Reduction
      real(4) vr_tmp1,vr_tmp2,misfit_tmp

      
      mfitmin = misfit(1)
      if (mfitmin/=mfitmin) then !Avoid NaN
          mfitmin = misfit(2)
      elseif (mfitmin/=mfitmin) then
          mfitmin = misfit(3)
      elseif (mfitmin/=mfitmin) then
          mfitmin = misfit(4)
      endif


      ns = ns1
      np = 0
      mopt = 1
      lu_mod = 14

      open(lu_mod,file='swfi.out',status='unknown')
      write(lu_mod,*) ns1,' Number of samples in starting pool'
      write(lu_mod,*) ns2,' Number of new samples per iteration'
      write(lu_mod,*) itrmax,' Number of iterations'

      !open(100,file='realmodels.dat')
      open(101,file='all_params.dat')
      ! loop over iterartions
      do it = 1, itrmax + 1
         mfitminc = misfit(np+1)
         mfitmean = 0.0
         mfitmincWAV = misfitWAV(np+1)
         mfitmeanWAV = 0.0
         write(*,*) "mfitminc,mfitmincWAV"
         write(*,*) mfitminc,mfitmincWAV
         ! find minimun and mean misfit
         do i = 1, ns
            jj = np + i
            write(*,*) 'it,i,np,jj',it,i,np,jj
            write(*,*) 'misfit: ',misfit(jj)
            if (misfit(jj)<mfitmin) then
               mfitmin = misfit(jj)
               mfitminWAV = misfitWAV(jj)
               mopt = jj
            endif
            mfitminc = min(mfitminc,misfit(jj))
            mfitmean = mfitmean + misfit(jj)
            mfitmincWAV = min(mfitminWAV,misfitWAV(jj))
            mfitmeanWAV = mfitmeanWAV + misfitWAV(jj)

      !      if (trM==1) then !trM free nd=12
      !         write(100,*) models(1:nd,jj), misfit(jj), misfitWAV(jj)
      !      elseif (trM==0) then !trM=0 nd=11
      !         write(100,*) models(1:6,jj),-(models(5,jj)+models(6,jj)),
      !&                      models(7:nd,jj), misfit(jj), misfitWAV(jj)
      !      else 
      !         write(*,*) 'Something wrong in writemodels'
      !      endif

            write(101,*) all_params(1:12,jj),misfit(jj),VR_array(jj),
     &  VR_array_cmp(jj,1),VR_array_cmp(jj,2),VR_array_cmp(jj,3)
         enddo
         mfitmean = mfitmean / ns
         mfitmeanWAV = mfitmeanWAV / ns

!         if (summary) then
!            write(lu_mod,801) it-1,mfitmin,mfitmean,mfitminc,
!     &  mfitminWAV, mfitmeanWAV, mfitmincWAV
!            do i = 1, ns
!              jj = np + 1
!            enddo
!         endif

         np = np + ns
         ns = ns2
      enddo
      !close(100) ! realmodels.dat
      close(101) ! all_params.dat


      ! Variance Reductions info
      !ns = ns1
      !open(102,file='VRs_info.txt')
      !np = 0
      !do it = 1, itrmax + 1
      !   do i = 1, ns
      !      jj = np + i
      !      do k = 1, rcvnum
      !         do icomp = 1, 3
      !            do iw = 1, nwave
      !               vr_tmp1 = VRs_array(iw,icomp,k,jj) !invalid VR
      !               vr_tmp2 = VR_array(jj) !Valid VR
      !               misfit_tmp = miss_array(iw,icomp,k,jj)
      !               write(102,*) it,jj,k,icomp,iw,vr_tmp1,vr_tmp2,
      !&                            misfit_tmp
      !            enddo
      !         enddo
      !      enddo
      !   enddo
      !   np = np + ns
      !   ns = ns2
      !enddo
      !close(102)



      nh = 24
      if(nh>nh_max) then
         write(*,*)
         write(*,*) ' Error - header array too small'
         write(*,*) ' current size  = ', nh_max
         write(*,*) ' required size = ', nh
         write(*,*)
         stop
      endif

      write(*,*) 'header(1:24)',header(1:24)
      write(header(1:24),fmt='(4i6)') ns1,ns2,itrmax,nd
 801  format('iteration:',i5,' misfit: min=',g12.6,
     & ',mean=',g12.6,', minc=',g12.6,
     & ',minWAV=',g12.6,' ,meanWAV=',g12.6,' mincWAV',g12.6)
      close(lu_mod)


      ! ==============================================
      ! write out synthetic waveform for optimal model
      write(*,*) 'mopt: ', mopt
      write(*,*) 'optimal model: ',models(1:nd,mopt)
      call forward_modelling(models(1:nd,mopt),nd,mopt)
      do k = 1, rcvnum
         do icomp = 1, 3 
           do iw = 1, nwave
              npts = n_data(k,iw) 
              !write(*,*) npts,icomp,iw,k
              !write(*,*) synthetic_data(1:npts,icomp,iw,k)
              syntmp_fil(1:npts) = synthetic_data(1:npts,icomp,iw,k) 
              call w_seis(syntmp_fil,iw,1,k,icomp,mopt)

              restmp(1:npts) = Res_data(1:npts,icomp,iw,k)
              call w_seis(restmp,iw,2,k,icomp,mopt)

           enddo
         enddo
      enddo
      ! =============================================

      endsubroutine


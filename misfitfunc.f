
      subroutine misfit_func(nd,rmisfitALL,rmisfitWAV)
      implicit none
      include "swfi_param.inc"
      include "wkbj_common.inc"
      include "swfi_common.inc"
      integer i, n, nd, k, icomp, iw, npts, npt
      real(4) rmisfitALL, rmisfitWAV
      real(4) tmpmis, tmprmis, rmis
      real(4) restmp(maxdata), weightval, weightsum, wei
      real(4) Resval
      real(4) vr_tmp1,vr_tmp2,vr_tmp3,vr_tmp4,vr_tmp5
      real(4) vr_tmp1_sum, vr_tmp2_sum
      real(4) correct_dist
      real(4) wei_cmp
      real(4) max_Z,max_R,max_T,max_now,max_cmp
      real(4) misfit_tmp
      real(4) vr_z1,vr_r1,vr_t1,vr_z2,vr_r2,vr_t2

      count_ns = count_ns + 1

      tmpmis = 0d0
      rmisfitALL = 0d0
      rmisfitWAV = 0d0
      rmis = 0d0
      wei = 0d0
      weightsum = 0d0
      do iw = 1, nwave
        weightsum = weightsum + real(weight_data(iw))
      enddo

      vr_tmp4 = 0d0
      vr_tmp1_sum = 0d0
      vr_tmp2_sum = 0d0
      vr_z1 = 0d0
      vr_r1 = 0d0
      vr_t1 = 0d0
      vr_z2 = 0d0
      vr_r2 = 0d0
      vr_t2 = 0d0
      do k = 1, rcvnum
         do icomp = 1, 3
            do iw = 1, nwave
               npts = n_data(k,iw)
               vr_tmp1 = 0d0
               vr_tmp2 = 0d0
               vr_tmp3 = 0d0

               Res_data(1:npts,icomp,iw,k) = 
     &         abs(synthetic_data(1:npts,icomp,iw,k) - 
     &         observed_data(1:npts,icomp,iw,k))

               vr_tmp1 = 
     &         sum((Res_data(1:npts,icomp,iw,k))**2d0) * tint
               vr_tmp2 = 
     &         sum(abs(observed_data(1:npts,icomp,iw,k))**2d0) * tint

               ! write out residue
               restmp(1:npts) = Res_data(1:npts,icomp,iw,k)
               !call w_seis(restmp,iw,2,k,icomp)

               ! old VR calculation: incorrect!! 
               ! calc Variance Reduction for each iw, comp, rcvs
               vr_tmp3 = (1d0 - (vr_tmp1 / vr_tmp2)) * 100d0
               VRs_array(iw,icomp,k,count_ns) = vr_tmp3

               ! New VR calculation (see Yamaya. 2025) 
               vr_tmp1_sum = vr_tmp1_sum + vr_tmp1
               vr_tmp2_sum = vr_tmp2_sum + vr_tmp2

               miss_array(iw,icomp,k,count_ns) = 
     &         sum(Res_data(1:npts,icomp,iw,k))  

               vr_tmp4 = vr_tmp4 + vr_tmp3

               if (icomp==1) then
                   vr_z1 = vr_z1 + vr_tmp1 
                   vr_z2 = vr_z2 + vr_tmp2 
               elseif (icomp==2) then
                   vr_r1 = vr_r1 + vr_tmp1 
                   vr_r2 = vr_r2 + vr_tmp2 
               elseif (icomp==3) then
                   vr_t1 = vr_t1 + vr_tmp1 
                   vr_t2 = vr_t2 + vr_tmp2 
               endif

            enddo
         enddo
      enddo

      ! old VR calculation: incorrect 
      !vr_tmp5 = dble(nwave)*3d0*dble(rcvnum)
      !VR_array(count_ns) = vr_tmp4 / vr_tmp5
      !write(*,*)  'VR old ',VR_array(count_ns)

      ! New VR calcutaion
      VR_array(count_ns) = (1d0 - (vr_tmp1_sum / vr_tmp2_sum)) * 100d0

      VR_array_cmp(count_ns,1) = (1d0-(vr_z1/vr_z2)) * 100d0 
      VR_array_cmp(count_ns,2) = (1d0-(vr_r1/vr_r2)) * 100d0 
      VR_array_cmp(count_ns,3) = (1d0-(vr_t1/vr_t2)) * 100d0 
      vr_z1 = VR_array_cmp(count_ns,1)
      vr_r1 = VR_array_cmp(count_ns,2)
      vr_t1 = VR_array_cmp(count_ns,3)
      write(*,*)  'VR ',VR_array(count_ns),'Z',vr_z1,'R',vr_r1,'T',vr_t1

      do k = 1, rcvnum

         if (geomet_wei == 0) then
            correct_dist = 1d0
         elseif (geomet_wei == 1) then
            correct_dist = sqrt(distance(k))
         else
            write(*,*) 'wrong in correctdist'
            stop
         endif

         do icomp = 1, 3
            if (icomp==1) then
               wei_cmp = Z_wei
            elseif (icomp==2) then
               wei_cmp = R_wei
            elseif (icomp==3) then
               wei_cmp = T_wei
            endif
            do iw = 1, nwave
              npts = n_data(k,iw)
              if (ZRT_ratio == 1) then
                 max_Z = maxval(observed_data(1:npts,1,iw,k))
                 max_R = maxval(observed_data(1:npts,2,iw,k))
                 max_T = maxval(observed_data(1:npts,3,iw,k))
                 max_now = maxval(observed_data(1:npts,icomp,iw,k))
                 max_cmp = max(max_Z,max_R,max_T)
                 wei_cmp = max_cmp/max_now
                 !write(*,*) icomp,iw,max_now,max_cmp,wei_cmp
              endif
              weightval = real(weight_data(iw))
              wei = weightval / weightsum

              ! devide by npts is not valid maybe. See evernote
              rmis = rmis + sum((Res_data(1:npts,icomp,iw,k)
     &                * wei * wei_cmp * correct_dist) ** rnrm)


            enddo
         enddo
      enddo

      rmis = rmis ** (1/rnrm)
      rmisfitALL = rmis
      rmisfitWAV = rmis
      write(*,*) 'misfit ',rmisfitALL 

      ! Variance Reductions info for the calc test
      open(102,file='./Waveforms/VRs_forward.txt')
      do k = 1, rcvnum
          do icomp = 1, 3
              do iw = 1, nwave
                  vr_tmp1 = VRs_array(iw,icomp,k,1) !invalid VR
                  vr_tmp2 = VR_array(1) !Valid VR
                  misfit_tmp = miss_array(iw,icomp,k,1)
                  write(102,*) k,icomp,iw,vr_tmp1,vr_tmp2,
     &                            misfit_tmp
               enddo
          enddo
      enddo
      close(102)


      !write(*,*) 'stop @misfit_func'
      !stop 

      endsubroutine
               

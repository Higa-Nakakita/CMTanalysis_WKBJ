
      subroutine fixterm_calc(icomp,k)
      ! Compute the fix-term (ray-path and receiver) for each station and each mode.
      ! For each modebranch, perform the following steps:
      ! 1st step:
      ! Make array of eigen‐parameters (frequency, Q, group velocity) as
      ! functions of the angular-order index "ntmp".
      ! 2nd step:
      ! Compute the complex spectra of fix-term in the frequency
      ! domain, and store it in the array "zAexp_fix".

      implicit none
      include "swfi_param.inc"
      include "wkbj_common.inc"
      include "swfi_common.inc"
      include "eigen_common.inc"
      include "new_eigen_common.inc"

      integer STmode, icomp, k
      integer n, j, l, ntmp, model_path, model_rcv

      integer narray(0:nmx,2),l1,l2,jjtmp
      real(8) omgrcv, frqtmp, omgtmp

      ! fixterm calc
      complex(16) calc_rcvterm, calc_pathterm
      complex(16) zAexp_rcv, zAexp_path

      ! calc rcv term
      real(8) omg_rcv(lmx)
      real(8) eif_rcvU(lmx), eif_rcvV(lmx), eif_rcvW(lmx)
      common /calc_rcv/omg_rcv, eif_rcvU, eif_rcvV, eif_rcvW

      ! calc path term
      real(8) dl(lmx), omg_path(lmx)
      real(8) Q_path(lmx), gv_path(lmx)
      common /calc_path/dl,omg_path,Q_path,gv_path
      real(8) del_tmp,sin_D,geosp_tmp
      common /calc_path2/del_tmp,sin_D,geosp_tmp

      !~~~~ local modes  ~~~~~~~~~~~~~ 
      integer div_idx, divnum
      integer Flag_invalid
      real(8) omg_local(maxdivs,lmx)
      real(8) Q_local(maxdivs,lmx),gv_local(maxdivs,lmx)
      common /calc_path_local/omg_local,Q_local,gv_local
      real(8) omg_ck, omg_prem, omg_sup, omg_inf, gv_ck, Q_ck 
      real(8) Q_prem,gv_prem


      if (icomp==1 .or. icomp==2) then ! Vertical or Radial
          STmode = 1  ! Speroidal mode
      elseif (icomp==3) then  ! Transverse
          Stmode = 2  ! Toroidal mode
      endif

      !write(*,*) 'STmode,icomp,k'
      !write(*,*) 'fixterm_calc',STmode,icomp,k

      model_path = 1
      model_rcv  = 2

! -------- constanr variables (of fixterm) ---------
      del_tmp= del_array(k) * d2r

      sin_D  = dsin(del_tmp)
      geosp_tmp = 8d0*pi*dabs(sin_D)
!      write(*,*) 'del_tmp,sin_D,geosp_tmp'
!      write(*,*) del_tmp,sin_D,geosp_tmp
! ------------------------------------


      !open(2,file='ck_zAexprcv.txt',position='append')
      !open(101,file='ck_local_eigp.txt',position='append')
      do n = nn1, nn2
         narray(n,STmode) = 0
         !l1 = lbeg(n,STmode,k)
         !l2 = lend(n,STmode,k)
         ! 1st step
         ! consider local l range
         l1 = lbeg_path(n,STmode,k)
         l2 = lend_path(n,STmode,k)
         do l = l1, l2
            narray(n,STmode) = narray(n,STmode) + 1
            ntmp = narray(n,STmode)
!            write(*,*) n,l,ntmp
            dl(ntmp) = dble(l)
            omg_rcv(ntmp) = dble(eig_para(model_rcv,STmode,n,l,1,k))
            omg_path(ntmp)= dble(eig_para(model_path,STmode,n,l,1,k))
            Q_path(ntmp)  = dble(eig_para(model_path,STmode,n,l,2,k))
            gv_path(ntmp) = dble(eig_para(model_path,STmode,n,l,3,k))
!            write(*,*) omg_path(ntmp),Q_path(ntmp),gv_path(ntmp)
          !  if (n<3) then
          !     if (icomp==1 .or. icomp==3) then
          !        write(101,*) 'PREM',STmode,n,l,k,real(omg_path(ntmp)),
          !&                      real(Q_path(ntmp)),real(gv_path(ntmp))
          !     endif
          !  endif

            omg_prem = omg_path(ntmp)
            omg_sup = omg_prem + omg_prem * 0.8d0 
            omg_inf = omg_prem - omg_prem * 0.8d0
            Q_prem = Q_path(ntmp)
            gv_prem = gv_path(ntmp)
            ! check omg_rcv
            if (omg_rcv(ntmp)<omg_inf .or. omg_rcv(ntmp)>omg_sup) then
                write(*,*) 'Wrong in rcv omg?'
                write(*,*) 'ST,n,l,k',STmode,n,l,k
                write(*,*) 'PREM path omg',omg_prem
                write(*,*) 'inf, sup',omg_inf,omg_sup
                write(*,*) 'rcv omg',omg_rcv(ntmp)
                stop
            endif

            !~~~~~~~~~~ local modes ~~~~~~~~~~~~~~~~~~~
            ! store local eigen parameters for path integral
            divnum = div_num_array(k)
            do div_idx = 1, divnum
               omg_local(div_idx,ntmp) = 
     &      dble(eig_para_local(div_idx,STmode,n,l,1,k))
               Q_local(div_idx,ntmp) = 
     &      dble(eig_para_local(div_idx,STmode,n,l,2,k))
               gv_local(div_idx,ntmp) = 
     &      dble(eig_para_local(div_idx,STmode,n,l,3,k))

               ! Detection of invalid eigpara
               Flag_invalid = 0
               omg_ck = omg_local(div_idx,ntmp)
               Q_ck = Q_local(div_idx,ntmp)
               gv_ck = gv_local(div_idx,ntmp)
               if (omg_ck<omg_inf .or. omg_ck>omg_sup) then
                     write(*,*) 'Wrong in local omg ??'
                     write(*,*) 'ST,n,l,k,idx',STmode,n,l,k,div_idx
                     write(*,*) 'PREM omg',omg_prem
                     write(*,*) 'inf, sup',omg_inf,omg_sup
                     write(*,*) omg_ck
                     Flag_invalid = 1
                     !stop
               endif
               if (Q_ck<0.01d0 .or. Q_ck>2000d0) then
                     write(*,*) 'Wrong in local Q ??'
                     write(*,*) 'ST,n,l,k,idx',STmode,n,l,k,div_idx
                     write(*,*) Q_ck
                     Flag_invalid = 2
                     !stop
               endif
               if (gv_ck<0.01d0 .or. gv_ck>20d0) then
                     write(*,*) 'Wrong in local gv ??'
                     write(*,*) 'ST,n,l,k,idx',STmode,n,l,k,div_idx
                     write(*,*) gv_ck
                     Flag_invalid = 3
                     !stop
               endif

               ! 251118
               ! Remedy invalid eigen parameter 
               ! Replace by PREM value
               if (Flag_invalid/=0) then
                     open(1,file='remedy_eigps.txt',position='append')
                     write(1,*) STmode,n,l,k,div_idx
                     write(1,*) omg_ck,Q_ck,gv_ck
                     write(1,*) omg_prem,Q_prem,gv_prem
                     omg_local(div_idx,ntmp) = omg_prem 
                     Q_local(div_idx,ntmp) = Q_prem 
                     gv_local(div_idx,ntmp) = gv_prem 
                     !write(*,*) omg_prem,Q_prem,gv_prem
                     !write(*,*) omg_ck,Q_ck,gv_ck
                     !stop
                     close(1)
               endif

          !     if (n<3) then
          !        if (icomp==1 .or. icomp==3) then
          !           write(101,*) div_idx,STmode,n,l,k,real(omg_ck),
          !&                      real(Q_ck),real(gv_ck)
          !        endif
          !     endif

            enddo
            !stop
            ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            ! arrange freq_min,max correspond to l for each n
            if(l == lbeg_path(n,STmode,k)) then
               freq_min(n,STmode,k) = omg_path(ntmp) / pi2
!               write(*,*) 'freq_min',freq_min(n,STmode,k)
            elseif (l == lend_path(n,STmode,k)) then
               freq_max(n,STmode,k) = omg_path(ntmp) / pi2
!               write(*,*) 'freq_mx',freq_max(n,STmode,k) 
            endif

            ! rearrange eigen functions of receiver 
            omgrcv = omg_rcv(ntmp)
            if (STmode==1) then ! Spheroidal mode
                eif_rcvU(ntmp) = dble(eif_rcvS(n,l,1,k)) * omgrcv ! U
                eif_rcvV(ntmp) = dble(eif_rcvS(n,l,3,k)) * omgrcv ! V
                !write(*,*) eif_rcvS(n,l,1,k),eif_rcvS(n,l,3,k)
            elseif (STmode==2) then ! Toroidal mode
                eif_rcvW(ntmp) = dble(eif_rcvT(n,l,1,k)) * omgrcv ! W
            endif
            !write(*,*) 'eif_rcvU',ntmp,eif_rcvU(ntmp)

        enddo  ! end loop od angularorder l
        !write(*,*) 'mode=',n,' l=',l1,'~',l2,' last_ntmp=',ntmp 

        ! 2nd step
        jjtmp = 0
        do j = 1, nyquist(k)
          if (n==nn1) then
             frqprt(j) = dfrq(k) * dble(j-1)
             omgprt(j) = pi2 * frqprt(j)
          endif
          frqtmp = frqprt(j)
          omgtmp = omgprt(j)
! ------ set zero for (freq < freqmin(n) or freqmax(n)) ------
          if (frqtmp<freq_min(n,STmode,k)) then
             j_skip1(n,STmode,k) = j
             zAexp_fix(n,j,STmode,icomp,k) = dcmplx(0d0,0d0)
          elseif (frqtmp>freq_max(n,STmode,k)) then
             zAexp_fix(n,j,STmode,icomp,k) = dcmplx(0d0,0d0)
          else
              j_skip2(n,STmode,k) = j
              jjtmp = jjtmp + 1
              j_frq(n,j,STmode) = jjtmp
              ! calc rcv term
              zAexp_rcv = calc_rcvterm(ntmp,STmode,icomp,omgtmp)
              ! calc path term
              zAexp_path= calc_pathterm(n,ntmp,jjtmp,STmode,omgtmp,k)
              !write(*,*)  k,n,j,STmode,icomp,zAexp_path
              ! make fix term of WKBJ
              zAexp_fix(n,j,STmode,icomp,k) = zAexp_rcv * zAexp_path
              !write(*,*) k,n,j,STmode,icomp
              !write(*,*) zAexp_fix(n,j,STmode,icomp,k)
              !write(*,*) zAexp_rcv,zAexp_path
              !write(*,*) zAexp_fix(n,j,STmode,icomp,k)

            !  write(2,*) n,real(frqtmp),STmode,icomp,k,
        !&        real(real(zAexp_rcv)),real(aimag(zAexp_rcv))

          endif
        enddo ! enddo for freq
      enddo ! enddo for modebranch
      !close(2)
      !close(101) !ck_local_eigp.txt


      ! check for eigen func of rcv
      !if (k==1) then
      !   open(1,file='ck_eifrcv_nl.txt',position='append')
      !   do n = nn1, nn2
      !     do l = lbeg(n,STmode,k), lend(n,STmode,k)
      !       write(1,*) n,l,real(eif_rcvS(n,l,1,k)),
      !&         real(eif_rcvS(n,l,3,k)),real(eif_rcvT(n,l,1,k))
      !     enddo
      !   enddo
      !   close(1)
      !endif


      ! writeout wkbj_fixterm
      !open(1,file='ck_zAexpfix.txt',position='append')
      !do n = nn1, nn2
      !   do j = j_skip1(n,STmode,k), j_skip2(n,Stmode,k)
      !     write(1,*) n,j,frqprt(j),STmode,icomp,k,
      !&     real(real(zAexp_fix(n,j,STmode,icomp,k))),
      !&     real(aimag(zAexp_fix(n,j,STmode,icomp,k)))
      !   enddo
      !enddo
      !close(1)

      endsubroutine fixterm_calc



      complex(16) function calc_rcvterm(ntmp,STmode,comp,omgtmp)
      implicit none
      include 'swfi_param.inc'
      ! calc rcv term
      real(8) omg_rcv(lmx)
      real(8) eif_rcvU(lmx), eif_rcvV(lmx), eif_rcvW(lmx)
      common /calc_rcv/omg_rcv, eif_rcvU, eif_rcvV, eif_rcvW
      ! ------------------
      integer ntmp, STmode, comp
      real(8) omgtmp, spline
      real(8) Ur, Vr, Wr
      complex(16) zAexp_rcv

      !write(*,*) 'ntmp'
      !write(*,*) ntmp
      !write(*,*) 'omg_rcv'
      !write(*,*) omg_rcv
      !write(*,*) 'eif_rcvU'
      !write(*,*) eif_rcvU
      !write(*,*) 'omgtmp'
      !write(*,*) omgtmp
      if (STmode==1) then ! Spheroidal
         if(comp==1) then  ! Vertical
            Ur = spline(ntmp,omg_rcv,eif_rcvU,omgtmp)
            zAexp_rcv = dcmplx(Ur,0d0)
         elseif (comp==2) then ! Radial
            Vr = spline(ntmp,omg_rcv,eif_rcvV,omgtmp)
            zAexp_rcv = dcmplx(0d0,-Vr)
         endif
      elseif (STmode==2) then ! Toroidal
         Wr = spline(ntmp,omg_rcv,eif_rcvW,omgtmp)
         zAexp_rcv = dcmplx(0d0,Wr)
      endif

      calc_rcvterm = zAexp_rcv
 
      return
      endfunction calc_rcvterm

      complex(16) function calc_pathterm(n,ntmp,jjtmp,STmode,omgtmp,k)
      ! Perform the line-integration along the ray-path, 
      ! for each ray-path, mode, frequency 
      ! (see textbook Dahlen&Tromp (1998) equation 16.166).
      ! Before the integration, Akima-interportaion will be performed
      ! to make the integrand the function of epicentaral distance.
      implicit none
      include 'swfi_param.inc'
      include 'wkbj_common.inc'
      include 'swfi_common.inc'
      include 'new_eigen_common.inc'
      ! calc path term
      real(8) dl(lmx), omg_path(lmx)
      real(8) Q_path(lmx), gv_path(lmx)
      common /calc_path/dl,omg_path,Q_path,gv_path
      real(8) del_tmp,sin_D,geosp_tmp
      common /calc_path2/del_tmp,sin_D,geosp_tmp
      ! ---------------
      integer n, j, k, ntmp, jjtmp, STmode
      real(8) omgtmp, del
      real(8) spline,dltmp,dkptmp,Qp,gvtmp
      real(8) atten,atten_tmp,geosp,phstmp,denrm_wkbj
      complex(16) zphs, zAexp_rcv

      !~~~~ local modes  ~~~~~~~~~~~~~ 
      integer divnum, div_idx
      real(8) omg_local(maxdivs,lmx)
      real(8) Q_local(maxdivs,lmx),gv_local(maxdivs,lmx)
      common /calc_path_local/omg_local,Q_local,gv_local
      real(8) f,df,ddf, omg_1D(lmx),Q_1D(lmx),gv_1D(lmx)
      real(8) ktmp_local,ltmp_local,Qtmp_local,gvtmp_local 
      real(8) k_del(maxdivs),Q_del(maxdivs),gv_del(maxdivs)
      real(8) atten_kernel(maxdivs)
      real(8) gcpdel, gcpdels(maxdivs)
      real(8) integra, k_integ, atte_integ, atte_integ_tmp
      complex(16) zphs_local, zAexp_path

      ! 2025 1130 
      real(8) gcpdels_shift(maxdivs)
      real(8) k_del_shift(maxdivs),atten_kernel_shift(maxdivs)


      ! ===== 1D layered path calculation ========
      ! obtain angular order
      !dl_path(n,jjtmp,STmode) 
      !& = spline(ntmp,omg_path,dl,omgtmp)
      !dltmp = dl_path(n,jjtmp,STmode)
      dltmp = spline(ntmp,omg_path,dl,omgtmp)

      ! calc wavenumber 
      !dkp(n,jjtmp,STmode) = dsqrt(dltmp*(dltmp+1d0))
      !dkptmp = dkp(n,jjtmp,STmode)
      dkptmp = dsqrt(dltmp*(dltmp+1d0))

      ! obatin Quality factor
      !Qp_array(n,jjtmp,STmode) 
      !& = spline(ntmp,omg_path,Q_path,omgtmp)
      !Qp = Qp_array(n,jjtmp,STmode) 
      Qp = spline(ntmp,omg_path,Q_path,omgtmp)

      ! obtain group vel
      !gv_array(n,jjtmp,STmode)
      !& = spline(ntmp,omg_path,gv_path,omgtmp)
      !gvtmp = gv_array(n,jjtmp,STmode) * dkm2d * d2r
      gvtmp = spline(ntmp,omg_path,gv_path,omgtmp)
      gvtmp = gvtmp  * dkm2d * d2r

      ! ===== local mode path calculation ========
      ! 1st step: make 1D function of eigen param and gcp 
      ! for each (n,omg).
      ! i.e. create integral kernel for path intergal calculation
      divnum = div_num_array(k)
      do div_idx = 1, divnum
         omg_1D(1:ntmp) = omg_local(div_idx,1:ntmp)
         call makima1p_modi(f,df,ddf,omgtmp,ntmp,dl,omg_1D) 
         ltmp_local = f 
         ktmp_local = dsqrt(ltmp_local*(ltmp_local+1d0))
         k_del(div_idx) = ktmp_local
         !write(*,*) 'l',ltmp_local - dltmp
         !write(*,*) 'k',k_del(div_idx), dkptmp
         !write(*,*) k_del(div_idx) - dkptmp

         Q_1D(1:ntmp) = Q_local(div_idx,1:ntmp)
         call makima1p_modi(f,df,ddf,omgtmp,ntmp,Q_1D,omg_1D) 
         Qtmp_local = f 
         Q_del(div_idx) = Qtmp_local
         !write(*,*) 'Q',Qtmp_local - Qp

         gv_1D(1:ntmp) = gv_local(div_idx,1:ntmp)
         call makima1p_modi(f,df,ddf,omgtmp,ntmp,gv_1D,omg_1D) 
         gvtmp_local = f 
         gvtmp_local = gvtmp_local  * dkm2d * d2r
         gv_del(div_idx) = gvtmp_local
         !write(*,*) 'gv',gvtmp_local - gvtmp

         atten_kernel(div_idx) = omgtmp/(2d0*gvtmp_local*Qtmp_local) 

         if (div_idx==1 .and. k==1) then
            !  ==== 20250418 =======
            ! eigpara of path only used for parturbation of source
            ! locaation,
            ! so array must be made for near-source region
            dkp(n,jjtmp,STmode) = ktmp_local
            Qp_array(n,jjtmp,STmode) = Qtmp_local 
            gv_array(n,jjtmp,STmode) = gvtmp_local

            j_frq_path_max(n,STmode) = jjtmp
         endif
      enddo

      ! for visualize image of local eigs along path
      ! do not use 101. 101 is already used.
      !if (n < 3) then
      !   open(102,file='ck_local_integ.txt',position='append')
      !   do div_idx = 1, divnum
      !      write(102,*) k,div_idx,real(gcp_dels(div_idx,k)),
      !&       STmode,n,
      !&       real(omgtmp),real(k_del(div_idx)),
      !&       real(Q_del(div_idx)),real(gv_del(div_idx)/(dkm2d*d2r)),
      !&       real(dkptmp), real(Qp), real(gvtmp/(dkm2d*d2r))
      !   enddo
      !   close(102)
      !endif
         

      ! It seems that unit of del_tmp, gcp_dels is radian
      !write(*,*) gcp_dels(divnum,k), del_tmp
      !write(*,*) 'gcp ',gcp_dels(divnum,k) - del_tmp

      gcpdel = gcp_dels(divnum,k)
      gcpdels(1:divnum) = gcp_dels(1:divnum,k)

      ! ----------------------------------
      ! Remedy 20251129
      ! Create shift array.
      ! It is necessary to pass the x,y array corresponding to 
      ! an epicentral distance 0 to the integration subrtoutine "integra".
      ! Add same value to integrand array for the distance=0 as idx=1.
      ! Add 0 to the gcpdels array.
      do div_idx = 1, divnum
         gcpdels_shift(div_idx+1) = gcpdels(div_idx)
         k_del_shift(div_idx+1) = k_del(div_idx)
         atten_kernel_shift(div_idx+1)= atten_kernel(div_idx)
      enddo
      gcpdels_shift(1) = 0d0
      k_del_shift(1) = k_del(1)
      atten_kernel_shift(1)=atten_kernel(1)
      !write(*,*) "gcpdels_shift: ",gcpdels_shift(1:divnum+1)
      !write(*,*) "k_del_shift: ",k_del_shift(1:divnum+1)
      !write(*,*) "atten_kernel_shift: ",atten_kernel_shift(1:divnum+1)
      ! ----------------------------------


      ! 'Tezuka method'
      ! 'define N_integ, not dx  in integrand subroutine'
      ! 'integ_drad will be ignored'
      !k_integ = integra(divnum,0d0,gcpdel,
      !&  gcpdels(1:divnum),k_del(1:divnum))
      k_integ = integra(divnum+1,0d0,gcpdel,
     &  gcpdels_shift(1:divnum+1),k_del_shift(1:divnum+1))
      zphs_local = dcmplx(dcos(k_integ),-dsin(k_integ))
      !write(*,*) k,n,jjtmp,STmode
      !write(*,*) zphs_local,k_integ
      !write(*,*) 'integ',k_integ, phstmp

      !atte_integ_tmp = integra(divnum,0d0,gcpdel,
      !&  gcpdels(1:divnum),atten_kernel(1:divnum))
      atte_integ_tmp = integra(divnum+1,0d0,gcpdel,
     &  gcpdels_shift(1:divnum+1),atten_kernel_shift(1:divnum+1))
      atte_integ = dexp(-atte_integ_tmp)
      !write(*,*) 'gcpdel',gcpdel
      !write(*,*) 'atte_integ_tmp',atte_integ_tmp
      !if (k==11) then
      !       stop
      !endif


      ! calc atten, phase 
      atten_tmp = omgtmp/(2d0*gvtmp*Qp) * del_tmp
      atten = dexp(-atten_tmp)
      geosp = 1d0 / dsqrt(dkptmp * geosp_tmp)
      phstmp = dkptmp * del_tmp
      zphs = dcmplx(dcos(phstmp),-dsin(phstmp))

      !write(*,*) k_integ, phstmp
      !write(*,*) 'diff k',k_integ - phstmp
      !write(*,*) atte_integ_tmp, atten_tmp
      !write(*,*) 'diff attetmp',atte_integ_tmp - atten_tmp
      !write(*,*) atte_integ, atten
      !write(*,*) 'diff atte',atte_integ - atten
      !write(*,*) 'n,omgtmp', n, omgtmp

      
      ! Maybe this correction is not included in pathterm..
!          * normalisation correction of WKBJ (Dahlen & Tromp, 1998)
      denrm_wkbj = dkptmp / (omgtmp*gvtmp)

      ! ======= check ===========
      !write(*,*) 'path '
      !write(*,*) n,jjtmp
      !write(*,*)  geosp, atten, zphs
      !write(*,*)  geosp*atten*zphs
      !write(*,*)  geosp, atte_integ, zphs_local
      !write(*,*)  geosp*atte_integ*zphs_local


      ! summary
      zAexp_rcv = geosp * atten * zphs
      !calc_pathterm = zAexp_rcv * denrm_wkbj

      zAexp_path = geosp * atte_integ * zphs_local
      calc_pathterm = zAexp_path * denrm_wkbj
      !write(*,*) geosp,atte_integ,zphs_local,zAexp_path, denrm_wkbj

      !open(2,file='ck_pathparam.txt',position='append')
      !   write(2,*) n,jjtmp,real(omgtmp),STmode,k,
      !&              real(Qp),real(gv_array(n,jjtmp,STmode)),real(atten),
      !&              real(geosp),real(phstmp),real(zphs)
      !close(2)



      return
      endfunction calc_pathterm
        





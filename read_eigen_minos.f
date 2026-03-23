
      subroutine rd_eigen_minos(iicomp,k)
      implicit none
      integer k, iicomp
      integer isph1, isph2, itrd, STmode 
      integer divnum, div_idx
      include "swfi_param.inc"
      include "wkbj_common.inc"
      include "swfi_common.inc"
      include "eigen_common.inc"
      include "new_eigen_common.inc"

      write(*,*) 'read_eigen(',iicomp,k,')'

      ! Initialize
      isph1 = 0
      isph2 = 0
      itrd  = 0
      STmode= 0

      ! check STmode type
      if (iicomp==1) then ! Z
         isph1 = 1
      elseif (iicomp==2) then ! R
         isph2 = 1
      elseif (iicomp==3) then ! T
         itrd = 1
      else
         write(*,*) 'Error in iicomp @read_eigen'
         stop
      endif

      ! set component
      if ((isph1==1 .or. isph2==1) .and. itrd==0) then
          STmode = 1 ! Rayleigh
      elseif ((isph1==0 .and. isph2==0) .and. itrd==1) then
          STmode = 2 ! Love
      else
          write(*,*) 'Error in isph @wkbj_reference'
          stop
      endif

      ! ===== Spheroidal mode ========
      ! read eigen parameters (and eigen function?)
      if (STmode==1) then
          write(*,*) 'STmode=',STmode
          call read_eigen(1,'S',k) ! for pathav 
          call read_eigen(2,'S',k) ! for src
          call read_eigen(3,'S',k) ! for rcv
          call determ_lrange('S',k)
      endif

      ! ===== Toroidal mode =========
      if (STmode==2) then
          write(*,*) 'STmode=',STmode
          call read_eigen(1,'T',k)
          call read_eigen(2,'T',k)
          call read_eigen(3,'T',k)
          call determ_lrange('T',k)
      endif



      ! ~~~~~ local mode ~~~~~~~~~~~~~~~
      write(*,*) 
      if (STmode==1) then
          divnum = div_num_array(k)
          write(*,*) divnum
          do div_idx = 1, divnum
            call read_eigen_local(div_idx,'S',k)
          enddo
          call determ_lrange_path('S',k)
      endif
      if (STmode==2) then
          divnum = div_num_array(k)
          do div_idx = 1, divnum
            call read_eigen_local(div_idx,'T',k)
          enddo
          call determ_lrange_path('T',k)
      endif

      !stop
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~

      endsubroutine rd_eigen_minos



! ===============================================================
      subroutine read_eigen(modeltype,cm,kk)
      implicit none
      include "swfi_param.inc"
      include "wkbj_common.inc"
      include "swfi_common.inc"
      include "eigen_common.inc"
      include "new_eigen_common.inc"

      character(len=1) cmd, cm, cha
      integer modeltype, kk
      integer lofw, lw0,lw1,lw2,lw3,lw4,STmode
      integer j, k, l, n, kmax, eigp, eigf, n0, l0
      integer l_flag, nord1, nord2, ll1, ll2
      integer lminS(0:nmx,3,maxstas),lmaxS(0:nmx,3,maxstas)
      integer lminT(0:nmx,3,maxstas),lmaxT(0:nmx,3,maxstas)
      common /rd_eigen/lminS,lmaxS,lminT,lmaxT
      character(len=150) cha_tmp
      character(len=300) eip_filtmp, eif_path_filtmp
      character(len=300) eip_fil, eif_path_fil
      real(8) Tprd


      write(cmd,'(a)') cm
      write(*,*) 'modeltype',modeltype 

      if (modeltype==1) then
         eip_filtmp = trim(prem_path)//'/eigen_pathav_'//cmd
      elseif (modeltype==2) then
         eip_filtmp = trim(prem_path)//'/eigen_pathav_'//cmd
      elseif (modeltype==3) then
         eip_filtmp = trim(prem_path)//'/eigen_pathav_'//cmd
      endif

      ! 251216
      ! for reading eif_rcvS (modified polarity)
      if (cmd=="S") then
          !eip_filtmp = trim(prem_path)//'/eigen_pathav_modi_'//cmd
          eip_filtmp = trim(prem_path)//'/eigen_pathav_modi2_'//cmd
      endif

      lw1 = lofw(eip_filtmp)
      eip_fil(1:lw1) = eip_filtmp(1:lw1)
      lw2 = lofw(eif_path_filtmp)
      eif_path_fil(1:lw2) = eif_path_filtmp(1:lw2)

      if (cm=='S') then
         kmax = 4
         STmode= 1
      elseif (cm=='T') then
         kmax = 2
         STmode= 2
      endif

      ! === open binary files ==============
      ! (eigp): eigen_pathav (n,l,w)
      ! (eigf): eigen_func_pathav (n,l,U,) 
      ! reading format is correspond to weigen (Prog3)
      ! STmode=1:S 2:T
      eigp = 41
      eigf = 42
!      write(*,*) eip_fil(1:lw1)
      open(eigp,file=eip_fil(1:lw1),status='unknown',form='unformatted')
! ------- mode S ------------------------------------------
      if (cmd=='S') then
         read(eigp) nord1, nord2,
     &    (lminS(n,modeltype,kk),lmaxS(n,modeltype,kk),n=nord1,nord2)
         if (nord2 > nmx) then
            nord2 = nmx
            write(*,*) 'nmx adjustment is applied..'
         endif

         do n = nord1, nord2
            ll1 = lminS(n,modeltype,kk)
            ll2 = lmaxS(n,modeltype,kk) 
            l_flag = 0 ! is used for arrangement of lminS
            if (ll2 > lmx) then
               ll2 = lmx
               write(*,*) 'lmx adjustment is applied..'
            endif

            do l = ll1, ll2
               if (modeltype==1) then ! eigen_pathav
                  read(eigp) n0,l0,
     &           (eig_para(modeltype,STmode,n,l,j,kk),j=1,3),
     &           (eif_rcvS(n,l,k,kk),k=1,kmax)
               elseif (modeltype==2) then ! read eigen_src
                  read(eigp) n0,l0,
     &           (eig_para(modeltype,STmode,n,l,j,kk),j=1,3)
         !&           (eif_srcS(n,l,k,kk),k=1,kmax)

               elseif (modeltype==3) then ! read eigen_rcv
                 read(eigp) n0,l0,
     &           (eig_para(modeltype,STmode,n,l,j,kk),j=1,3)
          !&           (eif_rcvS(n,l,k,kk),k=1,kmax)
               endif

               Tprd = pi2 / eig_para(modeltype,STmode,n,l,1,kk)
               if (Tprd<prdmax .and. l_flag==0) then
                  lminS(n,modeltype,kk) = l
                  l_flag = 1
               endif
               if (Tprd>prdmin) then
                  lmaxS(n,modeltype,kk) = l
               endif
             enddo
         enddo  
! ----- mode S -----------------------------------------

! ----- mode T -----------------------------------------
      elseif(cmd=='T') then
         write(*,*) 'T'
         read(eigp) nord1, nord2,
     &   (lminT(n,modeltype,kk),lmaxT(n,modeltype,kk),n=nord1,nord2)
         if(nord2>nmx) then
            nord2 = nmx
         endif

        do n = nord1, nord2
          ll1 = lminT(n,modeltype,kk)
          ll2 = lmaxT(n,modeltype,kk)
          l_flag = 0
          if (ll2 > lmx) then
             ll2 = lmx
          endif

          do l = ll1, ll2
             if (modeltype==1) then
                 read(eigp) n0,l0,
     &           (eig_para(modeltype,STmode,n,l,j,kk),j=1,3)

             elseif (modeltype==2) then
                 read(eigp) n0,l0,
     &           (eig_para(modeltype,STmode,n,l,j,kk),j=1,3)
           !&           (eif_srcT(n,l,k,kk),k=1,kmax)

             elseif (modeltype==3) then
                 read(eigp) n0,l0,
     &           (eig_para(modeltype,STmode,n,l,j,kk),j=1,3),
     &           (eif_rcvT(n,l,k,kk),k=1,kmax)
             endif

             Tprd = pi2 /eig_para(modeltype,STmode,n,l,1,kk)
             if (Tprd<prdmax .and. l_flag==0) then
                lminT(n,modeltype,kk) = l
                l_flag = 1
             endif
             if (Tprd>prdmin) then
                lmaxT(n,modeltype,kk) = l
             endif
          enddo
        enddo
      endif      
! ------ mode T over --------------------------------------
      close(eigp)

      endsubroutine read_eigen



      subroutine determ_lrange(cm,k)
      implicit none

      include "swfi_param.inc"
      include "wkbj_common.inc"
      include "swfi_common.inc"
      include "eigen_common.inc"

      integer k
      character(len=1) cm
      integer n,STmode,l1,l2
      real(4) omega_max, omega_min, Tmin, Tmax
      integer lminS(0:nmx,3,maxstas),lmaxS(0:nmx,3,maxstas)
      integer lminT(0:nmx,3,maxstas),lmaxT(0:nmx,3,maxstas)
      common /rd_eigen/lminS,lmaxS,lminT,lmaxT

      write(*,*) 'cm,k,rcvnum,nn1,nn2'
      write(*,*) cm,k,rcvnum,nn1,nn2

      ! set the range of angular order l in case of multi mode?
      if (cm=='S') then
         STmode = 1
      elseif (cm=='T') then
         STmode = 2
      endif

      open(1,file='modes_range_info.txt',position='append')
      write(1,*) 'PREM',k,rcvnm(k)
      do n = nn1, nn2
        !write(*,*) 'n=',n
        !write(*,*) 'Smin',lminS(n,1,k),lminS(n,2,k),lminS(n,3,k)
        !write(*,*) 'Smax',lmaxS(n,1,k),lmaxS(n,2,k),lmaxS(n,3,k)
        !write(*,*) 'Tmin',lminT(n,1,k),lminT(n,2,k),lminT(n,3,k)
        !write(*,*) 'Tmax',lmaxT(n,1,k),lmaxT(n,2,k),lmaxT(n,3,k)
        if(STmode==1) then
          lbeg(n,STmode,k) =
     &    max0(lminS(n,1,k),lminS(n,2,k),lminS(n,3,k))
          lend(n,STmode,k) =
     &    min0(lmaxS(n,1,k),lmaxS(n,2,k),lmaxS(n,3,k))
        elseif(STmode==2) then
          lbeg(n,STmode,k) =
     &    max0(lminT(n,1,k),lminT(n,2,k),lminT(n,3,k))
          lend(n,STmode,k) =
     &    min0(lmaxT(n,1,k),lmaxT(n,2,k),lmaxT(n,3,k))
        endif

         l1 = lbeg(n,STmode,k)
         l2 = lend(n,STmode,k)
         write(*,*) l1,l2

         omega_max = eig_para(1,STmode,n,l2,1,k)
         omega_min = eig_para(1,STmode,n,l1,1,k)
         Tmin = pi2 / omega_max
         Tmax = pi2 / omega_min

         write(*,'(A1,A,i4,A,i5,A,i5,A,f8.2,A,f8.2)')
     &   cm,'  n:',n,'  l:',lbeg(n,STmode,k),' -',lend(n,STmode,k),
     &   ';  periods:',Tmin,' -',Tmax

         write(1,'(A1,A,i4,A,i5,A,i5,A,f8.2,A,f8.2)')
     &   cm,'  n:',n,'  l:',lbeg(n,STmode,k),' -',lend(n,STmode,k),
     &   ';  periods:',Tmin,' -',Tmax

       enddo
       close(1)

      endsubroutine determ_lrange


      subroutine rd_eif_src
      ! Read the eigenfunction file that include depth information.
      implicit none
      include "swfi_param.inc"
      include "wkbj_common.inc"
      include "swfi_common.inc"
      include "eigen_common.inc"
      include "new_eigen_common.inc"

      character(len=200) cha_tmp1, cha_tmp2
      integer lofw, lw1, lw2
      integer n, l, m, nnn1, nnn2, llmin(0:nmx), llmax(0:nmx)
  

      !cha_tmp1 = eigendir_loca(1) 
      !lw1 =  lofw(cha_tmp1)
      !cha_tmp2 = cha_tmp1(1:lw1)//'/source_eif.bin'
      !cha_tmp2 = trim(prem_path)//'/source_eif_modi.bin'
      cha_tmp2 = trim(prem_path)//'/source_eif_modi_kai.bin'
      lw2 =  lofw(cha_tmp2)
      

      open(102,file=cha_tmp2(1:lw2),form='unformatted')
      read(102) numlay,(src_depth(m),m=1,numlay)
      read(102) nnn1, nnn2, (llmin(n),llmax(n),n=nnn1,nnn2)
      do n =nn1, nn2
        do l =llmin(n), llmax(n)
          do m = 1, numlay
            read(102) src_eifS_z(n,l,m,1),src_eifS_z(n,l,m,2),
     &      src_eifS_z(n,l,m,3),src_eifS_z(n,l,m,4), 
     &      src_eifT_z(n,l,m,1),src_eifT_z(n,l,m,2)
          enddo
        enddo
      enddo
      close(102)


      write(*,*) 'rd_eif_src; OVER'
      endsubroutine


      subroutine read_eigen_local(div_idx,cm,k)
      implicit none
      include "swfi_param.inc"
      include "wkbj_common.inc"
      include "swfi_common.inc"
      include "new_eigen_common.inc"
      integer i, j, n, n0, l0, l, nord1, nord2, ll1, ll2
      integer div_idx, k, kmax, STmode
      integer U_eigP, l_flag
      real(8) Tprd
      character(len=1) cm, cmd
      character(len=150) cha_tmp
      character(len=200) eip_filtmp, eip_fil

      write(cmd,'(a)') cm
      cha_tmp = cha_div_eig(div_idx,k)
      eip_filtmp = trim(cha_tmp)//'/eigen_pathav_'//cmd
      eip_fil = trim(eip_filtmp)
      !write(*,*) k,div_idx,trim(eip_fil)

      if (cm=='S') then
         kmax = 4
         STmode= 1
      elseif (cm=='T') then
         kmax = 2
         STmode= 2
      endif

      U_eigP = 41
      open(U_eigP,file=eip_fil,status='unknown',form='unformatted')
! ------- mode S begin ------------------------------------------
      if (cmd=='S') then
          ! read(U_eigP) nord1, nord2
          ! write(*,*) nord1, nord2
         read(U_eigP) nord1, nord2,
     &   (lminS_path(div_idx,n,k),lmaxS_path(div_idx,n,k),n=nord1,nord2)
      !       write(*,*) nord1, nord2,
      !&   (lminS_path(div_idx,n,k),lmaxS_path(div_idx,n,k),n=nord1,nord2)
         if (nord2 > nmx) then
            nord2 = nmx
            write(*,*) 'nmx adjustment is applied..'
         endif

         do n = nord1, nord2
            ll1 = lminS_path(div_idx,n,k)
            ll2 = lmaxS_path(div_idx,n,k)
            l_flag = 0 ! is used for arrangement of lminS
            if (ll2 > lmx) then
               ll2 = lmx
               write(*,*) 'lmx adjustment is applied..'
            endif

            do l = ll1, ll2
              read(U_eigp) n0,l0,
     &        (eig_para_local(div_idx,STmode,n,l,j,k),j=1,3)
         !         write(*,*) div_idx, n0,l0,
         !&        (eig_para_local(div_idx,STmode,n,l,j,k),j=1,3)

              Tprd = pi2 / eig_para_local(div_idx,STmode,n,l,1,k)
              if (Tprd<prdmax .and. l_flag==0) then
                 lminS_path(div_idx,n,k) = l
                 l_flag = 1
              endif
              if (Tprd>prdmin) then
                 lmaxS_path(div_idx,n,k) = l
              endif
            enddo
         enddo
! ----- mode S over  -----------------------------------------

! ----- mode T begin -----------------------------------------
      elseif(cmd=='T') then
         read(U_eigp) nord1, nord2,
     &   (lminT_path(div_idx,n,k),lmaxT_path(div_idx,n,k),n=nord1,nord2)
         if(nord2>nmx) then
            nord2 = nmx
         endif

        do n = nord1, nord2
          ll1 = lminT_path(div_idx,n,k)
          ll2 = lmaxT_path(div_idx,n,k)
          l_flag = 0
          if (ll2 > lmx) then
             ll2 = lmx
          endif

          do l = ll1, ll2
            read(U_eigp) n0,l0,
     &      (eig_para_local(div_idx,STmode,n,l,j,k),j=1,3)

            Tprd = pi2 /eig_para_local(div_idx,STmode,n,l,1,k)
            if (Tprd<prdmax .and. l_flag==0) then
               lminT_path(div_idx,n,k) = l
               l_flag = 1
            endif
            if (Tprd>prdmin) then
               lmaxT_path(div_idx,n,k) = l
            endif
          enddo
        enddo

      endif
! ------ mode T over --------------------------------------
      close(U_eigP)
      endsubroutine read_eigen_local


      subroutine determ_lrange_path(cm,k)
      implicit none
      include "swfi_param.inc"
      include "wkbj_common.inc"
      include "swfi_common.inc"
      include "eigen_common.inc" !to share lbeg, lend
      include "new_eigen_common.inc"

      character(len=1) cm
      integer i, j, k, n, l, STmode, l1, l2
      integer lmin_tmp1,lmin_tmp2,lmax_tmp1,lmax_tmp2
      integer div_idx, divnum, div_max_idx, div_min_idx
      real(4) omega_max, omega_min, Tmin, Tmax

      if (cm=='S') then
         STmode = 1
      elseif (cm=='T') then
         STmode = 2
      endif

      divnum = div_num_array(k)
      if (STmode==1) then
         do n = nn1, nn2
            lmin_tmp1 = 0 
            lmax_tmp1 = lmx+1 
            div_min_idx = 0
            div_max_idx = 0
            do div_idx = 1, divnum
               lmin_tmp2 = lminS_path(div_idx,n,k) 
               if (lmin_tmp2 > lmin_tmp1) then
                   lmin_tmp1 = lmin_tmp2
                   div_min_idx = div_idx
               endif
               lmax_tmp2 = lmaxS_path(div_idx,n,k)
               if (lmax_tmp2 < lmax_tmp1) then
                   lmax_tmp1 = lmax_tmp2
                   div_max_idx = div_idx
               endif
            enddo
            lbeg_path(n,STmode,k) = lmin_tmp1
            lend_path(n,STmode,k) = lmax_tmp1
         enddo
      elseif (STmode==2) then
         do n = nn1, nn2
            lmin_tmp1 = 0                 
            lmax_tmp1 = lmx+1
            div_min_idx = 0
            div_max_idx = 0
            do div_idx = 1, divnum
               lmin_tmp2 = lminT_path(div_idx,n,k)
               if (lmin_tmp2 > lmin_tmp1) then
                   !write(*,*) n,'lmin_tmp2, lmin_tmp1'
                   !write(*,*) lmin_tmp2, lmin_tmp1
                   lmin_tmp1 = lmin_tmp2
                   div_min_idx = div_idx
               endif
               lmax_tmp2 = lmaxT_path(div_idx,n,k)
               if (lmax_tmp2 < lmax_tmp1) then
                   !write(*,*) n,'lmax_tmp2, lmax_tmp1'
                   !write(*,*) lmax_tmp2, lmax_tmp1
                   lmax_tmp1 = lmax_tmp2
                   div_max_idx = div_idx
               endif
            enddo
            lbeg_path(n,STmode,k) = lmin_tmp1
            lend_path(n,STmode,k) = lmax_tmp1
         enddo
      else
          write(*,*) 'Wrong in STmode',STmode
          write(*,*) 'Stop @ determ_lrange_path',k
          stop
      endif

      do n = nn1, nn2
        if(STmode==1) then
          lbeg_path(n,STmode,k) =
     &    max0(lbeg_path(n,STmode,k),lbeg(n,STmode,k))
          lend_path(n,STmode,k) =
     &    min0(lend_path(n,STmode,k),lend(n,STmode,k))
        elseif(STmode==2) then
          lbeg_path(n,STmode,k) =
     &    max0(lbeg_path(n,STmode,k),lbeg(n,STmode,k))
          lend_path(n,STmode,k) =
     &    min0(lend_path(n,STmode,k),lend(n,STmode,k))
        endif
      enddo

      !Cut computational cost 1 (Cut the upper part of high freq)
      if (cut_cost1==1) then
        do n = nn1, nn2
           l1 = lbeg_path(n,STmode,k)
           l2 = lend_path(n,STmode,k)
           omega_max = eig_para_local(div_max_idx,STmode,n,l2,1,k)
           Tmin = pi2 / omega_max
           !write(*,*) STmode,k,n,lend_path(n,STmode,k),Tmin
           do l = l1,l2
              omega_max = eig_para_local(div_max_idx,STmode,n,l,1,k)
              Tmin = pi2 / omega_max
              if (Tmin < prd_inf) then
                 write(*,*) n,l,Tmin
                 lend_path(n,STmode,k) = l
                 exit
              endif
           enddo
        enddo
      endif



      write(*,*) 'Local l range'
      open(1,file='modes_range_info.txt',position='append')
      write(1,*) 'Local modes',k, rcvnm(k) 
      do n = nn1, nn2
         l1 = lbeg_path(n,STmode,k)
         l2 = lend_path(n,STmode,k)
         omega_max = eig_para_local(div_max_idx,STmode,n,l2,1,k)
         omega_min = eig_para_local(div_min_idx,STmode,n,l1,1,k)
         Tmin = pi2 / omega_max
         Tmax = pi2 / omega_min
         write(*,'(A1,A,i4,A,i5,A,i5,A,f8.2,A,f8.2)')
     &   cm,'  n:',n,'  l:',
     &   lbeg_path(n,STmode,k),' -',lend_path(n,STmode,k),
     &   ';  periods:',Tmin,' -',Tmax
         write(1,'(A1,A,i4,A,i5,A,i5,A,f8.2,A,f8.2)')
     &   cm,'  n:',n,'  l:',
     &   lbeg_path(n,STmode,k),' -',lend_path(n,STmode,k),
     &   ';  periods:',Tmin,' -',Tmax
      enddo
      write(1,*) 
      close(1)

      endsubroutine determ_lrange_path

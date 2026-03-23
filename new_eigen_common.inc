!    eigen dir loca
      character(len=150) cha_rcv_eif(maxstas)
      character(len=150) cha_src_eif
      character(len=300) cha_div_eig(maxdivs,maxstas)
      common /eif_dirpath/cha_rcv_eif,cha_src_eif,cha_div_eig

!    eigen parameters
      real(4) eig_para_local(maxdivs,2,0:nmx,lmx,3,maxstas)
      real(4) eig_para_rcv_local(2,0:nmx,lmx,3,maxstas)
      real(4) eig_para_src_local(2,0:nmx,lmx,3)
      common /new_eigpara/eig_para_local,eig_para_rcv_local
     &                  , eig_para_src_local

!     eigen functions of receivers
      real(4) eif_rcvS_local(0:nmx,lmx,4,maxstas)
      real(4) eif_rcvT_local(0:nmx,lmx,2,maxstas)
      common /eigfunc_rcv/eif_rcvS_local, eif_rcvT_local

!     eigen functions of source
!      integer numlay
!      real(4) src_depth(1000)
!      real(4) src_eifS_z(0:nmx,lmx,1000,4) 
!      real(4) src_eifT_z(0:nmx,lmx,1000,2)
!      common /source_depth/numlay, src_depth, src_eifS_z, src_eifT_z

!     local mode
      integer div_num_array(200)
      real(8) gcp_dels(maxdivs,maxstas)
      common /localmodes/div_num_array, gcp_dels


!     angular order
      integer lminS_path(maxdivs,0:nmx,maxstas)
      integer lmaxS_path(maxdivs,0:nmx,maxstas)
      integer lminT_path(maxdivs,0:nmx,maxstas)
      integer lmaxT_path(maxdivs,0:nmx,maxstas)
      common /l_path/lminS_path,lmaxS_path,lminT_path,lmaxT_path

!     l range (shape of array not changed actually... pls check)
      integer lbeg_path(0:nmx,2,maxstas), lend_path(0:nmx,2,maxstas)
      common /lrange_path/lbeg_path, lend_path

!    path to eigendir 251216
      character(len=200) prem_path
      common /eigendirpath/prem_path

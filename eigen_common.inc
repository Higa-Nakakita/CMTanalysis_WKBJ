
!    eigen parameters
      real(4) eig_para(3,2,0:nmx,lmx,3,maxstas)
      common /eigpara/eig_para

!    eigen functions
      real(4) eif_rcvS(0:nmx,lmx,4,maxstas)
      real(4) eif_rcvT(0:nmx,lmx,2,maxstas)
      real(4) eif_srcS(0:nmx,lmx,4,maxstas)   !eif_srcS/T is not used.
      real(4) eif_srcT(0:nmx,lmx,2,maxstas)   
      common /eigfunc/eif_rcvS, eif_rcvT, eif_srcS, eif_srcT

!    angular order
      integer lbeg(0:nmx,2,maxstas), lend(0:nmx,2,maxstas)
      common /l_arrange/lbeg, lend

!    eigen functions of source
      integer numlay
      real(4) src_depth(1000)
      real(4) src_eifS_z(0:nmx,lmx,1000,4) 
      real(4) src_eifT_z(0:nmx,lmx,1000,2)
      common /source_depth/numlay, src_depth, src_eifS_z, src_eifT_z


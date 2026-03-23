!    common file for SurfaceWaveForm Inversion

!  source information
      real(8) theSc,phiSc,depth
      common /src_info/theSc,phiSc,depth

!  receivers information
      integer rcvnum, nwave
      common /rcv_window/rcvnum, nwave
      character(len=10) rcvnm(maxstas)
      common /rcv_name/rcvnm
      real(8) rcv_latlon(maxstas,2)
      common /rcvs_coordinate/rcv_latlon
      real(8) az_array(maxstas),baz_array(maxstas)
     &   ,del_array(maxstas),distance(maxstas)
      common /azdel_info/az_array,baz_array,del_array,distance


!   normal_mode info
      integer nn1, nn2
      common /normalmode/nn1,nn2

!   fittingh information
      integer n_data(maxstas,maxwindow),n_fft(maxstas,maxwindow)
      integer geomet_wei, ZRT_ratio
      real(4) t_begin(maxstas,maxwindow),t_end(maxstas,maxwindow), tint
      real(4) t0_buffer
      real(8) rnrm
      real(8) weight_data(maxwindow),bp1(maxwindow),bp2(maxwindow)
      real(8) bpfrq_min(maxwindow), bpfrq_max(maxwindow)
      real(4) Z_wei, R_wei, T_wei
      common /fit_info/n_data,n_fft,t_begin,t_end,tint,t0_buffer
      common /fit_info2/weight_data,rnrm,bp1,bp2,bpfrq_min,bpfrq_max
      common /fit_info3/geomet_wei, ZRT_ratio, Z_wei, R_wei, T_wei

!   discrete freq-domain information (WKBJ calculation)
      integer nyquist(maxstas), ntmax(maxstas)
      real(8) dfrq(maxstas), domega(maxstas)
      common /freq_interval/dfrq,domega,nyquist,ntmax

      real(8) freq_min(0:nmx,2,maxstas), freq_max(0:nmx,2,maxstas)
      common /freq_range/freq_min, freq_max

      real(8) frqprt(maxdata), omgprt(maxdata)
      common /freq_part/frqprt,omgprt

      integer j_skip1(0:nmx,2,maxstas), j_skip2(0:nmx,2,maxstas)
      integer j_frq(0:nmx,lmx,2),j_frq_path_max(0:nmx,2)
      common /jrange/j_skip1,j_skip2,j_frq,j_frq_path_max

      integer j_skip1_adjst(0:nmx,2),j_skip2_adjst(0:nmx,2)
      common /jrange_adjst/j_skip1_adjst,j_skip2_adjst


!   constarin info, configure info for inversion
      integer MT_CMT,trM,conf_rcv,conf_path,conf_fix,conf_src
      integer origT
      common /conf_trM/trM,conf_rcv,conf_path,conf_fix,conf_src
      common /condi_time/origT
      common /inv_flow/MT_CMT

!   Path parameters  (created by calc_pathterm)
      real(8) dl_path(0:nmx,jmx,2), dkp(0:nmx,jmx,2)
      real(8) Qp_array(0:nmx,jmx,2), gv_array(0:nmx,jmx,2)
      common /path_perturb/dl_path, dkp, Qp_array, gv_array

!  srcterm 
      integer iobstype
      common /obstype/iobstype

!   essential data for invrsion
!   freq-domain waveforms
      complex(16) zAexp_fix(0:nmx,maxdata,2,numcomp,maxstas)
      complex(16) zAexp_src(0:nmx,maxdata,2,maxstas)
      complex(16) zAexp_frq(0:nmx,maxdata,2,numcomp,maxstas)
      common /freq_wave/zAexp_fix, zAexp_src, zAexp_frq
     
!   time-domain waveforms
      real(4) observed_data(maxdata,numcomp,maxwindow,maxstas)
      real(4) synthetic_data(maxdata,numcomp,maxwindow,maxstas)
      real(4) Res_data(maxdata,numcomp,maxwindow,maxstas)
      common /waveforms/observed_data,synthetic_data,Res_data


!    iso or aniso (20250420)
      integer iso_ani, src_iso_ani, which_struc
      common /iso_aniso/iso_ani, src_iso_ani, which_struc

!    Cut computational cost
      integer cut_cost1,mode_skip
      real(4) prd_inf
      common /cut_cost/cut_cost1,mode_skip,prd_inf

!   Variance Reduction
      integer count_ns
      real(4) VRs_array(maxwindow,numcomp,maxstas,70000) !nsample_max*(nitr_max+1)
      real(4) miss_array(maxwindow,numcomp,maxstas,70000)
      real(4) VR_array(70000), VR_array_cmp(70000,3)
      common /VR/count_ns,VRs_array,miss_array,VR_array,VR_array_cmp

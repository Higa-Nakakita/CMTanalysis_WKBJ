

      subroutine w_seis(wave,iw,iobssynres,k,comp,nitr)
      implicit none
      include 'swfi_param.inc'
      include 'swfi_common.inc'
      include 'wkbj_common.inc'
      real(4) wave(maxdata)
      integer nitr, iw, iobssynres, k, comp
      character(len=1) cmp, tw
      character(len=50) filnm_tmp,filnm
      real(4) ytb
      integer npts, lofw, lw, lw1, lw2, lw3, nerr
      character(len=50) dirnm_tmp, dirnm, filpath_tmp, filpath
      character(len=5) conf

      ! Output directory
      dirnm_tmp = "./Waveforms/"

      ! set component
      if (comp==1) then
              write(cmp,"(a1)") "Z"
      elseif (comp==2) then
              write(cmp,"(a1)") "R"
      elseif (comp==3) then
              write(cmp,"(a1)") "T"
      else
              write(*,*) "stop @ w_seis"
              stop
      endif

      ! set num of time window
      write(tw,"(i1)") iw


      !convert nitr to character
      if(nitr < 10) then
         write(conf,'("0000",i1)') nitr
      else if(nitr < 100) then
         write(conf,'("000",i2)') nitr
      else if(nitr < 1000) then
         write(conf,'("00",i3)') nitr
      else if(nitr < 10000) then
         write(conf,'("0",i4)') nitr
      else
         write(conf,'(i5)') nitr
      end if
      !write(*,*) conf

      ! set file name
      if (iobssynres==0) then
              filnm_tmp=
     &      "obs_"//trim(rcvnm(k))//"_tw"//tw//".LH"//cmp//".sac"
      elseif (iobssynres==1) then
              filnm_tmp=
     &      conf//".syn_"//trim(rcvnm(k))//"_tw"//tw//".LH"//cmp//".sac"
      elseif (iobssynres==2) then
              filnm_tmp=
     &      conf//".res_"//trim(rcvnm(k))//"_tw"//tw//".LH"//cmp//".sac"
      elseif (iobssynres==3) then ! For forward check
              filnm_tmp=
     &      "syn_"//trim(rcvnm(k))//"_tw"//tw//".LH"//cmp//".sac"
      else
              write(*,*) "stop @ w_seis"
              stop
      endif
      !write(*,*) filnm_tmp


      ! writeout seismogram
      ytb = t_begin(k,iw)
      npts= n_data(k,iw)

      lw1 = lofw(filnm_tmp)
      filnm(1:lw1) = filnm_tmp(1:lw1)
      lw2 = lofw(dirnm_tmp)
      dirnm(1:lw2) = dirnm_tmp(1:lw2)
      filpath_tmp = dirnm(1:lw2)//filnm(1:lw1)
      lw3 = lofw(filpath_tmp)
      filpath(1:lw3) = filpath_tmp(1:lw3)
      !write(*,*) filpath_tmp
      !call writedata(filpath_tmp,npts,wave,ytb,tint,lw3)
      call wsac1(filpath_tmp,wave,npts,ytb,tint,nerr)
      !sac func (2025, 4/2) 

      end

      subroutine w_seis_full(wave,iw,iobssyn,k,comp)
      implicit none
      include 'swfi_param.inc'
      include 'swfi_common.inc'
      include 'wkbj_common.inc'
      real(4) wave(maxdata)
      integer nitr, iw, iobssyn, k, comp
      character(len=1) cmp, tw
      character(len=30) filnm_tmp,filnm
      real(4) ytb
      integer npts, lofw, lw, lw1, lw2, lw3, nerr
      character(len=50) dirnm_tmp, dirnm, filpath_tmp, filpath

      dirnm_tmp = "./full_Waveforms/"
            ! set component
      if (comp==1) then
              write(cmp,"(a1)") "Z"
      elseif (comp==2) then
              write(cmp,"(a1)") "R"
      elseif (comp==3) then
              write(cmp,"(a1)") "T"
      else
              write(*,*) "stop @ w_seis (comp) full"
              stop
      endif

      ! set num of time window
      write(tw,"(i1)") iw

      ! set file name
      if (iobssyn==0) then
          filnm_tmp=
     &    "obs_"//trim(rcvnm(k))//"_tw"//tw//"_full.LH"//cmp//".sac"
      elseif (iobssyn==1) then
          filnm_tmp=
     &    "syn_"//trim(rcvnm(k))//"_tw"//tw//"_full.LH"//cmp//".sac"
      else
              write(*,*) "stop @ w_seis (obssyn) full"
              stop
      endif

      ! writeout seismogram
      ytb = 0
      npts= ntmax(k)

      lw1 = lofw(filnm_tmp)
      filnm(1:lw1) = filnm_tmp(1:lw1)
      lw2 = lofw(dirnm_tmp)
      dirnm(1:lw2) = dirnm_tmp(1:lw2)
      filpath_tmp = dirnm(1:lw2)//filnm(1:lw1)
      lw3 = lofw(filpath_tmp)
      filpath(1:lw3) = filpath_tmp(1:lw3)

      !call writedata(filpath_tmp,npts,wave,ytb,tint,lw3)
      call wsac1(filpath_tmp,wave,npts,ytb,tint,nerr)

      end




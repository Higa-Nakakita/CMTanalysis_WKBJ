
      ! create div_del_info (read by WKBJ code)
      ! calc delta (radian) of gcp cross each grids.
      ! and store reference coordinates.
      ! Necessary files is follows.
      ! output_srclonlat.txt (made by gcp_del.py)
      program calc_divdel

      implicit none
      integer i, j, n, k, maxrow, rcvs_num, num, del_idx
      integer ios
      integer STmode, nord1, nord2, nmx, lmx
      integer, allocatable ::lminS(:),lmaxS(:),lminT(:),lmaxT(:)
      real(8) pi, pi2, r2d, d2r, dkm2d, d2km, dkm2r, m2km
      real(8) lon, lat, slat, slon, sdep
      real(8) cross_lat1,cross_lon1,cross_lat2,cross_lon2,g_lat1,g_lat2,g_lon1,g_lon2,d_del
      real(8) ref_glat, ref_glon, dist_sum
      real(8) geocen_deg 
      real(8) thelat1_geoc,thelat2_geoc,az,baz,del,d_dist
      real(8), allocatable :: rcvs_lat(:), rcvs_lon(:)
      real(8) dist_obspy_sum, diff
      character(len=10), allocatable :: rcvs_name(:)
      character(len=10) rcvname
      character(len=50) fil, target_dir, target_path, path_eip, rcvs_fil
      character(len=50), allocatable :: target_array(:)
      character(len=6) cha_slon
      character(len=5) cha_slat
      character(len=5) cha_lon
      character(len=4) cha_lat
      character(len=1) cmd

      maxrow = 100000
      nmx = 15
      lmx = 560

      call ref_inc

      ! ==== compute constant number (wkbj_common.inc) ======
      pi  = 4d0 * datan(1d0)
      pi2 = 2d0 * pi
      r2d = 180d0 / pi
      d2r = pi / 180d0
      dkm2d = r2d / 6371d0 ! convert km  -> deg on the earth 
      d2km  = d2r * 6371d0 ! convert deg -> km  on the earth
      dkm2r = 1d0 / 6371d0 ! convert km -> rad  on the earth
      m2km = 1e-3

      open(1,file='div_del_info')
      close(1,status='delete')

      open(1,file='src_coordinates',action='read',status='old')
      read(1,*) slat,slon,sdep
      write(cha_slon,'(f6.2)') slon
      write(cha_slat,'(f5.2)') slat
      close(1)


      !rcvs_fil = 'rcv_coordinates'
      rcvs_fil = '/work94A/higa/Fnet_data/Fnet_stalist_ref'
      open(1,file=rcvs_fil,action='read',status='old')
      k = 0
      do i = 1, maxrow
         read(1,*,iostat=ios)
         if (ios/=0) then
            exit
         endif
         k = k + 1
      enddo
      close(1)
      rcvs_num = k

      allocate(rcvs_lat(rcvs_num))
      allocate(rcvs_lon(rcvs_num))
      allocate(rcvs_name(rcvs_num))
      open(1,file=rcvs_fil,action='read',status='old')
      do k = 1, rcvs_num
         read(1,*) rcvname, lat, lon
         write(*,*) rcvname, lat, lon
         rcvs_name(k) = trim(rcvname)
         rcvs_lat(k) = lat
         rcvs_lon(k) = lon
      enddo
      close(1)


      ! store d_del info and passed grids info into array
      fil = 'output_'//cha_slat//'_'//cha_slon//'.txt'
      write(*,*) fil
      open(1,file=fil,action='read',status='old')
      k = 0
      do i = 1, maxrow
         write(*,*) i
         read(1,*,iostat=ios)
         if (ios/=0) then
            exit
         endif
         k = k + 1
      enddo
      close(1)
      num = k
      write(*,*) num

      del_idx = 0
      dist_sum = 0d0
      dist_obspy_sum = 0d0  ! to detect error
      open(2,file='div_del_info')
      open(1,file=fil,action='read',status='old')
      do i = 1, num
         read(1,*) cross_lat1,cross_lon1,cross_lat2,cross_lon2,g_lat1,g_lat2,g_lon1,g_lon2,d_del
         !write(*,*) cross_lat1,cross_lon1,cross_lat2,cross_lon2,g_lat1,g_lat2,g_lon1,g_lon2,d_del
         if (cross_lat1==cross_lat2 .and. cross_lon1==cross_lon2) then
              cycle
         endif
         dist_obspy_sum = dist_obspy_sum + d_del 

         thelat1_geoc = geocen_deg(cross_lat1)
         thelat2_geoc = geocen_deg(cross_lat2)
         call azdel(thelat1_geoc,cross_lon1,thelat2_geoc,cross_lon2,az,baz,del)
         !write(*,*) thelat1_geoc,cross_lon1,thelat2_geoc,cross_lon2,az,baz,del
         d_dist = del * d2km
         dist_sum = dist_sum + d_dist
         del_idx = del_idx + 1

         ref_glat = (g_lat1 + g_lat2)/2d0
         ref_glon = (g_lon1 + g_lon2)/2d0
         !write(*,*) 'ref_glat,ref_glon',ref_glat,ref_glon
         write(cha_lon,'(f5.1)') ref_glon
         write(cha_lat,'(f4.1)') ref_glat
         !write(*,*) 'reflat:',cha_lat,' reflon:',cha_lon,' del',d_del, &
         !&   'rad ',del*d2r 
         write(*,*) cross_lat1,cross_lon1,cross_lat2,cross_lon2, &
     &              ref_glat,ref_glon,del*d2r,d_dist,dist_sum,del_idx
         write(2,*) cross_lat1,cross_lon1,cross_lat2,cross_lon2, &
     &   ref_glat,ref_glon,del*d2r,d_dist,dist_sum,del_idx

         ! To detect mismatch gcp distance from obspy
         diff = dist_obspy_sum - dist_sum
         if (abs(diff)>5d0) then
            write(*,*) 'Detect too large mismatch of dist form obspy'
            write(*,*) dist_obspy_sum, dist_sum, diff
            stop
         endif

         do k = 1, rcvs_num
            if (cross_lat2==rcvs_lat(k) .and. cross_lon2==rcvs_lon(k)) then
                write(*,*) rcvs_name(k)
                write(*,*) 'del_idx, dist_sum',del_idx, dist_sum
                del_idx = 0
                dist_sum = 0d0
                dist_obspy_sum = 0d0 
            endif
         enddo
      enddo
      close(1)
      close(2)


      !allocate(target_array(latlon_num))
      !do i = 1, latlon_num
      !   lon = rcv_latlon_list(1,i)
      !   lat = rcv_latlon_list(2,i)
      !   write(cha_lon,'(f5.1)') lon
      !   write(cha_lat,'(f4.1)') lat
      !   target_dir = 'eifp_'//trim(cha_lon)//'_'//trim(cha_lat)
      !   target_path = trim(target_dir)//'/'//'correct_'//trim(cha_lon)//'_'//trim(cha_lat)
      !   target_array(i) = trim(target_path)
      !   
         !write(*,*) target_path
      !enddo

      !allocate(lminS(0:nmx))
      !allocate(lminT(0:nmx))
      !do i = 1, latlon_num
      !  do STmode = 1, 2
      !    if (STmode == 1) then
      !       cmd = 'S'
      !    else 
      !       cmd = 'T'
      !    endif
      !    path_eip = trim(target_array(i))//'/'//'eigen_pathav_'//cmd
      !    write(*,*) path_eip

      !    open(1,file=trim(path_eip),status='unknown',form='unformatted')
      !       read(1) nord1, nord2, &
      !& (lminS(n),lmaxS(n),n=nord1,nord2)
      !       if (nord2>nmx) then
      !           nord2 = nmx
      !           write(*,*) 'nmx adjustment is applied..'
      !       endif
      !    close(1)
      !  enddo
      !enddo
      
      !deallocate(target_array)

      endprogram calc_divdel



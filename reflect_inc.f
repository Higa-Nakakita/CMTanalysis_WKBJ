
      subroutine ref_inc
      implicit none
      include 'wkbj_common.inc'

      ! ==== compute constant number (wkbj_common.inc) ======
      pi  = 4d0 * datan(1d0)
      pi2 = 2d0 * pi
      r2d = 180d0 / pi
      d2r = pi / 180d0
      dkm2d = r2d / 6371d0 ! convert km  -> deg on the earth
      d2km  = d2r * 6371d0 ! convert deg -> km  on the earth
      dkm2r = 1d0 / 6371d0 ! convert km -> rad  on the earth

      endsubroutine ref_inc

!  input:
!    x(npts) = input data array
!    y(jpts) = output data array
!    j1,j2   = first and last points of actutal data
!    nl,nr = left and right taper width
      subroutine cos_tap(x,y,j1,j2,nl,nr)
      implicit none
      real(4) x(1),y(1),pi,arg
      integer i,j1,j2,nl,nr,jnew
      pi = 3.141592654
!     taper left side
      jnew = 0
      do i = j1, j1+nl-1
          jnew = jnew + 1
          arg = pi * float(i+1-j1-nl) / float(nl)
          y(jnew) = x(i) * (1+cos(arg))/2
      enddo
!     re-arrange for flat portion
      do i = j1+nl, j2-nr
          jnew = jnew + 1
          y(jnew) = x(i)
      enddo
!     taper right side
      do i = j2-nr+1, j2
          jnew = jnew + 1
          y(jnew) = x(i) * (1+cos(arg))/2
      enddo
      return
      end

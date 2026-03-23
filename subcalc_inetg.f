        real(8) function integra(N,a,b,xarray_in,yarray_in)
        implicit none
        integer i, N, N_integ
        integer, parameter :: nmax_integ=100000
        real(8) a, b, dx, tmp, tmp2
        real(8) xarray_in(1:N), yarray_in(1:N)
        real(8) xarray(nmax_integ),yarray(nmax_integ)
        real(8) dyarray(nmax_integ)
        real(8) xtmp, xtmp_old, f, df, ddf
        real(8) trapezoid_integ
        real(8) x_ck1, x_ck2


        !write(*,*) 'x',xarray_in
        !write(*,*) 'y',yarray_in
        ! check
        x_ck2 = -1d0
        do i = 1, N
           x_ck1 = xarray_in(i)
           if (x_ck1 <= x_ck2) then
               write(*,*) 'Invalid xarray_in'
               write(*,*) x_ck1,x_ck2
               stop
           endif
           x_ck2 = x_ck1
        enddo


        ! =========================================
        ! By using akima interpolation
        ! make ideal xarray and yarray using integrand 
        xarray = 0d0
        yarray = 0d0

        N_integ = 1000
        dx = (b-a)/dble(N_integ)
        ! make divided xarray with accuracy of dx
        if (N_integ+1>nmax_integ) then
            write(*,*) 'too large N_integ',N_integ
            write(*,*) 'Increase the value of dx.'
            write(*,*) 'Please sacrifice the accuracy.'
            write(*,*) 'stop'
            stop
        endif

        !do i = 1, N_integ+1
        !    xtmp = a + dx * dble(i-1)
        !    xarray(i) = xtmp
        !enddo
        do i = 1, N_integ+1
             xtmp = a + dx * dble(i-1)
             xarray(i) = xtmp
        enddo

        ! make approx continuous function for integration
        xtmp_old = -1d0
        !write(*,*) xarray_in
        !write(*,*) yarray_in
        do i = 1, N_integ+1
            xtmp = xarray(i)
            if (xtmp <= xtmp_old) then
                write(*,*) 'Invalid xaxis',xtmp, xtmp_old
                stop
            endif
            xtmp_old = xtmp
            call makima1p_modi(f,df,ddf,xtmp,N,yarray_in,xarray_in)
            yarray(i) = f
            if (f/=f) then
                write(*,*) 'detect NaN'
                write(*,*) i,xtmp,f
                stop
            endif
            dyarray(i)= df
        enddo

        !open(2,file='interpo.txt',action='write')
        !do i = 1, N_integ+1
        !    write(2,*) xarray(i), yarray(i), dyarray(i)
        !enddo 
        !close(2)
        !stop
        ! ==========================================

        !    tmp = trapezoid_integ(dx,N_integ+1,a,b,
        ! &  xarray(1:N_integ+1),yarray(1:N_integ+1))

        !write(*,*) 'dx,N_integ+1,a,b',dx,N_integ+1,a,b
        !write(*,*) 'x',xarray(N_integ-100:N_integ+1)
        !write(*,*) 'y',yarray(N_integ-100:N_integ+1)

        tmp = trapezoid_integ(dx,N_integ+1,a,b,
     &  xarray(1:N_integ+1),yarray(1:N_integ+1))

        integra = tmp
        !write(*,*) 'subroutine integra',integra

        endfunction


        real(8) function trapezoid_integ(dx,N,a,b,xn,yn)
        implicit none
        integer i, N
        real(8) dx, a, b, sigma, sigma2, x1, x2, y1, y2
        real(8) xn(N), yn(N)
        
        sigma = 0d0
        sigma2 = 0d0
        do i = 1, N-1
           x1 = xn(i)
           x2 = xn(i+1)
           y1 = yn(i)
           y2 = yn(i+1)
           sigma = sigma + (y2+y1)
        enddo
        sigma2 = 0.5d0 * sigma * dx
        trapezoid_integ = sigma2

        endfunction



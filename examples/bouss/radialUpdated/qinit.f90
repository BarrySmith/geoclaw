
subroutine qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    
    use qinit_module, only: qinit_type,add_perturbation
    use geoclaw_module, only: sea_level, grav
    use starting_module, only: tstart, mx_starting_radial, &
                           x_radial,eta_radial,u_radial
    
    implicit none
    
    ! Subroutine arguments
    integer, intent(in) :: meqn,mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    
    ! Locals
    integer :: i,j,m
    real(kind=8) :: x0,y0,x,y,dx0,dy0,dsigma,r,Rearth,pi,eta,hU
    real(kind=8) :: cos_theta,sin_theta,a1,a2,u,theta,xx,yy
    integer :: k


    Rearth = 6367.5d3
    pi = 4.d0*atan(1.d0)
    
    ! Set flat state based on sea_level
    q = 0.d0
    forall(i=1:mx, j=1:my)
        q(1,i,j) = max(0.d0, sea_level - aux(1,i,j))
    end forall

    !write(6,*) '+++ mx_starting_radial = ', mx_starting_radial
    !write(6,*) '+++ x_radial(mx_starting_radial) = ',x_radial(mx_starting_radial)
    !write(57,*) '+++ qinit:'

    x0 = 0.d0
    y0 = 0.d0

    do i=1,mx
      x = (xlower + (i-0.5d0)*dx)
      dx0 = x - x0
      do j=1,my
          if (q(1,i,j) > 0.d0) then
              y = (ylower + (j-0.5d0)*dy)
              dy0 = y - y0
              !dsigma = 2.d0 * asin(sqrt(sin(0.5d0*dy0)**2 + cos(y0) * cos(y) &
              !          * sin(0.5*dx0)**2))
              !r = Rearth * dsigma
              r = sqrt(dx0**2 + dy0**2)  ! cartesian
              !eta = 0.d0
              !u = 0.d0

              !eta = 5.d0*exp(-((r-50.d3)/20.d3)**2)
              !eta = 2.d0*exp(-(r/20.d3)**2)
              !write(6,*) '+++ x,y,r: ',x,y,r

              if (r < x_radial(mx_starting_radial-1)) then
                   k = 2
                   eta = 0.d0
                   do while ((r>x_radial(k)) .and. (k<mx_starting_radial))
                      k = k+1
                      enddo
                   !write(57,*) dx0*180.d0/pi, dy0*180.d0/pi, dsigma,r,k
                   a1 = (x_radial(k) - r)/(x_radial(k) - x_radial(k-1))
                   a2 = (r - x_radial(k-1))/(x_radial(k) - x_radial(k-1))
                   eta = a1*eta_radial(k-1) + a2*eta_radial(k)
                   u = a1*u_radial(k-1) + a2*u_radial(k)
                   q(1,i,j) = q(1,i,j) + eta
                   endif
            
              hU = q(1,i,j) * u                 ! from input file

              !hU = sqrt(grav*q(1,i,j)) * eta     ! from SWE eigenvector
              !hU = 0.d0   ! stationary

              yy = y-y0
              xx = x-x0
              theta = atan2(yy,xx)
              q(2,i,j) = hU * cos(theta)
              q(3,i,j) = hU * sin(theta)

              endif
          enddo
       enddo

    
end subroutine qinit


subroutine setprob

    !> Copy this file to your directory and modify to set up problem
    !! parameters or read other data.
    !!
    !! This default version is for Boussinesq terms
    !! Eventually reading these should be moved into bouss_module

    use amr_module, only: outunit
    use bouss_module
    use starting_module, only: tstart, mx_starting_radial, &
                               x_radial,eta_radial,u_radial
    implicit none
    integer iunit,i
    character(len=25) fname

#ifdef WHERE_AM_I
    write(*,*) 'starting setprob'
#endif
    iunit = 7
    fname = 'setprob.data'
!   # open the unit with new routine from Clawpack 4.4 to skip over
!   # comment lines starting with #:
    call opendatafile(iunit, fname)

    read(7,*) ibouss
    if (ibouss .eq. 0) then
       write(*,*)" Using Shallow Water equations"
       write(outunit,*)" Using Shallow Water equations"
    else if (ibouss .eq. 1) then
       write(*,*)" Using Madsen Sorensen equations"
       write(outunit,*)" Using Madsen Sorensen equations"
    else if (ibouss .eq. 2) then
       write(*,*)" Using SGN equations"
       write(outunit,*)" Using SGN equations"
    else
       write(*,*)" Unrecognized option for equation set"
       write(outunit,*)" Unrecognized option for equation set"
       stop
    endif

    if (ibouss .gt. 0) then
      read(7,*) minLevelBouss
      read(7,*) maxLevelBouss
      read(7,*) boussMinDepth
      read(7,*) isolver
      read(7,*) startBoussTime

      write(*,900) minLevelBouss, maxLevelBouss
      write(outunit,900) minLevelBouss, maxLevelBouss
 900   format("==> Applying Bouss equations to selected grids between levels ",i3," and ",i3)

      write(*,*)"==> Use Bouss. in water deeper than ",boussMinDepth
      write(outunit,*)"==> Use Bouss. in water deeper than ",boussMinDepth

      if (isolver .eq. 1) then
         write(*,*)" No longer supporting GMRES solver"
         stop
      else if (isolver .eq. 2) then
#ifdef HAVE_PARDISO
         !write(*,*)" Using Pardiso solver"
         !write(outunit,*)" Using Pardiso solver"
         write(*,*)"Cannot use expired Pardiso solver"
         write(outunit,*)"Cannot use expired Pardiso solver"
         stop
#else
         write(*,*)"need to install Pardiso for this option"
         stop
#endif
      else if (isolver .eq. 3) then
#ifdef HAVE_PETSC
        write(*,*)"Using a PETSc solver"
        write(outunit,*)"Using PETSc solver"
#else
        write(*,*)"need to install PETSc for this option"
        stop
#endif
      else
        write(*,*)"Unknown solver",isolver," choose 1,2 or 3"
        write(outunit,*)"Unknown solver",isolver," choose 1,2 or 3"
        stop
     endif
     close(unit=iunit)


      if (startBoussTime .le. t0) then
         write(*,*)"Using Bouss equations from the start"
         write(outunit,*)"Using Bouss equations from the start"
         startWithBouss = .true.
      else
         write(*,*)"==> Wait until time ",startBoussTime," for before starting Bouss"
         write(*,*)"==> Using SWE until then."
         startWithBouss = .false.
      endif

    endif

    close(unit=iunit)

    open(unit=59, file='starting.data', status='old',form='formatted')
    read(59,*) tstart
    read(59,*) mx_starting_radial
    write(6,*) 'Reading starting.data at tstart = ',tstart

    do i=1,mx_starting_radial
        read(59,*) x_radial(i),eta_radial(i),u_radial(i)
        enddo
    close(unit=59)
              
#ifdef WHERE_AM_I
    write(*,*) 'ending   setprob'
#endif
              

end subroutine setprob


subroutine setprob

    ! for problem specific features, copy to your directory and modify

    use grid_module, only: mx_edge
    use setprob_module, only: eta0,u0
    implicit none
    integer :: mx_cell,i,iunit

    iunit = 51
    
    open(unit=iunit, file='eta_u.txt', status='old', form='formatted')

    mx_cell = mx_edge - 1
    write(6,*) 'Reading eta0,u0 with mx_cell = ',mx_cell
    do i=1,mx_cell
        read(iunit,*) eta0(i),u0(i)
    enddo

    close(unit=iunit)

end subroutine setprob

subroutine sumrho( nfile, coefflist, filelist, matsize, rho )
    implicit none
    integer(kind=4),                     intent(in) :: nfile
    real(kind=8),      dimension(nfile), intent(in) :: coefflist
    character(len=32), dimension(nfile), intent(in) :: filelist
    integer(kind=4),   dimension(3),     intent(in) :: matsize
    real(kind=8),      dimension(matsize(1),matsize(2),matsize(3)), intent(out) :: rho
    real(kind=8),      dimension(matsize(1),matsize(2),matsize(3)) :: irho
    integer(kind=4) :: i
    rho = rho * 0.0
    do i = 1, nfile
        call cube_read_scafield( filelist(i), matsize, irho )
        rho = rho + irho * coefflist(i)
    end do
end subroutine sumrho

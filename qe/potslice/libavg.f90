subroutine calc_dS( a, b, S )    
    implicit none
    real(kind=8), dimension(3), intent(in) :: a, b
    real(kind=8), intent(out) :: S
    real(kind=8) :: la, lb, dotab, cosab, sinab
    la = sqrt( dot_product( a, a ) )
    lb = sqrt( dot_product( b, b ) )
    dotab = dot_product( a, b )
    cosab = dotab / (la*lb)
    sinab = sqrt( 1.0 - cosab * cosab )
    S = la * lb * sinab
end subroutine calc_dS


subroutine calc_integral( matsize, fd, sym, X, dS, integ )
    implicit none
    integer(kind=4), dimension(3), intent(in) :: matsize
    real(kind=8), dimension(matsize(1),matsize(2),matsize(3)), intent(in) :: fd
    integer(kind=4), intent(in) :: sym, X
    real(kind=8), intent(in)  :: dS
    real(kind=8), intent(out) :: integ
    integer(kind=4) :: i, j
    integ = 0.0
    select case (sym)
        case (1)
            do i = 1, matsize(2) - 1
                do j = 1, matsize(3) - 1
                    integ = integ + 0.25 * ( fd( X, i,   j ) + fd( X, i,   j+1 ) &
                                           + fd( X, i+1, j ) + fd( X, i+1, j+1 ) )
                end do
            end do
        case (2)
            do i = 1, matsize(1) - 1
                do j = 1, matsize(3) - 1
                    integ = integ + 0.25 * ( fd( i,   X, j ) + fd( i,   X, j+1 ) &
                                           + fd( i+1, X, j ) + fd( i+1, X, j+1 ) )
                end do
            end do
        case (3)
            do i = 1, matsize(1) - 1
                do j = 1, matsize(2) - 1
                    integ = integ + 0.25 * ( fd( i,   j, X ) + fd( i,   j+1, X ) &
                                           + fd( i+1, j, X ) + fd( i+1, j+1, X ) )
                end do
            end do
        case default
            integ = 0.0
    end select
    integ = integ * dS 
end subroutine calc_integral

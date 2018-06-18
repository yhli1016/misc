subroutine sumldos (ikpoint, iband, iid, pct)

    use shareddata

    implicit none

    integer(kind=4), intent(in)  :: ikpoint, iband, iid
    real(kind=8), intent(out) :: pct
    integer(kind=4) :: i, j
    real(kind=8) :: sumproj, normfactor

    ! check for illegal parameters
    if ( ikpoint < 0 .or. ikpoint > nkpoint ) then
        write (*,*) "ILLEGAL KPOINT INDEX!"
    else if ( iband < 0 .or. iband > nband ) then
        write (*,*) "ILLEGAL BAND INDEX!"
    else if ( (trim(lnoncolin) == "T") .or. (trim(lspinorbit) == "T") ) then
        write (*,*) "SPIN NOT IMPLEMENTED!"
    else
        normfactor = sum( proj(:,ikpoint,iband) )
        sumproj = 0.0
        ! loop over atomic orbitals
        do i = 1, nwfc
            if (atomid(i) == iid) then
                sumproj = sumproj + proj( i, ikpoint, iband )
            end if
        end do
        pct = pct + sumproj / normfactor
    end if

end subroutine sumldos

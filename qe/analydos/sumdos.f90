module sumdos
    use projection
    implicit none

    contains

    subroutine sumldos(ikpoint, iband, iid, pct)
        implicit none

        integer(kind=4), intent(in)  :: ikpoint, iband, iid
        real(kind=8), intent(out) :: pct
        integer(kind=4) :: i, j
        real(kind=8) :: sumproj, normfactor
    
        ! check for illegal parameters
        if (ikpoint < 0 .or. ikpoint > nkpoint) then
            write (*,*) "ILLEGAL KPOINT INDEX!"
        else if (iband < 0 .or. iband > nband) then
            write (*,*) "ILLEGAL BAND INDEX!"
        else if ((trim(adjustl(lnoncolin)) == "T") .or. &
                 (trim(adjustl(lspinorbit)) == "T")) then
            write (*,*) "SPIN NOT IMPLEMENTED!"
        else
            normfactor = sum(proj(:,ikpoint,iband))
            sumproj = 0.0
            ! loop over atomic orbitals
            do i = 1, nwfc
                if (atomid(i) == iid) then
                    sumproj = sumproj + proj(i, ikpoint, iband)
                end if
            end do
            pct = pct + sumproj / normfactor
        end if
    end subroutine

    subroutine sumpdos(ikpoint, iband, atomic_species, nlm_cmp, pct)
        implicit none

        integer(kind=4),  intent(in)  :: ikpoint, iband
        character(len=2), intent(in)  :: atomic_species
        integer(kind=4), dimension(3), intent(in)  :: nlm_cmp
        real(kind=8),     intent(out) :: pct
        integer(kind=4) :: i, j
        real(kind=8) :: sumproj, normfactor
    
        ! check for illegal parameters
        if (ikpoint < 0 .or. ikpoint > nkpoint) then
            write (*,*) "ILLEGAL KPOINT INDEX!"
        else if (iband < 0 .or. iband > nband) then
            write (*,*) "ILLEGAL BAND INDEX!"
        else if ((trim(adjustl(lnoncolin)) == "T") .or. &
                 (trim(adjustl(lspinorbit)) == "T")) then
            write (*,*) "SPIN NOT IMPLEMENTED!"
        else
            normfactor = sum(proj(:,ikpoint,iband))
            sumproj = 0.0
            ! loop over atomic orbitals
            do i = 1, nwfc
                ! check the species and nlm of this atomic orbital
                if ((trim(adjustl(atomic_species)) == &
                     trim(adjustl(atomsym(i))))  &
                    .and. (nlm_cmp(1) == nlm(i,1)) &
                    .and. (nlm_cmp(2) == nlm(i,2)) &
                    .and. (nlm_cmp(3) == nlm(i,3))) then
                    sumproj = sumproj + proj(i, ikpoint, iband)
                end if
            end do
            pct = pct + sumproj / normfactor
        end if
    end subroutine
end module

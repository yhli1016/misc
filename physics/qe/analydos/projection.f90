!===============================================================================
!
! Common structure of trimmed proj.dat for QE < 6.x
!
! 1: number of atomic wavefunctions, number of kpoints, number of bands
! 2: is non-collinear magnetization enabled, is spin-orbit interaction enabled
! 3: atomic wavefunction id, atomid+atomsymbol, n, l, m
! 4: ikpoint, iband, projection
! ... ...    ... ...
!
! For QE >= 6.x, the 3rd line becomes:
! atomic wavefunction id, atomid, atomsymbol, wfcsymbol, n, l, m
!
!===============================================================================

module projection
    implicit none

    integer(kind=4)  :: nwfc, nkpoint, nband
    character(len=1) :: lnoncolin, lspinorbit
    integer(kind=4),  allocatable, dimension(:)      :: wfcid, atomid
    character(len=3), allocatable, dimension(:)      :: atomsym
#ifdef QE6
    character(len=2), allocatable, dimension(:)      :: wfcsym
#endif
    integer(kind=4),  allocatable, dimension(:,:)    :: nlm
    real(kind=8),     allocatable, dimension(:,:,:)  :: proj

    contains

    subroutine load_projection(projfn)
        implicit none

        character(len=32), intent(in) :: projfn
        integer(kind=4) :: i, j, k, tempa, tempb

        open(unit=7, file=projfn, status="old")
        read (7,*) nwfc, nkpoint, nband
        read (7,*) lnoncolin, lspinorbit
        allocate(wfcid(nwfc))
        allocate(atomid(nwfc))
        allocate(atomsym(nwfc))
#ifdef QE6
        allocate(wfcsym(nwfc))
#endif
        allocate(nlm(nwfc,3))
        allocate(proj(nwfc,nkpoint,nband))
        do i = 1, nwfc
#ifdef QE6
            read (7,"(2I5,2A4,I6,2I5)") wfcid(i), atomid(i), atomsym(i), &
                                        wfcsym(i), nlm(i,:)
#else
            read (7,"(2I5,A3,3I5)") wfcid(i), atomid(i), atomsym(i), nlm(i,:)
#endif
            do j = 1, nkpoint
                do k = 1, nband
                     read (7,*) tempa, tempb, proj(i,j,k)
                end do
            end do
        end do
        close(unit=7)
    end subroutine

    subroutine free_projection()
        implicit none

        deallocate(wfcid)
        deallocate(atomid)
        deallocate(atomsym)
#ifdef QE6
        deallocate(wfcsym)
#endif
        deallocate(nlm)
        deallocate(proj)
    end subroutine
end module

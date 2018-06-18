!===============================================================================
!
! Common structure of trimmed proj.dat
!
! number of atomic wavefunctions, number of kpoints, number of bands
! is non-collinear magnetization enabled, is spin-orbit interaction enabled
! atomic wavefunction id, atomid+atomsymbol, n, l, m
! ikpoint, iband, projection
! ... ...    ... ...
!
!===============================================================================

module shareddata

    implicit none

    integer(kind=4)  :: nwfc, nkpoint, nband
    character(len=1) :: lnoncolin, lspinorbit
    integer(kind=4),  allocatable, dimension(:)      :: wfcid, atomid
    character(len=3), allocatable, dimension(:)      :: atomsym
    integer(kind=4),  allocatable, dimension(:,:)    :: nlm
    real(kind=8),     allocatable, dimension(:,:,:)  :: proj

end module shareddata

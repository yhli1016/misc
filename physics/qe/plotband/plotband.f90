!===============================================================================
!
! This program reads in the outpout of bands.x and converts it to a plotable
! format.
!
! The length unit of k-path is in 1/bohr. This is particlularly convenient
! for calculating effective mass.
!
! This program does not require any input files. All needed parameters are read
! from standard input.
!
!===============================================================================


program plotband

    implicit none
    
    !
    ! controlling flags
    !
    character(len=32)                         :: iname, oname
    integer(kind=4)                           :: ounit
    real(kind=8)                              :: efermi, alat
    !
    ! band structure
    !
    integer(kind=4)                           :: nbnd, nks
    real(kind=8), allocatable, dimension(:,:) :: kpoints, bndstr
    real(kind=8), allocatable, dimension(:)   :: kpath
    real(kind=8),              dimension(3)   :: k1, k2
    real(kind=8)                              :: angle
    !
    ! constants
    !
    real(kind=8),              parameter      :: PI = 3.141592653589793
    real(kind=8),              parameter      :: HAR = 27.21138506
    !
    ! loop counters
    !
    integer(kind=4)                           :: i, j
    !
    ! namelist corresponding to the header of iname
    !
    namelist /plot/ nbnd, nks


    !
    ! read parameters from stdin
    !
    write (*,*) "Give the filename of input file"
    read  (*,*) iname
    write (*,*) "Give the finename of output file"
    read  (*,*) oname
    write (*,*) "Give the unit for energy ( 1-eV, 2-Hartree )"
    read  (*,*) ounit
    write (*,*) "Give the fermi energy in eV"
    read  (*,*) efermi
    write (*,*) "Give alat in bohr"
    read  (*,*) alat

    !
    ! read bandstructure file
    !
    open ( unit=7, file=iname, status="old" )
    read ( 7, nml=plot )
    allocate( kpoints( nks, 3    ) )
    allocate(  bndstr( nks, nbnd ) )
    allocate(   kpath( nks       ) )
    do i = 1, nks
        read (7,*) kpoints(i,:)
        read (7,*) bndstr(i,:)
    end do
    close ( unit=7 )
    
    !
    ! units conversion
    !
    kpoints = 2.0 * PI  / alat * kpoints
    if ( ounit == 2 ) then
        bndstr  = 1.0 / HAR * bndstr
        efermi  = 1.0 / HAR * efermi
    end if

    !
    ! determine kpath
    !
    kpath(1) = 0.0
    do i = 2, nks
        kpath(i) = kpath(i-1) + sqrt( dot_product( kpoints(i,:)-kpoints(i-1,:), kpoints(i,:)-kpoints(i-1,:) ) )
    end do
    
    
    !
    ! determine high symmetry points
    !
    ! algorithm from plotband.f90 of PWSCF
    !
    ! the first and the last kpoints are always high symmetric points
    !
    write (*,"(A36)") "list of high symmetric points"
    write (*,"(A44)") "---------------------------------------------"
    write (*,"(A4,4A10)") "i", "x", "kx", "ky", "kz"
    write (*,"(I4,4F10.6)") 1, kpath(1), kpoints(1,:)
    do i = 2, nks-1
        k1 = kpoints(i,:) - kpoints(i-1,:)
        k2 = kpoints(i+1,:) - kpoints(i,:)
        angle = dot_product( k1, k2 ) / sqrt(dot_product(k1,k1)*dot_product(k2,k2))
        if ( abs(angle-1.0) > 1.0e-4 ) then
            write (*,"(I4,4F10.6)") i, kpath(i), kpoints(i,:)
        end if
    end do
    write (*,"(I4,4F10.6)") nks, kpath(nks), kpoints(nks,:)
    write (*,"(A44)") "---------------------------------------------"
    
    !
    ! output
    !
    open ( unit=8, file=oname, status="replace" )
    do i = 1, nks
        write (8,"(F10.6)",advance="no") kpath(i)
        do j = 1, nbnd
            write (8,"(F12.6)",advance="no") bndstr(i,j) - efermi
        end do
        write (8,*)
    end do
    close ( unit=8 )
    
    
    !
    ! clean up
    !
    deallocate( kpoints )
    deallocate(  bndstr )
    deallocate(   kpath )
    
end program plotband

!===============================================================================
!
! REMEMBER TO REMOVE THE HEAD LINES BEFORE INVOKING THIS PROGRAM.
!
! THE PROJECTION FILE TO ANALYZE MUST BE IN THE FOLLOWING FORM FOR QE 5.x:
!     234       1     144
!    F    F
!    1    1C      1    0    1
!       1       1        0.0498992229
!       1       2        0.0000003860
!       1       3        0.0201991768
!       1       4        0.0000000008
!       1       5        0.0026582018
!       1       6        0.0456243229
!       1       7        0.0000000000
!       1       8        0.0005048665
!       1       9        0.0000052600
!
! FOR QE 6.x IT SHOULD BE:
!
!      64     121      16
!    F    F
!    1    1  In  5s     1    0    1
!       1       1        0.1095559867
!       1       2        0.0554159126
!       1       3        0.2140375616
!
! This program summates the projection of specified states \Psi_{ikpoint,iband}
! on the atomic states with specified atomic species and quantum number. It is
! helpful when calculating the contribution of specified atomic states to given
! state.
!
! This program requires two CLI parameters specifying the input and output file.
! Inputfile must contain the the following lines:
!   proj.dat       // filename of projection file
!   ikmin ikmax    // min and max of kpoint index to analyze
!   ibmin ibmax    // min and max of band   index to analyze
!   ntype           // number of atomic state types
!   species(1) n(1) l(1) m(1)  // atomic species and quantum number taken into
!                              // consideration
!   ... ...  ... ...
!
!===============================================================================

program analypdos
    use projection, only: load_projection, free_projection
    use sumdos,     only: sumpdos
    implicit none

    ! file names of input file, projection file and output file
    character(len=32) :: inpfn, projfn, outfn
    ! flags in the input file
    integer(kind=4) :: ikmin, ikmax, ibmin, ibmax, ntype
    character(len=2), allocatable :: atomic_species(:)
    integer(kind=4), allocatable :: nlm_inp(:,:)
    ! contribution collected from listed atomic species in percentage
    real(kind=8) :: pct
    ! loop counters
    integer(kind=4) :: i, j, k

    ! parse cli-parameters and read input
    call get_command_argument(1, inpfn)
    call get_command_argument(2, outfn)
    open(unit=7, file=inpfn, status="old")
    read (7,*) projfn
    read (7,*) ikmin, ikmax
    read (7,*) ibmin, ibmax
    read (7,*) ntype
    allocate(atomic_species(ntype))
    allocate(nlm_inp(ntype,3))
    do i = 1, ntype
        read (7,*) atomic_species(i), nlm_inp(i,:)
    end do
    close(unit=7)

    ! read projection file. 
    call load_projection(projfn)

    ! collect projection from specified atomic species
    open(unit=7, file=outfn, status="replace")
    write (7,"(2A8,A14)") "#  ikpnt", "iband", "Proj"
    do i = ikmin, ikmax
        do j = ibmin, ibmax
            pct = 0.0
            do k = 1, ntype
                call sumpdos(i, j, atomic_species(k), nlm_inp(k,:), pct)
            end do
            write (7,"(2I8,F14.9)") i, j, pct 
        end do
    end do
    close(unit=7)

    ! clean up
    deallocate(atomic_species)
    deallocate(nlm_inp)
    call free_projection()
end program

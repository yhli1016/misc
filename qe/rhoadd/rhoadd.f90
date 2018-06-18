!===============================================================================
!
! This program calculates the linear combination of scalar fields in cube files.
!
! This program is dedicated to calculating differential charge density.
!
! This program requires a CLI parameter specifying the input file which must be 
! in the following format.
!
! output.cube       ! filename of output file
! 3                 ! number of input files
!  1.0 file1.cube   ! coefficent(i), input filename(i)
! -1.0 file2.cube
! -1.0 file3.cube
!
!===============================================================================


program rhoadd

    implicit none

    !
    ! input
    !
    character(len=32)                                :: inpname, outname
    integer(kind=4)                                  :: nfile
    real(kind=8),      allocatable, dimension(:)     :: coefflist
    character(len=32), allocatable, dimension(:)     :: filelist
    !
    ! cube file information
    !
    character(len=32),              dimension(2)     :: header
    real(kind=8),                   dimension(3)     :: mesh_origin
    integer(kind=4),                dimension(3)     :: mesh_size
    real(kind=8),                   dimension(3,3)   :: mesh_unit
    integer(kind=4)                                  :: natom
    integer(kind=4),   allocatable, dimension(:)     :: atom_number
    real(kind=8),      allocatable, dimension(:)     :: atom_charge
    real(kind=8),      allocatable, dimension(:,:)   :: atom_pos
    real(kind=8),      allocatable, dimension(:,:,:) :: rho
    !
    ! meshgrid consitent flag
    !
    logical                                          :: iconsist
    !
    ! loop counter
    !
    integer(kind=4)                                  :: i


    !
    ! read input
    !
    call GET_COMMAND_ARGUMENT( 1, inpname )
    write (*,*) "Reading input from: ", inpname
    write (*,*)
    open  (unit = 8, file = inpname, status = "old")
    read  (8,*) outname
    read  (8,*) nfile
    allocate( coefflist(nfile) )
    allocate(  filelist(nfile) )
    do i = 1,nfile
        read (8,*) coefflist(i), filelist(i)
    end do
    close (unit = 8)

    !
    ! check meshgrid
    !
    write (*,*) "Checking meshgrid"
    write (*,*)
    call cube_check_meshinfo( nfile, filelist, 1.0E-6, .true., iconsist )

    !
    ! if iconsist == .true., continue computing; else quit
    !
    if ( iconsist ) then
        !
        ! set header
        !
        header(1) = "Combined charge density"
        header(2) = "Generated using rhoadd code"

        !
        ! read natom from cube file
        !
        call cube_read_natom( filelist(1), natom )

        !
        ! read mesh info
        !
        call cube_read_meshinfo( filelist(1), mesh_origin, mesh_size, mesh_unit )

        !
        ! read atoms info from cube file
        !
        allocate( atom_number(natom) )
        allocate( atom_charge(natom) )
        allocate( atom_pos(natom,3)  )
        call cube_read_atominfo( filelist(1), natom, atom_number, atom_charge, atom_pos )

        !
        ! allocate rho using meshinfo
        !
        allocate( rho(mesh_size(1),mesh_size(2),mesh_size(3) ) )

        !
        ! call sumrho
        !
        write (*,*) "Summing up cube files"
        write (*,*)
        call sumrho( nfile, coefflist, filelist, mesh_size, rho )

        !
        ! writing cube
        !
        write (*,*) "Writing final cube file"
        write (*,*)
        call cube_write( outname, header, natom, mesh_origin, mesh_size, mesh_unit, &
                        atom_number, atom_charge, atom_pos, rho )

        !
        ! deallocate
        !
        deallocate( rho )
        deallocate( atom_number )
        deallocate( atom_charge )
        deallocate( atom_pos    )
    else
        write (*,*) "ERROR: MESHGRID INCONSISTENT! CHECK YOUR CUBE FILES!"
        write (*,*)
    end if

    !
    ! deallocate
    !
    deallocate( coefflist )
    deallocate(  filelist )

end program rhoadd

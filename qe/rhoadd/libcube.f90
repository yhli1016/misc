!===============================================================================
!
! This file defines the following subroutines that are used to manipulate Gauss-
! ian98 cube files:
!
! (1) cube_read_natom
!     This subroutine reads number of atoms from cube file.
!
! (2) cube_read_meshinfo
!     This subroutine reads the position of (0,0,0) point, size and dx along ea-
!     ch direction of the grid from cube file.
!
! (3) cube_read_atominfo
!     This subroutine reads atom Zval, charge and position from cube file.
!
! (4) cube_read_scafield
!     This subroutine reads scalar field from cube file.
!
! (5) cube_write
!     This subroutine writes cube file.
!
! (6) cube_check_meshinfo
!     This subrtouine checkes if the meshgrid of cube files are consistent.
!
! For their interfaces, see the definitions below.
!
! Changelog:
!
! 20150927: code style unified
!
! 20150928: added subroutine cube_check_meshgrid
!
!===============================================================================


subroutine cube_read_natom( filename, natom )

    implicit none

    character(len=32), intent(in)   :: filename
    integer(kind=4),   intent(out)  :: natom
    character(len=32), dimension(2) :: header

    open ( unit=8, file=filename, status="old" )
    read (8,*) header(1)
    read (8,*) header(2)
    read (8,*) natom
    close ( unit=8 )

end subroutine cube_read_natom


subroutine cube_read_meshinfo( filename, mesh_origin, mesh_size, mesh_unit )

    implicit none
    
    character(len=32),               intent(in)  :: filename
    real(kind=8),    dimension(3),   intent(out) :: mesh_origin
    integer(kind=4), dimension(3),   intent(out) :: mesh_size
    real(kind=8),    dimension(3,3), intent(out) :: mesh_unit
    character(len=32), dimension(2) :: header
    integer(kind=4) :: natom
    integer(kind=4) :: i

    open ( unit=8, file=filename, status="old" )
    read (8,*) header(1)
    read (8,*) header(2)
    read (8,*) natom, mesh_origin(:)
    do i = 1, 3
        read (8,*) mesh_size(i), mesh_unit(i,:)
    end do
    close ( unit=8 )

end subroutine cube_read_meshinfo


subroutine cube_read_atominfo( filename, natom, atom_number, atom_charge, atom_pos )

    implicit none

    character(len=32),                   intent(in)  :: filename
    integer(kind=4),                     intent(in)  :: natom
    integer(kind=4), dimension(natom),   intent(out) :: atom_number
    real(kind=8),    dimension(natom),   intent(out) :: atom_charge
    real(kind=8),    dimension(natom,3), intent(out) :: atom_pos
    character(len=32), dimension(2) :: header
    integer(kind=4) :: shader_natom
    real(kind=8), dimension(3) :: mesh_origin
    integer(kind=4), dimension(3) :: mesh_size
    real(kind=8), dimension(3,3) :: mesh_unit
    integer(kind=4) :: i

    open ( unit=8, file=filename, status="old" )
    read (8,*) header(1)
    read (8,*) header(2)
    read (8,*) shader_natom, mesh_origin(:)
    do i = 1, 3
        read (8,*) mesh_size(i), mesh_unit(i,:)
    end do
    do i = 1, natom
        read (8,*) atom_number(i), atom_charge(i), atom_pos(i,:)
    end do
    close ( unit=8 )

end subroutine cube_read_atominfo


subroutine cube_read_scafield( filename, mesh_size, scafield )

    implicit none

    character(len=32),                                                  intent(in)  :: filename
    integer(kind=4), dimension(3),                                      intent(in)  :: mesh_size
    real(kind=8),    dimension(mesh_size(1),mesh_size(2),mesh_size(3)), intent(out) :: scafield
    character(len=32), dimension(2)   :: header
    integer(kind=4)                   :: natom
    real(kind=8),      dimension(3)   :: mesh_origin
    integer(kind=4),   dimension(3)   :: shader_mesh_size
    real(kind=8),      dimension(3,3) :: mesh_unit
    integer(kind=4)                   :: atom_number_i
    real(kind=8)                      :: atom_charge_i
    real(kind=8),      dimension(3)   :: atom_pos_i
    integer(kind=4) :: i, j

    open ( unit=8, file=filename, status="old" )
    read (8,*) header(1)
    read (8,*) header(2)
    read (8,*) natom, mesh_origin(:)
    do i = 1, 3
        read (8,*) shader_mesh_size(i), mesh_unit(i,:)
    end do
    do i = 1, natom
        read (8,*) atom_number_i, atom_charge_i, atom_pos_i(:)
    end do
    do i = 1, mesh_size(1)
        do j = 1, mesh_size(2)
            read (8,*) scafield(i,j,:)
        end do
    end do
    close ( unit=8 )

end subroutine cube_read_scafield


subroutine cube_write( filename, header, natom, mesh_origin, mesh_size, mesh_unit, &
                       atom_number, atom_charge, atom_pos, scafield )

    implicit none

    character(len=32),                     intent(in) :: filename
    character(len=32), dimension(2),       intent(in) :: header
    integer(kind=4),                       intent(in) :: natom
    real(kind=8),      dimension(3),       intent(in) :: mesh_origin
    integer(kind=4),   dimension(3),       intent(in) :: mesh_size
    real(kind=8),      dimension(3,3),     intent(in) :: mesh_unit
    integer(kind=4),   dimension(natom),   intent(in) :: atom_number
    real(kind=8),      dimension(natom),   intent(in) :: atom_charge
    real(kind=8),      dimension(natom,3), intent(in) :: atom_pos
    real(kind=8),      dimension(mesh_size(1),mesh_size(2),mesh_size(3)), intent(in) :: scafield
    integer(kind=4) :: i, j, k

    open ( unit=8, file=filename, status="replace" )
    write (8,*) header(1)
    write (8,*) header(2)
    write (8,"(I5,3F12.6)") natom, mesh_origin(:)
    do i = 1,3
        write (8,"(I5,3F12.6)") mesh_size(i), mesh_unit(i,:)
    end do
    do i = 1,natom
        write (8,"(I5,4F12.6)") atom_number(i), atom_charge(i), atom_pos(i,:)
    end do
    do i = 1, mesh_size(1)
        do j = 1, mesh_size(2)
            do k = 1, mesh_size(3)
                write(8,"(ES13.5)",advance="no") scafield(i,j,k)
                if (mod(k,6)==0.or.k==mesh_size(3)) then
                    write(8,*)
                end if
            end do
        end do
    end do
    close ( unit=8 )

end subroutine cube_write


subroutine cube_check_meshinfo( nfile, filelist, thr, iverbose, iconsist )

    implicit none

    integer(kind=4),                     intent(in)  :: nfile
    character(len=32), dimension(nfile), intent(in)  :: filelist
    real(kind=8),                        intent(in)  :: thr
    logical,                             intent(in)  :: iverbose
    logical,                             intent(out) :: iconsist
    real(kind=8),    dimension(3)   :: origin_ref, origin_comp, origin_diff
    integer(kind=4), dimension(3)   :: size_ref,   size_comp, size_diff
    real(kind=8),    dimension(3,3) :: unit_ref,   unit_comp, unit_diff
    integer(kind=4)                 :: i, j


    if (iverbose) then
        write (*,*) "--------------------------------------------------"
        write (*,*) "Verbose output from subroutine cube_check_meshinfo"
        write (*,*)
    end if

    !
    ! initialize iconsist and read in reference file
    ! filelist(1) is always treated as the reference
    !
    call cube_read_meshinfo( filelist(1), origin_ref, size_ref, unit_ref )
    iconsist = .true.

    do i = 2, nfile

        if (iverbose) then
            write (*,*) "Checking ", filelist(i), " against ", filelist(1)
            write (*,*)
        end if

        call cube_read_meshinfo( filelist(i), origin_comp, size_comp, unit_comp )

        !
        ! calculate difference terms
        !
        origin_diff = abs( origin_ref - origin_comp )
        size_diff = size_ref - size_comp
        unit_diff = abs( unit_ref - unit_comp )

        !
        ! check mesh_size
        !
        do j = 1, 3
            if ( size_diff(j) .ne. 0 ) then
                iconsist = .false.
                if (iverbose) then
                    write (*,*) "ERROR! mesh_size is inconsistent!"
                    write (*,*) "mesh_size in ", filelist(1), " is ", size_ref(:)
                    write (*,*) "mesh_size in ", filelist(i), " is ", size_comp(:)
                    write (*,*)
                end if
            end if
        end do

        !
        ! check mesh_origin
        !
        if ( sum(origin_diff) > thr ) then
            iconsist = .false.
            if (iverbose) then
                write (*,*) "ERROR! mesh_origin is inconsistent!"
                write (*,*) "mesh_origin in ", filelist(1), " is ", origin_ref(:)
                write (*,*) "mesh_origin in ", filelist(i), " is ", origin_comp(:)
                write (*,*)
            end if
        end if

        !
        ! check mesh_unit
        !
        if ( sum(unit_diff) > thr ) then
            iconsist = .false.
            if (iverbose) then
                write (*,*) "ERROR! mesh_unit is inconsistent!"
                write (*,*) "mesh_unit in ", filelist(1), " is "
                do j = 1, 3
                    write (*,*) unit_ref(j,:)
                end do
                write (*,*) "mesh_unit in ", filelist(i), " is "
                do j = 1, 3
                    write (*,*) unit_comp(j,:)
                end do
                write (*,*)
            end if
        end if

        !
        ! If iconsist == .true. everything is OK.
        ! Else there is no need to continue checking.
        !
        if (iverbose) then
            if ( iconsist ) then
                write (*,*) "Passed"
                write (*,*)
            else
                write (*,*) "Meshgrid inconsistency detected"
                write (*,*) "Quitting subroutine cube_check_meshinfo"
                write (*,*)
                write (*,*) "--------------------------------------------------"
                return
            end if
        end if

    end do
    
    if (iverbose) then
        write (*,*) "--------------------------------------------------"
    end if

end subroutine cube_check_meshinfo

!===============================================================================
!
! This prpgram calculates the average and integrated scalar field along c axis.
!
! This program is handy when calculating the average potential or integrated
! charge density along z-axis.
!
! This program does not require any input files. All needed parameters are read
! from standard input.
!
!===============================================================================


program potslice

    implicit none

    !
    ! controlling flags
    !
    character(len=32)                           :: cubefile
    character(len=32)                           :: potfile
    real(kind=8)                                :: scale_factor
    integer(kind=4)                             :: ialgo
    !
    ! scalar field
    !
    real(kind=8),              dimension(3)     :: mesh_origin
    integer(kind=4),           dimension(3)     :: mesh_size
    real(kind=8),              dimension(3,3)   :: mesh_unit
    real(kind=8), allocatable, dimension(:,:,:) :: scafield
    !
    ! integrated and averaged potential
    !
    real(kind=8), allocatable, dimension(:)     :: rz, potint, potavg
    real(kind=8)                                :: dz, dS
    !
    ! S is only used when ialgo = 1
    !
    real(kind=8)                                :: S
    !
    ! NGXY is only used for ialgo = 2
    !
    integer(kind=4)                             :: NGXY
    !
    ! loop counters
    !
    integer(kind=4)                             :: i


    !
    ! read parameters from stdin
    !
    write (*,*) "Give the file name of cube file"
    read  (*,*) cubefile
    write (*,*) "Give the file name of pot file"
    read  (*,*) potfile
    write (*,*) "Give the scale factor"
    read  (*,*) scale_factor
    write (*,*) "Give the algorithm"
    read  (*,*) ialgo

    !
    ! array allocation
    !
    call cube_read_meshinfo( cubefile, mesh_origin, mesh_size, mesh_unit )
    allocate( scafield(mesh_size(1),mesh_size(2),mesh_size(3)) )
    allocate(       rz(mesh_size(3)) )
    allocate(   potint(mesh_size(3)) )
    allocate(   potavg(mesh_size(3)) )

    !
    ! read scalar field and rescale
    !
    call cube_read_scafield( cubefile, mesh_size, scafield )
    scafield = scafield * scale_factor

    !
    ! calculate potavg and potint
    !
    dz = sqrt( dot_product( mesh_unit(3,:), mesh_unit(3,:) ) )
    call calc_dS( mesh_unit(1,:), mesh_unit(2,:), dS )
    if ( ialgo == 1 ) then
        S  = dS * (mesh_size(1)-1) * (mesh_size(2)-1)
        do i = 1, mesh_size(3)
            rz(i) = (i-1) * dz
            call calc_integral( mesh_size, scafield, 3, i, dS, potint(i) )
            potavg(i) = potint(i) / S
        end do
    else
        NGXY = mesh_size(1) * mesh_size(2)
        do i = 1, mesh_size(3)
            rz(i) = (i-1) * dz
            potint(i) = sum( scafield(:,:,i) ) * dS
            potavg(i) = sum( scafield(:,:,i) ) / NGXY
        end do
    end if


    !
    ! write potfile
    !
    open ( unit=7, file=potfile, status="replace" )
    do i = 1, mesh_size(3)
        write (7,"(3ES13.5)") rz(i), potint(i), potavg(i)
    end do
    close ( unit=7 )

    !
    ! array deallocation
    !
    deallocate( scafield )
    deallocate(       rz )
    deallocate(   potint )
    deallocate(   potavg )

end program potslice

!===============================================================================
!
! This program calculates the average potential ona set of planes and can be us-
! ed to determine ionization potential and electron affinity.
!
! This program requires an input file called potavg.inp which
! must be in the form below.
!   pot.cube    # name of cube file
!   6           # number of planes
!   1 0.0       # dir of plane, 1 for a, 2 for b and 3 for c
!   1 1.0       # the following float stands for fractional coordinate
!   2 0.0
!   2 1.0
!   3 0.0
!   3 1.0
!
!===============================================================================


program potavg

    implicit none
    
    !
    ! vars to store information from potavg.in
    !
    character(len=32) :: cubename
    integer(kind=4)   :: nplane
    integer(kind=4), allocatable, dimension(:) :: plane_sym, plane_ind
    real(kind=8),    allocatable, dimension(:) :: plane_pos
    
    !
    ! vars to store cube file information
    !
    real(kind=8),    dimension(3)   :: mesh_origin
    integer(kind=4), dimension(3)   :: mesh_size
    real(kind=8),    dimension(3,3) :: mesh_unit
    real(kind=8), allocatable, dimension(:,:,:) :: field
    
    !
    ! vars to calculate average potential
    !
    real(kind=8), dimension(3) :: dS, S
    real(kind=8) :: integi, total_integ, total_S, avg_pot
    
    !
    ! loop counters
    !
    integer(kind=4) :: i
    
    !
    ! misc
    !
    integer(kind=4) :: tmpi
    real(kind=8) :: si
    
    !
    ! read potavg.in
    !
    open ( unit=8, file="potavg.inp", status="old" )
    read (8,*) cubename
    read (8,*) nplane
    allocate( plane_sym(nplane) )
    allocate( plane_ind(nplane) )
    allocate( plane_pos(nplane) )
    do i = 1, nplane
        read (8,*) plane_sym(i), plane_pos(i)
    end do
    
    !
    ! read cube
    !
    call cube_read_meshinfo( cubename, mesh_origin, mesh_size, mesh_unit )
    allocate( field(mesh_size(1),mesh_size(2),mesh_size(3) ) )
    call cube_read_scafield( cubename, mesh_size, field )
    
    !
    ! plane_pos -> plane_ind
    !
    do i = 1, nplane
        tmpi = nint( plane_pos(i) * mesh_size(plane_sym(i)) )
        if ( tmpi < 1 ) then
            tmpi = 1
        else
            if ( tmpi > mesh_size(plane_sym(i)) ) then
                tmpi = mesh_size(plane_sym(i))
            end if
        end if
        plane_ind(i) = tmpi
    end do
    
    !
    ! calculate dS and S
    !
    call calc_dS( mesh_unit(2,:), mesh_unit(3,:), dS(1) )
    call calc_dS( mesh_unit(3,:), mesh_unit(1,:), dS(2) )
    call calc_dS( mesh_unit(1,:), mesh_unit(2,:), dS(3) )
    S(1) = dS(1) * (mesh_size(2)-1) * (mesh_size(3)-1)
    S(2) = dS(2) * (mesh_size(3)-1) * (mesh_size(1)-1)
    S(3) = dS(3) * (mesh_size(1)-1) * (mesh_size(2)-1)
    
    !
    ! calculate integral
    !
    total_integ = 0.0
    total_S     = 0.0
    do i = 1, nplane
        call calc_integral( mesh_size, field, plane_sym(i), plane_ind(i), &
                            dS(plane_sym(i)), integi )
        total_integ = total_integ + integi
        total_S     = total_S     + S(plane_sym(i))
    end do
    avg_pot = total_integ / total_S
    
    !
    ! output
    !
    write (*,"(A20,3I4)") "Grid dimensions :", mesh_size(:)
    write (*,"(A20,I4)")  "Number of planes :", nplane
    do i = 1, nplane
        select case (plane_sym(i))
        case (1)
            write (*,"(A20,3I4)") "Plane indice :", plane_ind(i), 0, 0
        case (2)
            write (*,"(A20,3I4)") "Plane indice :", 0, plane_ind(i), 0
        case (3)
            write (*,"(A20,3I4)") "Plane indice :", 0, 0, plane_ind(i)
        case default
            write (*,"(A20,3I4)") "Plane indice :", 0, 0, 0
        end select
    end do
    write (*,"(A20,ES12.5)") "Total integral :", total_integ
    write (*,"(A20,ES12.5)") "Total area :", total_S
    write (*,"(A20,ES12.5)") "Average potential :", avg_pot
    
    deallocate( plane_sym )
    deallocate( plane_ind )
    deallocate( plane_pos )
    deallocate(     field )
    
end program potavg

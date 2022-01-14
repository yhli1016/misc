!-------------------------------------------------------------------------------
!
! Sometimes GW and BSE calculations are performed on different hosts. In this
! case, two WFN_co files will be created. However, small numerical errors are
! unavoidable due to the difference in hardware and software configurations of
! different hosts, which will cause inteqp and absorption codes refuse to work.
!
! This program solves the problem by replacing the meanfield energies in eqp1.dat
! obtained on another host with the accurate energies in 02-wfn/xxx.save/K0000X/
! eigenval.xml obtained on the current host.
!
! This program requires one input file called repemf.inp which must contain the
! following lines.
! 
! emf.dat        # file name of meanfield energies produced by getemf.sh
! 45 40          # number of kpoints and bands of meanfield energies
! eqp1.dat       # file name of eqp1.dat
! 45 18          # number of kppints and bands of eqp1.dat
! out.dat        # filename of output file
!
!-------------------------------------------------------------------------------

program repemf

        implicit none

        ! repemf.inp
        character(len=32) :: emfname, eqpname, outname
        integer(kind=4) :: nkemf, nbemf, nkeqp, nbeqp
        ! emf.dat
        real(kind=8), allocatable :: emf(:,:)
        ! eqp.dat
        ! emf2 is read but never used
        real(kind=8), allocatable :: kpoints(:,:), emf2(:,:), eqp(:,:)
        integer(kind=4), allocatable :: spin(:,:), bandind(:,:)
        ! loop counters
        integer(kind=4) :: ik, ib
        ! temporary variables
        integer(kind=4) :: foo

        ! open repemf.inp
        open (unit=7, file="repemf.inp", status="old")
        read (7,*) emfname
        read (7,*) nkemf, nbemf
        read (7,*) eqpname
        read (7,*) nkeqp, nbeqp
        read (7,*) outname
        close (unit=7)

        ! memory allocatation
        allocate(emf(nbemf, nkemf))
        allocate(kpoints(3, nkeqp))
        allocate(emf2(nbeqp, nkeqp))
        allocate(eqp(nbeqp, nkeqp))
        allocate(spin(nbeqp, nkeqp))
        allocate(bandind(nbeqp, nkeqp))

        ! read emf.dat
        open (unit=8, file=emfname, status="old")
        do ik = 1, nkemf
                do ib = 1, nbemf
                        read (8,*) emf(ib,ik)
                end do
        end do
        close (unit=8)
        ! emf is in Hartree so we have to transform it to eV
        emf = emf * 27.21138506
        
        ! read eqp.dat
        ! emf2 is read but never used
        open (unit=9, file=eqpname, status="old")
        do ik = 1, nkeqp
                read (9,*) kpoints(:,ik), foo
                do ib = 1, nbeqp
                        read (9,*) spin(ib,ik), bandind(ib,ik), emf2(ib,ik), eqp(ib,ik)
                end do
        end do
        close (unit=9)

        ! write to out.dat
        open (unit=10, file=outname, status="replace")
        do ik = 1, nkeqp
                write (10,"(3F13.9,I8)") kpoints(:,ik), nbeqp
                do ib = 1, nbeqp
                        write (10,"(2I8,2F15.9)") spin(ib,ik), bandind(ib,ik), emf(bandind(ib,ik),ik), eqp(ib,ik)
                end do
        end do
        close (unit=10)

        ! free memory
        deallocate(emf)
        deallocate(kpoints)
        deallocate(emf2)
        deallocate(eqp)
        deallocate(spin)
        deallocate(bandind)

end program repemf

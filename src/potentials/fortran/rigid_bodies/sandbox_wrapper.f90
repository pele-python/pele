module sandbox_wrapper
    use sandbox_commons
    use sandbox_module

    implicit none

    contains

    subroutine setup_commons(nmol, nsites, r, s, os, as)
    use sandbox_commons
    implicit none

    integer, intent(in) :: nmol, nsites
    double precision, intent(in) :: r !radius
    double precision, intent(in) :: s, os, as !step sizes
    integer, parameter :: npar = 1

    write(*,*) "setting up sandbox_commons"

    natoms = nmol*2 !I think

    allocate(step(npar))
    allocate(ostep(npar))
    allocate(astep(npar))

    allocate(tmove(npar))
    allocate(omove(npar))

    allocate( coords(3*natoms, npar) )
    allocate( vt(nsites) )

    myunit = 1 !should be set later
    !natoms !already set
    nsave = 1 !should be set later
    mynode = 0

    perccut = 1.d0 !only used with percolatet used
    radius = R !size of system should be large

    percolatet = .false.
    twod = .false.
    mpit = .false. !should be set later

    step(1)  = s
    ostep(1) = os
    astep(1) = as
    !vt(:) !will be set later
    !coords(:,:) !will be set later

    tmove(1) = .true. !will be set later
    omove(1) = .false. !will be set later

    end subroutine setup_commons

    subroutine input(fname)
    implicit none
    character(len=200) :: fname
    integer myun
    myun = 77
    write(*,*) "reading ", fname
    OPEN(UNIT=myun,FILE=fname,STATUS='UNKNOWN')
    call sandbox_input(myun)
    close(myun)
    end subroutine input


    subroutine getEnergyGradient(natoms, mycoords, mygrad, energy)
    !use sandbox_commons, only : natoms
    implicit none
    integer, intent(in) :: natoms
    double precision, intent(in) :: mycoords(3*natoms)
    double precision, intent(out) :: mygrad(3*natoms), energy
    double precision mynewcoords(3*natoms, 1)
    mynewcoords(:,1) = mycoords(:)
    call sandbox( mycoords, mygrad, energy, .true. )
    end subroutine getEnergyGradient

end module sandbox_wrapper


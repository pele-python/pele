module sandbox_commons
    implicit none

    integer          :: myunit, natoms, nsave, mynode
    double precision :: perccut, radius
    double precision, allocatable :: step(:), ostep(:), astep(:), vt(:), coords(:,:)
    logical          :: percolatet, twod, mpit
    logical, allocatable     :: tmove(:), omove(:)


    contains

    subroutine setup(nmol, nsites, r, s, os, as)
    implicit none

    integer, intent(in) :: nmol, nsites
    double precision, intent(in) :: r !radius
    double precision, intent(in) :: s, os, as !step sizes
    integer, parameter :: npar = 1

    natoms = nmol*2 !I think

    allocate(step(npar))
    allocate(ostep(npar))
    allocate(astep(npar))

    allocate(tmove(npar))
    allocate(omove(npar))

    allocate( coords(natoms, npar) )
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

    end subroutine setup

end module sandbox_commons

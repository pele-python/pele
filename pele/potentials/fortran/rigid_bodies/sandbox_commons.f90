module sandbox_commons
    implicit none

    integer          :: myunit, natoms, nsave, mynode
    double precision :: perccut, radius
    double precision, allocatable :: step(:), ostep(:), astep(:), vt(:), coords(:,:)
    logical          :: percolatet, twod, mpit
    logical, allocatable     :: tmove(:), omove(:)


    contains


end module sandbox_commons

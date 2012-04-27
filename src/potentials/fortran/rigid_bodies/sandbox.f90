module sandbox_module
    implicit none

    type sandbox_site
        type(sandbox_molecule), pointer :: molecule
        integer :: class    ! interaction class
        double precision :: r_molframe(3), p_molframe(3), rotmat_molframe(3, 3) ! molecule frame position and orientation
        double precision :: r(3), p(3), rotmat(3, 3)  ! lab frame position and orientation
        double precision :: rotmat_deriv(3, 3, 3)   ! dR^i/dp^I = deriv with resp to molecule angle-axis

        logical :: is_aniso

        logical :: has_pole
        double precision :: mu_hat0(3), mu_hat(3), mu_hat_deriv(3, 3)   ! mu_hat derivatives with resp to p^i
    end type sandbox_site

    type sandbox_molecule
        type(sandbox_site), allocatable :: sites(:)
        double precision :: r(3), p(3)
        double precision :: rotmat(3, 3), rotmat_deriv(3, 3, 3)
        integer :: r_index, p_index
    end type sandbox_molecule
    
    type interaction
        logical :: exists
        integer :: func_id  ! see table below for potential <-> func_id correspondence
        double precision, allocatable :: params(:)
    end type interaction

    type interaction_class
        logical :: has_pole, is_aniso
    end type interaction_class

    double precision, parameter :: z_hat(3) = (/ 0.D0, 0.D0, 1.D0 /)

    type(sandbox_molecule), allocatable, target :: molecules(:)
    type(interaction), allocatable, target :: interactions(:, :)
    type(interaction_class), allocatable :: classes(:)

    integer :: num_mols, x_size   ! length of x array

    ! to add a new potential, add a sandbox_[potential name] subroutine, update pairwise_sandbox, and update
    ! sandbox_input

    ! table of potentials and their integer IDs func_id:
    ! lj 1
    ! coulomb 2
    ! dipole 3
    ! chiro 4

contains
    subroutine initialize_sandbox_molecule(mol, mol_index)
        implicit none
        type(sandbox_molecule), target, intent(inout) :: mol
        integer, intent(in) :: mol_index

        integer :: site_index
        type(sandbox_site), pointer :: site
        double precision :: dummy(3, 3)

        ! initialize molecule values
        mol%r_index = 3 * mol_index - 2
        mol%p_index = mol%r_index + 3 * num_mols

        do site_index = 1, size(mol%sites)
            site => mol%sites(site_index)
            site%molecule => mol    ! give all sites a pointer to their parent

            ! if the site is anisotropic, need to do some extra bookkeeping
            if(classes(site%class)%is_aniso) then
                site%is_aniso = .true.

                ! compute rotation matrix for site
                call rmdrvt(site%p_molframe, site%rotmat_molframe, dummy, dummy, dummy, .false.)

                ! if site has a pole, turn on a marker to keep track of mu_hat
                if(classes(site%class)%has_pole) then
                    site%has_pole = .true.
                    site%mu_hat0 = matmul(site%rotmat_molframe, z_hat)
                else
                    site%has_pole = .false.
                end if
            end if
        end do
    end subroutine initialize_sandbox_molecule

    subroutine update_sandbox_molecule(gtest, mol, r, p)
        implicit none
        logical, intent(in) :: gtest
        type(sandbox_molecule), target, intent(inout) :: mol
        double precision, intent(inout) :: r(3), p(3)

        double precision, parameter :: pi = 3.14159265359
        double precision :: pmod
        type(sandbox_site), pointer :: site
        integer :: site_index, k
        double precision :: dummy(3, 3), dR_dp(3, 3, 3)

        double precision :: swo(4, 3, 3)

        mol%r(:) = r(:)

        ! make sure that 0 < |p| < 2*pi
        pmod = sqrt(dot_product(p, p))
        if(pmod > 2 * pi) p(:) = p(:) / pmod * mod(pmod, 2 * pi)
        mol%p(:) = p(:)

        ! recompute the rotation matrices and derivatives
        call rmdrvt(mol%p, mol%rotmat, &
            & mol%rotmat_deriv(1, :, :), mol%rotmat_deriv(2, :, :), mol%rotmat_deriv(3, :, :), gtest)

        ! loop over sites
        do site_index = 1, size(mol%sites)
            site => mol%sites(site_index)

            site%r(:) = mol%r(:) + matmul(mol%rotmat, site%r_molframe(:))   ! reset site position

            if(site%is_aniso) then
                ! if site has a pole, update the pole
                if(site%has_pole) then
                    ! mu_hat = R^I . mu_hat0 = R^I . R^i0 . z_hat
                    site%mu_hat = matmul(mol%rotmat, site%mu_hat0)

                    ! dmu/dp^I_k = dR^I/dp^I_k . mu_hat0
                    do k = 1, 3
                        site%mu_hat_deriv(k, :) = matmul(mol%rotmat_deriv(k, :, :), site%mu_hat0)
                    end do
                end if
            end if
        end do

    end subroutine update_sandbox_molecule

    subroutine pairwise_sandbox(gtest, sitei, sitej, energy_contrib, grad_contrib)
        use rotations
        use vec3
        implicit none
        logical, intent(in) :: gtest
        type(sandbox_site), intent(in) :: sitei, sitej
        double precision, intent(out) :: energy_contrib, grad_contrib(12)

        integer :: classi, classj
        type(interaction), pointer :: this_interaction
        integer, pointer :: func_id
        double precision, pointer :: params(:)
        double precision :: dr_dp(3, 3)
        integer :: k

        energy_contrib = 0.D0
        grad_contrib(:) = 0.D0

        classi = min(sitei%class, sitej%class)
        classj = max(sitei%class, sitej%class)

        if(interactions(classi, classj)%exists) then

        func_id => interactions(classi, classj)%func_id
        params => interactions(classi, classj)%params

            if(func_id == 1) then
                ! lj
                ! sigma_0 eps_0 rep, att
                call sandbox_lj(gtest, params(1), params(2), params(3), params(4), sitei%r - sitej%r, energy_contrib, grad_contrib)
            elseif(func_id == 2) then
                ! coulomb
                ! sigma_0 k_c*q*q
                call sandbox_coulomb(gtest, params(1), params(2), sitei%r - sitej%r, energy_contrib, grad_contrib)
            elseif(func_id == 3) then
                ! dipole
                ! sigma_0 mu^2
                call sandbox_dipole(gtest, params(1), params(2), sitei%r - sitej%r, sitei%mu_hat, &
                    & sitej%mu_hat, sitei%mu_hat_deriv, sitej%mu_hat_deriv, energy_contrib, grad_contrib)
            elseif(func_id == 4) then
                ! chiro
                ! sigma_0 mu^2 cos_gamma sin_gamma
                call sandbox_chiro(gtest, params(1), params(2), params(3), params(4), sitei%r - sitej%r, sitei%mu_hat, &
                    & sitej%mu_hat, sitei%mu_hat_deriv, sitej%mu_hat_deriv, energy_contrib, grad_contrib)
            else
                print *, "pairwise_sandbox> don't have potential with func_id='", func_id, &
                    & "' for interaction of classes", classi, classj
                return
            end if

            if(gtest) then
                ! if site is not at molecule origin, need to compute dr^i/dp^I
                if(any(sitei%r_molframe /= 0.D0)) then
                    do k = 1, 3
                        ! dr^i/dp^I_k = R^I_k . r^i0
                        dr_dp(k, :) = matmul(sitei%molecule%rotmat_deriv(k, :, :), sitei%r_molframe)
                    end do

                    ! assign dU/dp^I = dr^i/dp^I . dU/dr^i
                    grad_contrib(7: 9) = grad_contrib(7: 9) + matmul(dr_dp, grad_contrib(1: 3))
                end if

                ! repeat for molecule J
                if(any(sitej%r_molframe /= 0.D0)) then
                    do k = 1, 3
                        dr_dp(k, :) = matmul(sitej%molecule%rotmat_deriv(k, :, :), sitej%r_molframe)
                    end do

                    grad_contrib(10: 12) = grad_contrib(10: 12) + matmul(dr_dp, grad_contrib(4: 6))
                end if
            end if
        end if

    end subroutine pairwise_sandbox


! --- begin sandbox potentials
    subroutine sandbox_lj(gtest, sigma_0, eps_0, rep, att, rij, energy, grad)
        ! gtest: whether gradients should be calculated
        ! sigma_0: length scale
        ! eps_0: energy_scale
        ! rep: repulsive coefficient, i.e., A in A / r^12 (for normal LJ, this is 1.0)
        ! att: attractive coefficient, i.e., B in B / r^6 (for normal LJ, this is 1.0)
        ! rij: r_j - r_i

        logical, intent(in) :: gtest
        double precision, intent(in) :: sigma_0, eps_0, rep, att, rij(3)
        double precision, intent(out) :: energy, grad(12)

        double precision :: rijmod, rij_hat(3) ! |rij| and normal vector pointing rij

        ! r = (|rij| / sigma_0), r_pm1 = r to power minus 1 = r^-1, etc
        double precision :: r_pm1, r_pm6, r_pm7, r_pm12, r_pm13

        energy = 0.D0
        if(gtest) grad(:) = 0.D0

        rijmod = sqrt(dot_product(rij, rij))

        r_pm1 = sigma_0 / rijmod
        r_pm6 = r_pm1 * r_pm1 * r_pm1 * r_pm1 * r_pm1 * r_pm1 
        r_pm12 = r_pm6 * r_pm6

        energy = 4.0 * eps_0 * (rep * r_pm12 - att * r_pm6)

        if(gtest) then
            rij_hat(:) = rij / rijmod
            r_pm7 = r_pm6 * r_pm1
            r_pm13 = r_pm12 * r_pm1

            grad(1: 3) = 24.0 * eps_0 / sigma_0 * (-2.0 * rep * r_pm13 + att * r_pm7) * rij_hat
            grad(4: 6) = -grad(1: 3)
            grad(7: 9) = 0.D0
            grad(10: 12) = 0.D0
        end if

    end subroutine sandbox_lj

    subroutine sandbox_coulomb(gtest, sigma_0, kqq, rij, energy, grad)
        ! gtest: if gradients are computed
        ! sigma_0: length scale for rij
        ! kqq: Coulomb's constant * charge_1 * charge_2
        ! rij: separation between sites

        ! energy: pairwise energy
        ! gradient: (1:3) = radial for site i, (4:6) = radial j, (7:9) = angular for i, (10:12), angular j

        logical, intent(in) :: gtest
        double precision, intent(in) :: sigma_0, kqq, rij(3)
        double precision, intent(out) :: energy, grad(12)

        double precision :: rijmod, rij_hat(3) ! |rij| and normal vector pointing rij

        ! r = (|rij| / sigma_0), r_pm1 = r to power minus 1 = r^-1, etc
        double precision :: r_pm1, r_pm2

        energy = 0.D0
        if(gtest) grad(:) = 0.D0

        rijmod = sqrt(dot_product(rij, rij))

        r_pm1 = sigma_0 / rijmod

        energy = kqq * r_pm1

        if(gtest) then
            rij_hat(:) = rij / rijmod
            r_pm2 = r_pm1 * r_pm1

            grad(1: 3) = -kqq / sigma_0 * r_pm2 * rij_hat(:)
            grad(4: 6) = -grad(1: 3)
            grad(7: 9) = 0.D0
            grad(10: 12) = 0.D0
        end if

    end subroutine sandbox_coulomb

    subroutine sandbox_dipole(gtest, sigma_0, mu_p2, rij, mu_hati, mu_hatj, mu_hat_derivi, mu_hat_derivj, energy, grad)
        use vec3, only: vec_cross, vec_dyad, identity3x3
        logical, intent(in) :: gtest
        double precision, intent(in) :: sigma_0, mu_p2, rij(3), mu_hati(3), mu_hatj(3), &
            & mu_hat_derivi(3, 3), mu_hat_derivj(3, 3)
        double precision, intent(out) :: energy, grad(12)

        double precision :: rijmod, rij_hat(3)
        double precision :: r_pm1, r_pm2, r_pm3, r_pm6
        double precision :: mu_dot_ri, mu_dot_rj, mu_dot_mu
        integer :: k

        energy = 0.D0
        if(gtest) grad(:) = 0.D0

        rijmod = sqrt(dot_product(rij, rij))
        if(rijmod == 0.D0) then
            rij_hat(:) = 0.D0
        else
            rij_hat(:) = rij / rijmod
        end if

        r_pm1 = sigma_0 / rijmod
        r_pm2 = r_pm1 * r_pm1
        r_pm3 = r_pm2 * r_pm1
        r_pm6 = r_pm2 * r_pm2 * r_pm2

        mu_dot_ri = dot_product(mu_hati, rij_hat)
        mu_dot_rj = dot_product(mu_hatj, rij_hat)
        mu_dot_mu = dot_product(mu_hati, mu_hatj)

        energy = mu_p2 * r_pm3 * (mu_dot_mu - 3.0 * mu_dot_ri * mu_dot_rj)

        if(gtest) then
            grad(1: 3) = -3.0 * mu_p2 * r_pm3 / rijmod * &
                & ((mu_dot_mu - 5.0 * mu_dot_ri * mu_dot_rj) * rij_hat(:) + mu_dot_ri * mu_hatj + mu_dot_rj * mu_hati)
            grad(4: 6) = -grad(1: 3)

            do k = 1, 3
                grad(6 + k) = mu_p2 * r_pm3 * (dot_product(mu_hat_derivi(k, :), mu_hatj) &
                    & - 3.0 * dot_product(mu_hat_derivi(k, :), rij_hat) * mu_dot_rj)
                grad(9 + k) = mu_p2 * r_pm3 * (dot_product(mu_hat_derivj(k, :), mu_hati) &
                    & - 3.0 * dot_product(mu_hat_derivj(k, :), rij_hat) * mu_dot_ri)
            end do
        end if

    end subroutine sandbox_dipole

    subroutine sandbox_chiro(gtest, sigma_0, mu_p2, cos_gamma, sin_gamma, rij, mu_hati, mu_hatj, &
      & mu_hat_derivi, mu_hat_derivj, energy, grad)
        use vec3, only: vec_cross
        implicit none
        logical, intent(in) :: gtest
        double precision, intent(in) :: sigma_0, mu_p2, sin_gamma, cos_gamma, rij(3), mu_hati(3), mu_hatj(3), &
            & mu_hat_derivi(3, 3), mu_hat_derivj(3, 3)
        double precision, intent(out) :: energy, grad(12)

        double precision :: rijmod, rij_hat(3)
        double precision :: r_pm1, r_pm3, r_pm4 ! r is in reduced units
        double precision :: mu_dot, mu_cross(3)  ! mu_hati . mu_hatj and mu_hati x mu_hatj

        integer :: k

        energy = 0.D0
        if(gtest) grad(:) = 0.D0

        rijmod = sqrt(dot_product(rij, rij))
        if(rijmod == 0.D0) then
            rij_hat(:) = 0.D0
        else
            rij_hat(:) = rij / rijmod
        end if
        
        r_pm1 = sigma_0 / rijmod
        r_pm3 = r_pm1 * r_pm1 * r_pm1

        mu_dot = dot_product(mu_hati, mu_hatj)
        mu_cross(:) = vec_cross(mu_hati, mu_hatj)

        energy = -mu_p2 * r_pm3 * (cos_gamma * mu_dot + sin_gamma * dot_product(mu_cross, rij_hat))

        if(gtest) then
            r_pm4 = r_pm3 * r_pm1

            grad(1: 3) = -mu_p2 * r_pm4 / sigma_0 * (cos_gamma * (-3.0 * mu_dot * rij_hat(:)) + &
                & sin_gamma * (mu_cross(:) - 4.0 * dot_product(mu_cross, rij_hat) * rij_hat(:)))
            grad(4: 6) = -grad(1: 3)

            do k = 1, 3
                grad(6 + k) = -mu_p2 * r_pm3 * (cos_gamma * dot_product(mu_hat_derivi(k, :), mu_hatj) & 
                    & + sin_gamma * dot_product(vec_cross(mu_hat_derivi(k, :), mu_hatj), rij_hat))

                grad(9 + k) = -mu_p2 * r_pm3 * (cos_gamma * dot_product(mu_hati, mu_hat_derivj(k, :)) & 
                    & + sin_gamma * dot_product(vec_cross(mu_hati, mu_hat_derivj(k, :)), rij_hat))
            end do
        end if

    end subroutine sandbox_chiro
! --- end sandbox potentials


! --- begin sandbox utilities
    function lower_case(word)
        character(len=*), intent(in) :: word
        character(len=len(word)) :: lower_case

        character(len=26), parameter :: uc = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        character(len=26), parameter :: lc = 'abcdefghijklmnopqrstuvwxyz'

        integer :: i, ind

        do i = 1, len(word)
            ind = index(uc, word(i: i))
            if (ind /= 0) then
                lower_case(i: i) = lc(ind: ind)
            else
                lower_case(i: i) = word(i: i)
            end if
        end do

    end function lower_case

    function first_char(u) result(c)
        implicit none
        integer, intent(in) :: u
        character :: c
        integer :: stat

        do
            read(u, '(a)', advance='no', iostat=stat) c
            if(stat /= 0) then
                c = ''
                backspace(u)
                return
            elseif(c == ' ') then
                cycle
            else
                backspace(u)
                return
            end if
        end do
    end function first_char

    subroutine print_matrix(m)
        implicit none
        double precision, intent(in) :: m(3, 3)

        integer :: i
        
        do i = 1, 3
            print *, m(i, 1), m(i, 2), m(i, 3)
        end do

    end subroutine print_matrix

    function max_distance(mols) result(d)
        ! compute the maximum distance of a molecule from the origin

        use vec3, only: vec_len

        implicit none
        type(sandbox_molecule), intent(in) :: mols(:)
        double precision :: d

        double precision :: this_d
        integer :: i

        d = 0.D0

        do i = 1, size(mols)
            this_d = vec_len(mols(i)%r)
            if(this_d > d) d = this_d
        end do

    end function max_distance
! --- end sandbox utilities

end module sandbox_module

subroutine sandbox_input(out_unit)
    use sandbox_commons, only: natoms
    use sandbox_module
    implicit none    

    integer, intent(in) :: out_unit ! probably GMIN_out

    integer :: sand_unit, stat
    type(sandbox_molecule), pointer :: mol
    type(sandbox_site), pointer :: site
    integer :: arity, num_these_mols, num_sites, num_classes
    character(len=20) :: nary_check, label(3)
    double precision :: number_fraction
    integer :: getunit, mol_index, mol_type_index, site_index
    character(len=20) :: buffer

    integer :: classi, classj, tmp
    character(len=20) :: potential_name
    type(interaction), pointer :: this_int
    character(len=20), pointer :: func
    double precision, pointer :: params(:)

    character :: c
    integer :: i, j

    x_size = 3 * natoms
    num_mols = natoms / 2
    allocate(molecules(num_mols))

    num_classes = 0

    sand_unit = getunit()
    open(unit=sand_unit, file='rbdata', status='old')

    ! determine if first line is an integer num_sites or a string 'n-ary'
    read(sand_unit, *) buffer
    read(buffer, *, iostat=stat) num_sites
    if(stat == 0) then
        ! it was an integer, so do single site
        arity = 1
        num_these_mols = num_mols
    else
        ! it was hopefully a string
        read(buffer, *, iostat=stat) nary_check

        ! sanity check: should have copied a string 'n-ary'
        if(.not.(stat == 0 .and. nary_check == 'n-ary')) then
            print *, 'woah, son: got to put in some proper input for rb!'
            return
        end if

        ! it was a string 'n-ary', so do n-ary mixture
        read(sand_unit, *) arity
        write(out_unit, *) 'sandbox_input> arity of mixture:', arity
    end if

    mol_index = 1
    do mol_type_index = 1, arity
        read(sand_unit, *)  ! read blank line
        if(arity == 1) then
            write(out_unit, *) 'sandbox_input>', num_these_mols, 'rigid bodies'
        elseif(arity > 1) then
            read(sand_unit, *) number_fraction  ! number fraction of this molecule type
            read(sand_unit, *) num_sites    ! number of sites in this molecule

            num_these_mols = nint(num_mols * number_fraction)
            write(out_unit, *) 'adding', num_these_mols, 'molecules of type', mol_type_index
        end if

        if(num_these_mols == 0) then    ! the number fraction is too small to give even one molecule
            read(sand_unit, *)  ! just skip the line
        else
            allocate(molecules(mol_index)%sites(num_sites)) ! prepare space for the sites in the molecule

            site_index = 1
            do
                if(site_index > num_sites) exit

                ! skip blank lines and lines beginning with #
                c = first_char(sand_unit)
                if(c == '#' .or. c == '') then
                    read(sand_unit, *, iostat=stat)
                    if(stat /= 0) exit
                    cycle
                end if

                ! if the line wasn't blank or comment, start assigning values
                site => molecules(mol_index)%sites(site_index)
                read(sand_unit, *) label(1), site%class, label(2), site%r_molframe(:), label(3), site%p_molframe(:)

                ! update number of classes
                if(site%class > num_classes) num_classes = site%class

                ! increment the site index
                site_index = site_index + 1
            end do
            mol_index = mol_index + 1

            ! copy the first molecule into the rest of the mols of this type
            do mol_index = mol_index, mol_index + (num_these_mols - 2)
                allocate(molecules(mol_index)%sites(num_sites))
                molecules(mol_index) = molecules(mol_index - 1)
            end do
        end if
    end do

    ! read in blank lines until reach 'matrix'
    do
        read(sand_unit, *, iostat=stat) buffer
        if(stat /= 0) then
            print *, "sandbox> where's the interation matrix?"
            return
        end if
        if(buffer == 'matrix') exit
    end do

    ! assign the interaction matrix
    allocate(classes(num_classes))
    classes%is_aniso = .false.
    classes%has_pole = .false.
    allocate(interactions(num_classes, num_classes))
    interactions%exists = .false.
    do
        ! skip blank lines and lines beginning with #
        c = first_char(sand_unit)
        if(c == '#' .or. c == '') then
            read(sand_unit, *, iostat=stat)
            if(stat /= 0) then
                exit
            end if
            cycle
        end if

        read(sand_unit, *, iostat=stat) classi, classj, potential_name
        if(stat /= 0) exit

        if(classi > num_classes .or. classj > num_classes) then
            print *, 'sandbox_input> youve specified interactions for classes with no sites'
            exit
        end if

        backspace(sand_unit)    ! get ready to read same line again

        ! keep interaction information in upper-right matrix
        if(classi > classj) then
            tmp = classi
            classi = classj
            classj = tmp
        end if

        ! assign existence
        interactions(classi, classj)%exists = .true.
        !call locase(potential_name) !js850> this is supposes to make the string lower case

        ! read in input params and do some internal processing if necessary
        if(potential_name == 'lj') then
            ! lj sigma_0 eps_0 rep att
            interactions(classi, classj)%func_id = 1
            allocate(interactions(classi, classj)%params(4))
            read(sand_unit, *, iostat=stat) tmp, tmp, buffer, interactions(classi, classj)%params(1: 4)
        else if(potential_name == 'coulomb') then
            ! coulomb sigma_0 k_c q q -> sigma_0 kcqq
            interactions(classi, classj)%func_id = 2
            allocate(interactions(classi, classj)%params(4))
            read(sand_unit, *, iostat=stat) tmp, tmp, buffer, interactions(classi, classj)%params(1: 4)
            interactions(classi, classj)%params(2) = interactions(classi, classj)%params(2) * &
                & interactions(classi, classj)%params(3) * interactions(classi, classj)%params(4)
        else if(potential_name == 'dipole') then
            ! dipole sigma_0 mu -> sigma_0 mu^2
            interactions(classi, classj)%func_id = 3
            allocate(interactions(classi, classj)%params(2))
            read(sand_unit, *, iostat=stat) tmp, tmp, buffer, interactions(classi, classj)%params(1: 2)
            interactions(classi, classj)%params(2) = interactions(classi, classj)%params(2) ** 2
            classes(classi)%is_aniso = .true.
            classes(classj)%is_aniso = .true.
            classes(classi)%has_pole = .true.
            classes(classj)%has_pole = .true.
        else if(potential_name == 'chiro') then
            ! chiro sigma_0 mu gamma -> sigma_0 mu^2 cos_gamma sin_gamma
            interactions(classi, classj)%func_id = 4
            allocate(interactions(classi, classj)%params(4))
            read(sand_unit, *, iostat=stat) tmp, tmp, buffer, interactions(classi, classj)%params(1: 3)
            interactions(classi, classj)%params(2) = interactions(classi, classj)%params(2) ** 2
            interactions(classi, classj)%params(4) = sin(interactions(classi, classj)%params(3))
            interactions(classi, classj)%params(3) = cos(interactions(classi, classj)%params(3))
            classes(classi)%is_aniso = .true.
            classes(classj)%is_aniso = .true.
            classes(classi)%has_pole = .true.
            classes(classj)%has_pole = .true.
        else
            print *, 'sandbox_input> entered unknown interaction name!'
            exit
        end if
    end do

    ! initialize molecules
    do mol_index = 1, size(molecules)
        call initialize_sandbox_molecule(molecules(mol_index), mol_index)
    end do

end subroutine sandbox_input

subroutine sandbox(x, grad, energy, gtest)
    use sandbox_commons, only: vt, twod
    use sandbox_module
    implicit none
    double precision, intent(inout) :: x(x_size)
    logical, intent(in) :: gtest
    double precision, intent(out) :: energy, grad(x_size)

    type(sandbox_molecule), pointer :: moli, molj
    type(sandbox_site), pointer :: sitei, sitej
    integer moli_index, molj_index, sitei_index, sitej_index
    double precision :: energy_contrib, grad_contrib(12)

    energy = 0.D0
    vt(:) = 0.D0
    if(gtest) grad(:) = 0.D0

    ! update the values for all the molecules
    do moli_index = 1, num_mols
        moli => molecules(moli_index)
        call update_sandbox_molecule(gtest, moli, x(moli%r_index: moli%r_index + 2), &
            & x(moli%p_index: moli%p_index + 2))
    end do

    ! outer loop over molecules
    do moli_index = 1, num_mols - 1
        moli => molecules(moli_index)

        ! inner loop over molecules
        do molj_index = moli_index + 1, num_mols
            molj => molecules(molj_index)

            ! loop over sites in outer molecule
            do sitei_index = 1, size(moli%sites)
                sitei => moli%sites(sitei_index)

                ! loop over sites in inner molecule
                do sitej_index = 1, size(molj%sites)
                    sitej => molj%sites(sitej_index)

                    energy_contrib = 0.D0
                    grad_contrib(:) = 0.D0

                    ! compute the energy and gradient for this pair of sites
                    call pairwise_sandbox(gtest, sitei, sitej, energy_contrib, grad_contrib)

                    ! add the energy to total and pairwise
                    energy = energy + energy_contrib
                    vt(moli_index) = vt(moli_index) + energy_contrib
                    vt(molj_index) = vt(molj_index) + energy_contrib

                    if(gtest) then  ! add gradient contributions
                        grad(moli%r_index: moli%r_index + 2) = grad(moli%r_index: moli%r_index + 2) + grad_contrib(1: 3)
                        grad(molj%r_index: molj%r_index + 2) = grad(molj%r_index: molj%r_index + 2) + grad_contrib(4: 6)
                        grad(moli%p_index: moli%p_index + 2) = grad(moli%p_index: moli%p_index + 2) + grad_contrib(7: 9)
                        grad(molj%p_index: molj%p_index + 2) = grad(molj%p_index: molj%p_index + 2) + grad_contrib(10: 12)
                    end if

                end do
            end do
        end do
    end do

    ! 2D: set all z gradients to zero
    if(twod) then
        do moli_index = 1, num_mols
            grad(molecules(moli_index)%r_index + 2) = 0.D0
        end do
    end if

end subroutine sandbox

subroutine sandbox_output
    use sandbox_commons, only: natoms, nsave, mpit, mynode
    use qmodule, only: qmin, qminp, ff
    use sandbox_module
    implicit none

    integer :: sandout_unit, coords_unit
    integer :: min_num, mol_num, site_num, atom_num, atom_index
    type(sandbox_molecule), pointer :: mol
    type(sandbox_site), pointer :: site
    character(len=25) :: sandout_name, coords_name, min_name, node
    double precision :: rotmat(3, 3), dummy(3, 3)

    integer :: getunit

    sandout_unit = getunit()

    ! open sandout file for writing
    if(mpit) then
        write(node, *) mynode + 1
        sandout_name = 'sandout.' // trim(adjustl(node)) // '.xyz'
    else
        sandout_name = 'sandout.xyz'
    end if
    open(unit=sandout_unit, file=sandout_name, status='unknown')

    ! loop over saved minima
    do min_num = 1, nsave
        ! put number of atoms and comment line. for now this we'll do just one atom per molecule
        write(sandout_unit, *) num_mols
        write(sandout_unit, *) 'energy of minimum', min_num, '=', qmin(min_num), 'first found at step', ff(min_num)

        ! loop over molecules
        do mol_num = 1, size(molecules)
            mol => molecules(mol_num)

            call rmdrvt(qminp(min_num, mol%p_index: mol%p_index + 2), rotmat, dummy, dummy, dummy, .false.)

            write(sandout_unit, *) 'O', qminp(min_num, mol%r_index: mol%r_index + 2), &
                & 'atom_vector', matmul(rotmat, z_hat)
        end do

        ! output coords. files
        coords_unit = getunit()
        write(min_name, '(i3)') min_num

        if(mpit) then
            coords_name = 'coords.' // trim(adjustl(min_name)) // '.' // trim(adjustl(node))
        else
            coords_name = 'coords.' // trim(adjustl(min_name))
        end if

        open(coords_unit, file=coords_name, status='unknown')
        do atom_num = 1, natoms
            atom_index = 3 * (atom_num - 1) + 1
            write(coords_unit, *) qminp(min_num, atom_index), qminp(min_num, atom_index + 1), qminp(min_num, atom_index + 2)
        end do
        close(coords_unit)

    end do

    close(sandout_unit)
            
end subroutine sandbox_output

subroutine sandbox_takestep(np)
    use sandbox_commons, only: myunit, coords, radius, percolatet, perccut, tmove, omove, step, ostep, astep, vt, twod
    use sandbox_module
    use vec3, only: vec_len, vec_random
    use rotations, only: rot_takestep_aa

    implicit none
    integer, intent(in) :: np

    double precision, target :: x(size(coords(:, np)))
    double precision, pointer :: r(:), p(:)
    double precision :: step_size, min_vt, dice_roll, twod_step(3)
    type(sandbox_molecule), pointer :: moli, molj
    integer :: moli_index, molj_index
    logical :: stop_pulling ! for angular moves

    x(:) = coords(:, np)
    min_vt = minval(vt)

    ! update ellipsoids before taking any steps
    do moli_index = 1, num_mols
        moli => molecules(moli_index)
        r => x(moli%r_index: moli%r_index + 2)
        p => x(moli%p_index: moli%p_index + 2)

        if(.not. percolatet) then
            if(vec_len(moli%r) > radius) then
                write(myunit, *) 'sandbox_takestep> initial coord outside container, brining in'
                r(:) = r(:) - sqrt(radius) * nint(r(:) / sqrt(radius))
            end if
        end if

        call update_sandbox_molecule(.false., moli, r, p)
    end do

    do moli_index = 1, num_mols
        moli => molecules(moli_index)
        r => x(moli%r_index: moli%r_index + 2)
        p => x(moli%p_index: moli%p_index + 2)

        if(twod) then
            ! translational step
            do
                call random_number(twod_step(1))
                call random_number(twod_step(2))
                twod_step(3) = 0.D0
                if(vec_len(twod_step) <= 1.0) exit
            end do
            r(:) = r(:) + (step(np) * twod_step(:))
        else
            ! angular move: accept step with probability 1 - exp(-beta * delta_energy), where beta = 1/kT = astep
            if(astep(np) > 0.D0) then
                call random_number(dice_roll)
                if(dice_roll > exp(-astep(np) * (vt(moli_index) - min_vt))) then
                    r = max_distance(molecules) * vec_random()

                    ! if using percolate, then pull the molecule in until it is within percolation distance of another
                    if(percolatet) then
                        stop_pulling = .false.
                        do
                            do molj_index = 1, num_mols
                                if(moli_index == molj_index) cycle
                                molj => molecules(molj_index)

                                if(vec_len(moli%r - molj%r) < perccut) then
                                    stop_pulling = .true.
                                    exit
                                end if
                            end do

                            if(stop_pulling) exit
                            r(:) = 0.95 * r(:)
                        end do
                    else ! if using radius, pull in until the molecule is inside the container
                        do
                            if(vec_len(r) < radius) exit
                            r(:) = 0.95 * r(:)
                        end do
                    end if

                end if
            end if

            ! translational move: uniformly distributed in a sphere of radius step(np)
            if(tmove(np)) then
                call random_number(step_size)
                step_size = sqrt(step(np) * step_size)  ! need to take square root to get uniform volume distribution
                r(:) = r(:) + step_size * vec_random()
            end if
        end if ! 2d

        ! orientational move
        if(omove(np)) then
            call random_number(step_size)
            call rot_takestep_aa(p(:), step_size * ostep(np))
        end if

    end do

    ! copy local results back
    coords(:, np) = x(:)

end subroutine sandbox_takestep

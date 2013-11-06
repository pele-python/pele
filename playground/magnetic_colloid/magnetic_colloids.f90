subroutine get_energy_forces(r,f,wpot,num_particles,box,theta)
implicit none
!calculation of forces
integer, intent(in) :: num_particles
real(8), intent(IN), dimension (3) :: box
real(8), intent(in) :: theta
real(8), intent(out) :: wpot
real(8), intent(IN), dimension (num_particles,3) :: r
real(8), intent(out), dimension (num_particles,3) :: f
integer :: i, j, nphi
real(8), dimension (4,num_particles,num_particles) :: rij_array
real(8), dimension (num_particles,3) :: m,mx,my,mz
real(8), dimension (3) :: fij, b
real(8), dimension (3,3) :: brms
real(8), dimension (3) :: rij, unit_rij

real(8) :: u, wpotm
real(8) :: xi, yi, zi
real(8) :: xij, yij, zij
real(8) :: rsq, rsqi, r6, r12
real(8) :: temp_set, d_time
real(8) :: r_cut, wca_cut
real(8) :: m1r, m2r, m12
real(8) :: rij_length
real(8) :: fric, eps

r_cut= 7.d0
eps= 1.d0
wca_cut= 2.d0**(1.d0/6.d0)
u= 0.d0
f= 0.d0
wpot= 0.0d0
forall (i=1:num_particles) mx(i,:)= (/ 1.d0, 0.d0, 0.d0 /)
forall (i=1:num_particles) my(i,:)= (/ 0.d0, 1.d0, 0.d0 /)
forall (i=1:num_particles) mz(i,:)= (/ 0.d0, 0.d0, 1.d0 /)

call Create_rij_array(num_particles,box,r,rij_array)
call generate_brms(theta,brms)

do nphi=1,3
        wpotm= 0.0d0
        b= brms(nphi,:)

        select case (nphi)
        case(1)
        m= mx
        case(2)
        m= my
        case(3)
        m= mz
        end select

        call SelfConsistentEnergy(num_particles,box,rij_array,b,m,wpotm)
        wpot= wpot + wpotm/3.0d0

!$OMP PARALLEL DEFAULT (PRIVATE) &
!$OMP SHARED (rij_array,m,f,num_particles,box,r_cut,b,eps,wca_cut,wpot,nphi) &
!$OMP FIRSTPRIVATE (u)
!$OMP DO

do i= 1, (num_particles-1)
do j= (i+1), num_particles
        rij_length= rij_array(4,i,j)
        if (rij_length <= r_cut) then
                rij(1)= rij_array(1,i,j)
                rij(2)= rij_array(2,i,j)
                rij(3)= rij_array(3,i,j)

                !force due to magnetic dipoles
                unit_rij= rij/rij_length
                m1r= dot_product(m(i,:),unit_rij)
                m2r= dot_product(m(j,:),unit_rij)
                m12= dot_product(m(i,:),m(j,:))
                rsq= dot_product(rij,rij)
                rsqi= 1.0d0/rsq
                fij= m1r*m(j,:) + m2r*m(i,:) + m12*unit_rij - 5.0d0*m1r*m2r*unit_rij
        !       fij= fij*rsqi*rsqi*3.0d0
        !       fij= fij/3.0d0
                fij= fij*rsqi*rsqi
                f(j,:)= f(j,:) + fij
                f(i,:)= f(i,:) - fij

                if (nphi==1) then
                        if (rij_length <= wca_cut) then
                        r6= rsqi*rsqi*rsqi
                        r12= r6**2
                        fij= rij*(24.0d0*r6 - 48.0d0*r12)*rsqi*eps
                        f(i,:)= f(i,:) + fij
                        f(j,:)= f(j,:) - fij
                        u= u + 4.0d0*(r12 - r6) + 1.0d0
                        endif
                endif
        endif
enddo
enddo
!$OMP END DO NOWAIT
!$OMP ATOMIC
wpot= wpot + u*eps
!$OMP END PARALLEL
        select case (nphi)
        case(1)
        mx= m
        case(2)
        my= m
        case(3)
        mz= m
        end select

enddo
end subroutine get_energy_forces


subroutine SelfConsistentEnergy(num_particles,box,rij_array,b0,m,wpot)
!use ifport
implicit none
integer :: num_particles
integer :: i,j
integer :: counter, maxcounter
real(8), dimension (4,num_particles,num_particles) :: rij_array
real(8), dimension (num_particles,3) :: m, m1, m2
real(8), dimension (num_particles,3) :: blocal, mb
real(8), dimension (3) :: b0, box
real(8) :: wpot, reduced_chi, threshold
real(8) :: ne, e
real(8) :: randomx
real(8) :: temp_set, d_time,t1,t2
logical :: converged

reduced_chi= 0.893d0/24.d0
threshold= 10.d0**(-8)

!initializiation
forall (i= 1:num_particles) m(i,:)= b0
blocal= 0.0d0
converged= .false.
e= 1.0d10
counter= 0

!loop until convergence
do while(.not. converged)
        ne= 0.0d0
        call SelfConsistentEnergy_update_moments(num_particles,box,m,b0,rij_array)
        mb= m*m !because m = b0 + blocal the m*m is the same as m*(b0 + blocal)
        ne= 0.5d0*sum(mb)
        if(abs(ne - e) .lt. threshold) then
                converged= .true.
        else
                if (counter>20) then
                        m1(:,:)= m(:,:)
                        call SelfConsistentEnergy_update_moments(num_particles,box,m,b0,rij_array)
                        m2(:,:)= m(:,:)
                        randomx= rand()
                        m= randomx*m1 + (1.0d0 - randomx)*m2
                        ne= 0.0d0
                        mb= m*m !because m = b0 + blocal the m*m is the same as m*(b0 + blocal)
                        ne= 0.5d0*sum(mb)
                        counter= 0
                endif
                e= ne
        end if
        counter= counter+1
end do
wpot= e

end subroutine SelfConsistentEnergy


subroutine SelfConsistentEnergy_update_moments(num_particles,box,m,b0,rij_array)
implicit none
integer :: num_particles
integer :: j, k
real(8), dimension(4,num_particles,num_particles) :: rij_array
real(8), dimension(num_particles,3) :: m
real(8), dimension(num_particles,3) :: blocal
real(8), dimension(3) :: b0, bj, box
real(8), dimension(3) :: rij, unit_rij
real(8) :: reduced_chi, threshold
real(8) :: rij_length
real(8) :: temp_set, d_time
real(8) :: xij,yij,zij,xi,yi,zi,wca_cut

reduced_chi= 0.893d0/24.d0
threshold= 10.d0**(-8)
wca_cut= 2.d0**(1.d0/6.d0)
blocal= 0.0d0

!$OMP PARALLEL DEFAULT (PRIVATE) &
!$OMP FIRSTPRIVATE (m) &
!$OMP SHARED (rij_array,b0,blocal,num_particles,box,reduced_chi,wca_cut)
!$OMP DO
do j= 1, num_particles
        !moment due to external field
        !m(j, :)= b0                            !single-thread code
        !field due to other dipoles
        do k= 1, num_particles
                if(j .eq. k) cycle
                rij(1)= rij_array(1,j,k)
                rij(2)= rij_array(2,j,k)
                rij(3)= rij_array(3,j,k)
                rij_length= rij_array(4,j,k)
                unit_rij= rij / rij_length

                !limit the dipolar field
                if (rij_length < wca_cut) then
                        rij_length= wca_cut
                endif

                bj= 3.0d0*dot_product(m(k, :), unit_rij)*unit_rij
                bj= bj - m(k, :)
                bj= bj / (rij_length**3)
                blocal(j, :)= blocal(j, :) + bj
        end do
        !m(j, :)= m(j, :) + blocal(j, :)        !single-thread code
end do
!$OMP END DO NOWAIT
!$OMP END PARALLEL
forall (j=1:num_particles) m(j,:)= b0(:) + reduced_chi*blocal(j,:) !parallel
end subroutine SelfConsistentEnergy_update_moments


subroutine Create_rij_array(num_particles,box,r,rij_array)
implicit none
integer :: i,j,num_particles
real(8), dimension (4,num_particles,num_particles) :: rij_array
real(8), dimension (num_particles,3) :: r
real(8), dimension (3) :: box, ri, rij
real(8) :: temp_set,d_time
real(8) :: xi,yi,zi,xij,yij,zij,rij_length

!$OMP PARALLEL DEFAULT (PRIVATE) &
!$OMP SHARED (rij_array,r,num_particles,box)
!$OMP DO
do i= 1, num_particles-1
        xi= r(i,1)
        yi= r(i,2)
        zi= r(i,3)
        do j= i+1, num_particles
                xij= r(j,1)-xi
                yij= r(j,2)-yi
                zij= r(j,3)-zi
                call minimum_image(xij,yij,zij,box)
                rij_array(1,i,j)= xij
                rij_array(2,i,j)= yij
                rij_array(3,i,j)= zij
                rij_array(1,j,i)= -xij
                rij_array(2,j,i)= -yij
                rij_array(3,j,i)= -zij
                rij_length= sqrt(xij**2 + yij**2 + zij**2)
                rij_array(4,i,j)= rij_length
                rij_array(4,j,i)= rij_length

        enddo
enddo
!$OMP END DO NOWAIT
!$OMP END PARALLEL
end subroutine Create_rij_array


subroutine minimum_image(xij,yij,zij,box)
implicit none
real(8) :: xij,yij,zij
real(8), dimension (3) :: box,boxl

boxl= 0.50d0*box
if (xij > boxl(1)) then
        xij= xij - box(1)
elseif (xij < -boxl(1)) then
        xij= xij + box(1)
endif
if (yij > boxl(2)) then
        yij= yij - box(2)
elseif (yij < -boxl(2)) then
        yij= yij + box(2)
endif
end subroutine minimum_image

subroutine generate_brms(theta,brms)
implicit none
integer :: nphi, i
real(8) :: phi, theta, pi
real(8) :: bxsq, bysq, bzsq
real(8), dimension (3) :: b0
real(8), dimension (3,3) :: brms

nphi= 100
bxsq= 0.0d0
bysq= 0.0d0
bzsq= 0.0d0
pi= 4.0d0*atan(1.0d0)

do i= 1, nphi
        phi= (i-1)*2.0d0*pi/real(nphi)
        b0= (/ sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta) /)
        bxsq= bxsq + b0(1)**2
        bysq= bysq + b0(2)**2
        bzsq= bzsq + b0(3)**2
enddo

bxsq= sqrt(bxsq/real(nphi))
bysq= sqrt(bysq/real(nphi))
bzsq= sqrt(bzsq/real(nphi))

brms(1,:)= (/ bxsq, 0.0d0, 0.0d0 /)
brms(2,:)= (/ 0.0d0, bysq, 0.0d0 /)
brms(3,:)= (/ 0.0d0, 0.0d0, bzsq /)
end subroutine generate_brms

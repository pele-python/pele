module rotations
    double precision epsilon
    parameter (epsilon = 1d-6)

contains

    ! multiply 2 quaternions q1, q2
    function rot_q_multiply(q1, q2) result(q3)
        implicit none
        double precision, intent(in) :: q1(4), q2(4)
        double precision :: q3(4)
    
        q3(1) = q1(1)*q2(1)-q1(2)*q2(2)-q1(3)*q2(3)-q1(4)*q2(4)
        q3(2) = q1(1)*q2(2)+q1(2)*q2(1)+q1(3)*q2(4)-q1(4)*q2(3)
        q3(3) = q1(1)*q2(3)-q1(2)*q2(4)+q1(3)*q2(1)+q1(4)*q2(2)
        q3(4) = q1(1)*q2(4)+q1(2)*q2(3)-q1(3)*q2(2)+q1(4)*q2(1)
    end function

    ! change a given angle axis rotation p1 by the
    ! rotation p2
    function rot_rotate_aa(p1, p2) result(p3)
        implicit none
        double precision, intent(in) :: p1(3), p2(3)
        double precision :: p3(3)

        p3 = rot_q2aa(rot_q_multiply(rot_aa2q(p2), rot_aa2q(p1)))
    end function

    ! convert angle axis to quaternion
    function rot_aa2q(p) result(q)
        use vec3
        implicit none

        double precision, intent(in) :: p(3)
        double precision :: q(4)
        double precision :: thetah

        thetah = 0.5d0 * vec_len(p)
        q(1) = cos(thetah)

        ! do linear expansion for small epsilon
        if(thetah < epsilon) then
            q(2:4) = 0.5d0 * p
        else
            q(2:4) = 0.5d0 * sin(thetah) * p / thetah
        endif
        ! make sure to have normal form
        if(q(1) < 0d0) q = -q
    end function

    function rot_aa2mx(p) result(m)
        implicit none
        double precision, intent(in) :: p
        double precision m(3,3), tmp(3,3)
!        call rmdrvt(p, m, tmp, tmp, tmp, .false.)
    end function

    ! convert quaternion to angle axis
    function rot_q2aa(qin) result(p)
        use vec3
        implicit none
        double precision, intent(in) :: qin(4)
        double precision :: p(3)
        double precision :: theta, s
        double precision :: q(4)

        q=qin
        if (q(1) < 0d0) q=-q
        ! avoid rounding errors
        if (q(1) > 1.0d0) q=q/sqrt(dot_product(q,q))
        theta = 2d0 * acos(q(1))
        s = sqrt(1d0-q(1)*q(1))
        if (s < 1d-6) then
            p = 2d0 * q(2:4)
        else
            p = q(2:4) / s * theta
        endif
    end function

    function rot_q2mx(qin) result(m)
        implicit none
        double precision, intent(in) :: qin(4)
        double precision :: m(3,3)
        double precision :: q(4)
        double precision :: sq(4), tmp1, tmp2
        integer i

        q = qin / sqrt(dot_product(qin,qin))

        do i=1,4
            sq(i) = q(i)*q(i)
        enddo

        m(1,1) = ( sq(2) - sq(3) - sq(4) + sq(1))
        m(2,2) = (-sq(2) + sq(3) - sq(4) + sq(1))
        m(3,3) = (-sq(2) - sq(3) + sq(4) + sq(1))

        tmp1 = q(2)*q(3)
        tmp2 = q(1)*q(4)
        m(2,1) = 2.0d0 * (tmp1 + tmp2)
        m(1,2) = 2.0d0 * (tmp1 - tmp2)

        tmp1 = q(2)*q(4)
        tmp2 = q(3)*q(1)
        m(3,1) = 2.0d0 * (tmp1 - tmp2)
        m(1,3) = 2.0d0 * (tmp1 + tmp2)
        tmp1 = q(3)*q(4)
        tmp2 = q(1)*q(2)
        m(3,2) = 2.0d0 * (tmp1 + tmp2)
        m(2,3) = 2.0d0 * (tmp1 - tmp2)
    end function

    ! convert matrix to quaternionfunction
    function rot_mx2q(mi) result(q)
        implicit none
        double precision, intent(in) :: mi(3,3)
        double precision :: q(4), m(3,3)
        double precision :: trace, s

        m = transpose(mi)
        trace=m(1,1)+m(2,2)+m(3,3)

        if (trace > 0d0) then
            s = sqrt(trace+1.0d0) * 2d0;
            q(1) = 0.25d0 * S;
            q(2) = (m(2,3) - m(3,2)) / s
            q(3) = (m(3,1) - m(1,3)) / s
            q(4) = (m(1,2) - m(2,1)) / s
        else if ((m(1,1) > m(2,2)).and.(m(1,1) > m(3,3))) then
            s=sqrt(1.0 + m(1,1) - m(2,2) - m(3,3)) * 2d0
            q(1) = (m(2,3) - m(3,2)) / s
            q(2) = 0.25d0 * s
            q(3) = (m(2,1) + m(1,2)) / s
            q(4) = (m(3,1) + m(1,3)) / s
        else if (m(2,2) > m(3,3)) then
            s = sqrt(1.0 + m(2,2) - m(1,1) - m(3,3)) * 2d0
            q(1) = (m(3,1) - m(1,3)) / s
            q(2) = (m(2,1) + m(1,2)) / s
            q(3) = 0.25d0 * s
            q(4) = (m(3,2) + m(2,3)) / s
        else
            s = sqrt(1.0 + m(3,3) - m(1,1) - m(2,2)) * 2d0
            q(1) = (m(1,2) - m(2,1)) / s
            q(2) = (m(3,1) + m(1,3)) / s
            q(3) = (m(3,2) + m(2,3)) / s
            q(4) = 0.25d0 * s
        endif

        if(q(1) < 0) q = -q
    end function

    ! convert matrix to angle axis
    function rot_mx2aa(m) result(p)
        implicit none
        double precision, intent(in) :: m(3,3)
        double precision :: p(3)

        p = rot_q2aa(rot_mx2q(m))
    end function

    ! uniform random rotation in angle axis formulation
    ! input: 3 uniformly distributed random numbers
    ! uses the algorithm given in
    !  K. Shoemake, Uniform random rotations, Graphics Gems III, pages 124-132. Academic, New York, 1992.
    ! This first generates a random rotation in quaternion representation. We should substitute this by
    ! a direct angle axis generation, but be careful: the angle of rotation in angle axis representation
    ! is NOT uniformly distributed
    function rot_random_q() result(q)
        implicit none
        double precision :: q(4)
        double precision :: u(3)
        double precision :: PI
        double precision :: dprand
        parameter (pi=3.141592654d0)

        u = (/dprand(),dprand(),dprand()/)
        q(1) = sqrt(1d0-u(1)) * sin(2d0*PI*u(2))
        q(2) = sqrt(1d0-u(1)) * cos(2d0*PI*u(2))
        q(3) = sqrt(u(1)) * sin(2d0*PI*u(3))
        q(4) = sqrt(u(1)) * cos(2d0*PI*u(3))
    end function

    ! generate a random rotation in angle axis
    function rot_random_aa() result(p)
        implicit none
        double precision :: p(3)
        p = rot_q2aa(rot_random_q())
    end function

    ! generate a small random rotation
    function rot_small_random_aa(maxtheta) result(p)
        use vec3
        implicit none
        double precision, intent(in) :: maxtheta
        double precision :: p(3)
        double precision :: u, dprand
        double precision :: s

        ! first choose a random unit vector
        p = vec_random()

        ! linear for too small steps
        ! this is not completely right but should be ok
        if(maxtheta < epsilon) then
            p = p*dprand()*maxtheta
            return
        endif

        s = 1d0 / (sin(0.5d0*maxtheta)**2)
        ! now choose the angle theta in range 0:step
        ! with distribution sin(0.5*theta)**2
        u=dprand()*maxtheta
        do while( dprand() > s*sin(0.5d0*u)**2)
            u=dprand()*maxtheta
        end do
        p = p*u
    end function

    ! change an angle axis vector by a small rotation
    subroutine rot_takestep_aa(p, maxtheta)
        implicit none
        double precision, intent(inout) :: p(3)
        double precision, intent(in) :: maxtheta

        p = rot_rotate_aa(p, rot_small_random_aa(maxtheta))
    end subroutine

    subroutine rot_construct_unit_system(x, e)
        use vec3
        implicit none
        double precision x(3,3)
        double precision e(3,3)
        double precision tmp(3)
        integer i

        ! construct coordinate system of current
        e(:,1) = x(:,2) - x(:,1)
        tmp = x(:,3) - x(:,1)
        e(:,2) = vec_cross(e(:,1), tmp)
        e(:,3) = vec_cross(e(:,1), e(:,2))

        ! now normalize the vectors
        do i=1,3
            e(:,i) = e(:,i) / vec_len(e(:,i))
        end do
    end subroutine

    ! calculate orientation in angle axis vector based
    ! on 3 point x(:,1), x(:,2), x(:,3) and 3 unrotated
    ! reference points in xref
    function rot_get_orientation_aa(x, xref) result(p)
        use vec3
        implicit none
        double precision, intent(in) :: x(3,3), xref(3,3)
        double precision p(3)
        double precision eref(3,3), e(3,3), minv(3,3)

        call rot_construct_unit_system(xref, eref)
        call rot_construct_unit_system(x, e)
        call invert3x3(eref, minv)

        ! construct rotation based on
        ! e = R*eref => R = e*eref^-1
        p = rot_mx2aa(matmul(e,minv))
    end function

end module

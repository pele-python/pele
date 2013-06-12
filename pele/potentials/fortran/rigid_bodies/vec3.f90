! vr274> this module contains various helper function for handling 
!        3 dimensional vectors + matrices. This can be useful for
!        coordinate handling.
module vec3
contains

    ! length of vector
    function vec_len(v)
        implicit none
        double precision, intent(in) :: v(3)
        double precision vec_len
        vec_len = sqrt(dot_product(v,v))
    end function

    ! angle between 2 vectors in radians
    function vec_angle(v1, v2)
        implicit none
        double precision, intent(in) :: v1(3), v2(3)
        double precision vec_angle
        vec_angle = acos(dot_product(v1, v2) / (vec_len(v1) * vec_len(v2)))
    end function

    ! cross product of 2 vectors
    function vec_cross(v1, v2) result(v3)
        implicit none
        double precision, intent(in) :: v1(3), v2(3)
        double precision v3(3)
        v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
        v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
        v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
    end function

    ! dyadic product of 2 vectors
    function vec_dyad(v1, v2) result(m)
        implicit none
        double precision, intent(in) :: v1(3), v2(3)
        double precision m(3,3)
        integer i,j
        m=0d0
        do i=1,3
          do j=i,3
             m(i,j) = v1(i)*v2(j)
          end do
        end do
    end function

    ! uniform random unit vector 
    function vec_random() result(p)
        double precision :: p(3)
        double precision :: PI
        double precision :: dprand
        parameter (pi=3.141592654d0)
        double precision z
        double precision u(2)
        u = (/dprand(),dprand()/)
        z = 2*u(1) - 1d0
        p(1) = sqrt(1-z*z) * cos(2d0*PI*u(2))
        p(2) = sqrt(1-z*z) * sin(2d0*PI*u(2))
        p(3) = z
    end function

    ! invert a 3x3 matrix
    subroutine invert3x3 (A, Ainv)
        double precision, intent(in) :: A(3,3)
        double precision, intent(out) ::Ainv(3,3)
        double precision Adet

        Ainv = 0d+0;
        Ainv(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
        Ainv(1,2) = A(3,2)*A(1,3) - A(1,2)*A(3,3)
        Ainv(1,3) = A(1,2)*A(2,3) - A(1,3)*A(2,2)

        Adet = Ainv(1,1)*A(1,1) + Ainv(1,2)*A(2,1) + Ainv(1,3)*A(3,1)

        Ainv(1,1) = Ainv(1,1)/Adet
        Ainv(1,2) = Ainv(1,2)/Adet
        Ainv(1,3) = Ainv(1,3)/Adet

        Ainv(2,1) = (A(2,3)*A(3,1) - A(2,1)*A(3,3))/Adet
        Ainv(2,2) = (A(1,1)*A(3,3) - A(3,1)*A(1,3))/Adet
        Ainv(2,3) = (A(2,1)*A(1,3) - A(1,1)*A(2,3))/Adet

        Ainv(3,1) = (A(2,1)*A(3,2) - A(2,2)*A(3,1))/Adet
        Ainv(3,2) = (A(3,1)*A(1,2) - A(1,1)*A(3,2))/Adet
        Ainv(3,3) = (A(1,1)*A(2,2) - A(1,2)*A(2,1))/Adet
        return
    end subroutine invert3x3

    function det3x3(A) result(Adet)
        implicit none
        double precision, intent(in) :: A(3,3)
        double precision Adet, tmp(3)

        tmp(1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
        tmp(2) = A(3,2)*A(1,3) - A(1,2)*A(3,3)
        tmp(3) = A(1,2)*A(2,3) - A(1,3)*A(2,2)

        Adet = tmp(1)*A(1,1) + tmp(2)*A(2,1) + tmp(3)*A(3,1)
    end function


    subroutine identity3x3(A)
        double precision, intent(out) :: A(3,3)
        A=0d0
        A(1,1)=1d0
        A(2,2)=1d0
        A(3,3)=1d0
    end subroutine
end module vec3

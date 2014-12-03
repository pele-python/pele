!
!js850> WARNING: f2py by default copies every multidimentional array passed to
!fortran.  This is because array indexing is different in fortran and python.
!This is a huge waste of resources
!

!     --------------------------------------------------------------------------

!     RMDRVT = rotation matrix derivative
!     P(3) = rotation vector
!     RM(3,3) = rotation matrix
!     DRMk(3,3) = derivative of the rotation matrix with respect to the
!                 kth component of the rotation vector
!     GTEST = true if derivatives are to be found

SUBROUTINE RMDRVT(P, RM, DRM1, DRM2, DRM3, GTEST)

    IMPLICIT NONE

    !     PN(3) = the unit vector parallel to P
    !     THETA = the modulus of the rotation vector P, equivalent to the
    !             angle of rotation
    !     THETA2 = THETA squared
    !     THETA3 = THETA**-3
    !     CT = cos(THETA)
    !     ST = sin(THETA)
    !     I3(3,3) = 3x3 identity matrix
    !     E(3,3) = the skew-symmetric matrix obtained from a unit vector
    !              parallel to P (equation (2) in the paper)
    !     ESQ(3,3) = the square of E
    DOUBLE PRECISION, INTENT(IN) :: P(3)
    DOUBLE PRECISION, INTENT(OUT) :: RM(3,3), DRM1(3,3), DRM2(3,3), DRM3(3,3)
    DOUBLE PRECISION :: PN(3), THETA, THETA2, THETA3, CT, ST, I3(3,3), E(3,3), ESQ(3,3)

    !     DEk = derivate of E with respect to the kth component of P
    DOUBLE PRECISION :: DE1(3,3), DE2(3,3), DE3(3,3)
    LOGICAL, INTENT(IN)   :: GTEST

    !     Set the values of the identity matrix I3
    I3(:,:) = 0.D0
    I3(1,1) = 1.D0; I3(2,2) = 1.D0; I3(3,3) = 1.D0

    !     Calculate the value of THETA2 as the square modulus of P
    THETA2  = DOT_PRODUCT(P,P)

    IF (THETA2 < 1.0D-12) THEN
        !        Execute if the angle of rotation is zero
        !        In this case the rotation matrix is the identity matrix
        RM(:,:) = I3(:,:)

        ! vr274> first order corrections to rotation matrix
        RM(1,2) = -P(3)
        RM(2,1) = P(3)
        RM(1,3) = P(2)
        RM(3,1) = -P(2)
        RM(2,3) = -P(1)
        RM(3,2) = P(1)

        !        If derivatives do not need to be found, we're finished
        IF (.NOT. GTEST) RETURN

        !        This is the special case described in the paper, where DRMk is
        !        equal to E(k), which is the skew-symmetric matrix obtained from
        !        P with P(k) equal to 1 and other components equal to zero
        !         PN        = (/1.D0, 0.D0, 0.D0/)
        !         E(:,:)    = 0.D0
        !         E(2,3)    = -PN(1)
        !         E(3,2)    = -E(2,3)
        !         DRM1(:,:) = E(:,:)

        !         PN        = (/0.D0, 1.D0, 0.D0/)
        !         E(:,:)    = 0.D0
        !         E(1,3)    =  PN(2)
        !         E(3,1)    = -E(1,3)
        !         DRM2(:,:) = E(:,:)

        !         PN        = (/0.D0, 0.D0, 1.D0/)
        !         E(:,:)    = 0.D0
        !         E(1,2)    = -PN(3)
        !         E(2,1)    = -E(1,2)
        !         DRM3(:,:) = E(:,:)

        ! hk286 - now up to the linear order in theta
        E(:,:)    = 0.D0
        E(1,1)    = 0.0D0
        E(1,2)    = P(2)
        E(1,3)    = P(3)
        E(2,1)    = P(2)
        E(2,2)    = -2.0D0*P(1)
        E(2,3)    = -2.0D0
        E(3,1)    = P(3)
        E(3,2)    = 2.0D0
        E(3,3)    = -2.0D0*P(1)
        DRM1(:,:) = 0.5D0*E(:,:)

        E(:,:)    = 0.D0
        E(1,1)    = -2.0D0*P(2)
        E(1,2)    = P(1)
        E(1,3)    = 2.0D0
        E(2,1)    = P(1)
        E(2,2)    = 0.0D0
        E(2,3)    = P(3)
        E(3,1)    = -2.0D0
        E(3,2)    = P(3)
        E(3,3)    = -2.0D0*P(2)
        DRM2(:,:) = 0.5D0*E(:,:)

        E(:,:)    = 0.D0
        E(1,1)    = -2.0D0*P(3)
        E(1,2)    = -2.0D0
        E(1,3)    = P(1)
        E(2,1)    = 2.0D0
        E(2,2)    = -2.0D0*P(3)
        E(2,3)    = P(2)
        E(3,1)    = P(1)
        E(3,2)    = P(2)
        E(3,3)    = 0.0D0
        DRM3(:,:) = 0.5D0*E(:,:)

    ELSE
        !       Execute for the general case, where THETA dos not equal zero
        !       Find values of THETA, CT, ST and THETA3
        THETA   = SQRT(THETA2)
        CT      = COS(THETA)
        ST      = SIN(THETA)
        THETA3  = 1.D0/(THETA2*THETA)

        !        Set THETA to 1/THETA purely for convenience
        THETA   = 1.D0/THETA

        !        Normalise P and construct the skew-symmetric matrix E
        !        ESQ is calculated as the square of E
        PN(:)   = THETA*P(:)
        E(:,:)  = 0.D0
        E(1,2)  = -PN(3)
        E(1,3)  =  PN(2)
        E(2,3)  = -PN(1)
        E(2,1)  = -E(1,2)
        E(3,1)  = -E(1,3)
        E(3,2)  = -E(2,3)
        ESQ     = MATMUL(E,E)

        !        RM is calculated from Rodrigues' rotation formula (equation (1)
        !        in the paper)
        RM      = I3(:,:) + (1.D0-CT)*ESQ(:,:) + ST*E(:,:)

        !        If derivatives do not need to found, we're finished
        IF (.NOT. GTEST) RETURN

        !        Set up DEk using the form given in equation (A4) in the paper
        DE1(:,:) = 0.D0
        DE1(1,2) = P(3)*P(1)*THETA3
        DE1(1,3) = -P(2)*P(1)*THETA3
        DE1(2,3) = -(THETA - P(1)*P(1)*THETA3)
        DE1(2,1) = -DE1(1,2)
        DE1(3,1) = -DE1(1,3)
        DE1(3,2) = -DE1(2,3)

        DE2(:,:) = 0.D0
        DE2(1,2) = P(3)*P(2)*THETA3
        DE2(1,3) = THETA - P(2)*P(2)*THETA3
        DE2(2,3) = P(1)*P(2)*THETA3
        DE2(2,1) = -DE2(1,2)
        DE2(3,1) = -DE2(1,3)
        DE2(3,2) = -DE2(2,3)

        DE3(:,:) = 0.D0
        DE3(1,2) = -(THETA - P(3)*P(3)*THETA3)
        DE3(1,3) = -P(2)*P(3)*THETA3
        DE3(2,3) = P(1)*P(3)*THETA3
        DE3(2,1) = -DE3(1,2)
        DE3(3,1) = -DE3(1,3)
        DE3(3,2) = -DE3(2,3)

        !        Use equation (A3) in the paper to find DRMk
        DRM1(:,:) = ST*PN(1)*ESQ(:,:) + (1.D0-CT)*(MATMUL(DE1,E) + MATMUL(E,DE1)) &
            + CT*PN(1)*E(:,:) + ST*DE1(:,:)

        DRM2(:,:) = ST*PN(2)*ESQ(:,:) + (1.D0-CT)*(MATMUL(DE2,E) + MATMUL(E,DE2)) &
            + CT*PN(2)*E(:,:) + ST*DE2(:,:)

        DRM3(:,:) = ST*PN(3)*ESQ(:,:) + (1.D0-CT)*(MATMUL(DE3,E) + MATMUL(E,DE3)) &
            + CT*PN(3)*E(:,:) + ST*DE3(:,:)

    ENDIF
END SUBROUTINE RMDRVT



function sitedist(drij, p1, p2, S, W, cog) result(dist)
    implicit none
    double precision, intent(in) :: drij(3), p1(3), p2(3)
    double precision :: DR1(3,3), DR2(3,3), DR3(3,3)
    double precision, intent(in) :: S(3,3), W, cog(3)
    double precision dist

    double precision R1(3,3), R2(3,3), dR(3,3)
    double precision d_M, d_P, d_mix

    ! Note, we don't need to calculate the derivatives this time
    call RMDRVT(p1, R1, DR1, DR2, DR3, .FALSE.)
    call RMDRVT(p2, R2, DR1, DR2, DR3, .FALSE.)

    dR = R2 - R1

    d_M = W*sum((drij)**2)
    DR1 = matmul(dR, matmul(S, transpose(dR)))
    d_P = DR1(1,1) + DR1(2,2) + DR1(3,3)
    d_mix = 2. * W * dot_product(drij, matmul(dR, cog))

    dist = d_M + d_P + d_mix
end function

function trace3(M) result(tr)
    implicit none
    double precision, intent(in) :: M(3,3)
    double precision tr
    tr = M(1,1) + M(2,2) + M(3,3)
end function

subroutine sitedist_grad(drij, p1, p2, S, W, cog, g_M, g_P)
    implicit none
    double precision, intent(in) :: drij(3), p1(3), p2(3)
    double precision :: R11(3,3), R12(3,3), R13(3,3)
    double precision, intent(in) :: S(3,3), W, cog(3)
    double precision, intent(out) :: g_M(3), g_P(3)

    double precision R1(3,3), R2(3,3), dR(3,3)
    double precision trace3

    call RMDRVT(p2, R2, R11, R12, R13, .FALSE.)
    call RMDRVT(p1, R1, R11, R12, R13, .TRUE.)

    dR = R2 - R1

    g_M = -2.*W*(drij)

    g_P(1) = -2.*trace3(matmul(R11, matmul(S, transpose(dR))))
    g_P(2) = -2.*trace3(matmul(R12, matmul(S, transpose(dR))))
    g_P(3) = -2.*trace3(matmul(R13, matmul(S, transpose(dR))))

    g_M = g_M - 2.*W *  matmul(dR, cog)
    g_P(1) = g_P(1) - 2.*W * dot_product(drij, matmul(R11, cog))
    g_P(2) = g_P(2) - 2.*W * dot_product(drij, matmul(R12, cog))
    g_P(3) = g_P(3) - 2.*W * dot_product(drij, matmul(R13, cog))


end subroutine

!function aadist(coords1, coords2, nrigid, S, W, cog)
!    implicit none
!    double precision, intent(in) :: coords1(6*nrigid)
!    double precision, intent(in) :: coords2(6*nrigid)
!
!    double precision aadist
!    integer, intent(in) :: nrigid
!    double precision sitedist
!
!    double precision dist
!    integer i, icom, ip
!    double precision, intent(in) :: S(3,3), W, cog(3)
!    dist=0.
!
!    do i = 1,nrigid
!        icom = 3*i-2
!        ip = 3*i-2 + 3*nrigid
!
!        dist = dist + sitedist( coords1(icom:icom+3), coords1(ip:ip+3), &
!                                coords2(icom:icom+3), coords2(ip:ip+3), &
!                                S, W, cog )
!    enddo
!    aadist = dist
!end function


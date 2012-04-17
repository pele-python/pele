!module commons
!  double precision stockmu,efield
!  logical efieldt
!end module

SUBROUTINE STOCKAA (natoms,X, G, ENERGY,stockmu)

      ! USE COMMONS, ONLY: STOCKMU, EFIELDT, EFIELD
      IMPLICIT NONE
integer natoms
      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, J7, J8, REALNATOMS, OFFSET 
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: ENERGY, R2, R4, R6, R12, ABSRIJ, RIJSQ, DPFCT, DVDR
      DOUBLE PRECISION :: RI(3), RJ(3), RIJ(3), NR(3), P(3), DU(3), EI(3), EJ(3)
      DOUBLE PRECISION :: DR1(NATOMS/2,3), DR2(NATOMS/2,3), DR3(NATOMS/2,3) 
      DOUBLE PRECISION :: E(NATOMS/2,3), DE1(NATOMS/2,3), DE2(NATOMS/2,3), DE3(NATOMS/2,3)
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      DOUBLE PRECISION :: ALP, BET, GAM, VR, VA, VB, VG, FIJN, FIJEI, FIJEJ, FIJ(3) 
      LOGICAL          :: GTEST = .true.
      double precision :: stockmu,efield
      logical :: efieldt
      
      CALL DEFSTOCK(STOCKMU, DU, DPFCT)

      ENERGY  = 0.D0
      efieldt=.false.
      stockmu=0.5
      IF (GTEST) G(:) = 0.D0

      REALNATOMS = NATOMS/2
      OFFSET     = 3*REALNATOMS
  
      DO J1 = 1, REALNATOMS

         J3 = 3*J1
         J5 = OFFSET + J3
         RI = X(J3-2:J3)
         P  = X(J5-2:J5)

         CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, GTEST)

         E(J1,:)   = MATMUL(RMI(:,:),DU(:))
 
         IF (GTEST) THEN

            DE1(J1,:) = MATMUL(DRMI1(:,:),DU(:))
            DE2(J1,:) = MATMUL(DRMI2(:,:),DU(:))
            DE3(J1,:) = MATMUL(DRMI3(:,:),DU(:))

         ENDIF

      ENDDO

      DO J1 = 1, REALNATOMS  

         J3 = 3*J1
         J5 = OFFSET + J3
 
         RI(:)  = X(J3-2:J3)
         EI(:)  = E(J1,:)

         DO J2 = J1 + 1, REALNATOMS

            J4 = 3*J2
            J6 = OFFSET + J4

!     LJ CONTRIBUTION
             
            RJ(:)  = X(J4-2:J4)
            RIJ(:) = RI(:) - RJ(:)
            RIJSQ  = DOT_PRODUCT(RIJ(:),RIJ(:))
            R2     = 1.D0/RIJSQ
            R6     = R2*R2*R2
            R12    = R6*R6
            ABSRIJ = DSQRT(RIJSQ)
            NR(:)  = RIJ(:)/ABSRIJ
            
            ENERGY = ENERGY + 4.D0*(R12 - R6)
            
            IF (GTEST) THEN

               DVDR = 4.D0*(-12.D0*R12 + 6.D0*R6)*R2

               G(J3-2:J3)  = G(J3-2:J3) + DVDR*RIJ(:)
               G(J4-2:J4)  = G(J4-2:J4) - DVDR*RIJ(:)

            ENDIF

!     DIPOLAR CONTRIBUTION

            R4     = R2*R2
            EJ(:)  = E(J2,:)
            ALP    = DOT_PRODUCT(NR(:),EI(:))
            BET    = DOT_PRODUCT(NR(:),EJ(:))
            GAM    = DOT_PRODUCT(EI(:),EJ(:))

            ENERGY = ENERGY + DPFCT*R2*(GAM/3.D0 - ALP*BET)/ABSRIJ

            IF (GTEST) THEN

               VR     = -DPFCT*R4*(GAM - 3.D0*ALP*BET)
               VA     = -DPFCT*BET*R2/ABSRIJ
               VB     = -DPFCT*ALP*R2/ABSRIJ
               VG     =  DPFCT*R2/(3.D0*ABSRIJ)

               FIJN   = VR - (VA*ALP+VB*BET)/ABSRIJ
               FIJEI  = VA/ABSRIJ
               FIJEJ  = VB/ABSRIJ
               FIJ(:) = FIJN*NR(:) + FIJEI*EI(:) + FIJEJ*EJ(:)

               G(J3-2:J3) = G(J3-2:J3) + FIJ(:)
               G(J4-2:J4) = G(J4-2:J4) - FIJ(:)

               G(J5-2) = G(J5-2) + VA*DOT_PRODUCT(NR(:),DE1(J1,:)) + VG*DOT_PRODUCT(DE1(J1,:),EJ(:))
               G(J5-1) = G(J5-1) + VA*DOT_PRODUCT(NR(:),DE2(J1,:)) + VG*DOT_PRODUCT(DE2(J1,:),EJ(:))
               G(J5)   = G(J5)   + VA*DOT_PRODUCT(NR(:),DE3(J1,:)) + VG*DOT_PRODUCT(DE3(J1,:),EJ(:))

               G(J6-2) = G(J6-2) + VB*DOT_PRODUCT(NR(:),DE1(J2,:)) + VG*DOT_PRODUCT(EI(:),DE1(J2,:))
               G(J6-1) = G(J6-1) + VB*DOT_PRODUCT(NR(:),DE2(J2,:)) + VG*DOT_PRODUCT(EI(:),DE2(J2,:))
               G(J6)   = G(J6)   + VB*DOT_PRODUCT(NR(:),DE3(J2,:)) + VG*DOT_PRODUCT(EI(:),DE3(J2,:))

            ENDIF 
 
         ENDDO

         IF (EFIELDT) THEN

            ENERGY = ENERGY - STOCKMU*EFIELD*EI(3)

            IF (GTEST) THEN

               G(J5-2) = G(J5-2) - STOCKMU*EFIELD*DE1(J1,3)
               G(J5-1) = G(J5-1) - STOCKMU*EFIELD*DE2(J1,3)
               G(J5)   = G(J5)   - STOCKMU*EFIELD*DE3(J1,3)

            ENDIF

         ENDIF

      ENDDO

      END SUBROUTINE STOCKAA
 
!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFSTOCK(STOCKMU, DU, DPFCT)
    
      IMPLICIT NONE

      DOUBLE PRECISION :: STOCKMU, DPFCT, DU(3)

      DU(:)   = (/0.D0, 0.D0, 1.D0/)
      DPFCT   = 3.D0*STOCKMU*STOCKMU

      END SUBROUTINE DEFSTOCK

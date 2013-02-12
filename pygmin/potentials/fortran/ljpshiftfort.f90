!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of GMIN.
!
!   GMIN is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   GMIN is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!*************************************************************************
!
!  Subroutine LJPSHIFT calculates the energy, cartesian gradient and second
!  derivative matrix analytically for Lennard-Jones in reduced units
!  (epsilon=sigma=1) using a shifted, truncated potential.
!
!  Adapted for the binary LJ glass described by Sastry, Debenetti and
!  Stillinger, Nature, 393, 554, 1998. Atom types are A and B. The first
!  NTYPEA are A, the next NBTYPE=NATOMS-NTYPEA are B. epsilon and sigma for A are the
!  units of energy and distance, so we also need EPSAB, EPSAB, SIGAB and
!  SIGAA in these units. Sastry et al. density is 1.2 i.e. a box length
!  of 5.975206 for 256 atoms. 
!
!
!*************************************************************************
!*******************************************************************
!
      SUBROUTINE LJPSHIFT(X, V, POTEL, GTEST, STEST,&
         NATOMS, BOXLX, BOXLY, BOXLZ, CUTOFF, PERIODIC, NTYPEA,&
         EPSAB, EPSBB, SIGAB, SIGBB)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(INOUT) :: X(3*NATOMS) 
      DOUBLE PRECISION, INTENT(OUT) :: V(3*NATOMS), POTEL 
      LOGICAL, INTENT(IN) :: GTEST, STEST, periodic
      DOUBLE PRECISION, INTENT(IN) :: BOXLX, BOXLY, BOXLZ
      DOUBLE PRECISION, INTENT(IN) :: EPSAB, EPSBB, SIGAB, SIGBB
      DOUBLE PRECISION, INTENT(IN) :: CUTOFF
      INTEGER, INTENT(IN) :: NATOMS, NTYPEA
      DOUBLE PRECISION SIGAB6, SIGAB12, SIGBB6, SIGBB12, &
                       IRCUT2AA, IRCUT2AB, IRCUT2BB, SIGRCAA6, &
                       SIGRCAA12, &
                       SIGRCAB6, SIGRCAB12, SIGRCBB6, SIGRCBB12, &
                       CONSTAA, CONSTBB, CONSTAB, &
                       RCONSTAA, RCONSTAB, RCONSTBB,  &
                       CUTAA, CUTAB, CUTBB, &
                       EPSAA, SIGAA6
      INTEGER N, J1, J2

      N=NATOMS
      CUTAA=CUTOFF
      CUTAB=CUTOFF*SIGAB
      CUTBB=CUTOFF*SIGBB
      IRCUT2AA = 1.D0/CUTAA**2
      IRCUT2AB = 1.D0/CUTAB**2
      IRCUT2BB = 1.D0/CUTBB**2
      SIGAB6=SIGAB**6
      SIGAB12=SIGAB6**2
      SIGBB6=SIGBB**6
      SIGBB12=SIGBB6**2
      SIGRCAA6= 1.0D0/CUTAA**6
      SIGRCAA12=SIGRCAA6**2
      SIGRCAB6=SIGAB6/CUTAB**6
      SIGRCAB12=SIGRCAB6**2
      SIGRCBB6=SIGBB6/CUTBB**6
      SIGRCBB12=SIGRCBB6**2
      CONSTAA=4.0D0*SIGRCAA6-7.0D0*SIGRCAA12
      CONSTAB=4.0D0*SIGRCAB6-7.0D0*SIGRCAB12
      CONSTBB=4.0D0*SIGRCBB6-7.0D0*SIGRCBB12
      RCONSTAA=(6.0D0*SIGRCAA12-3.0D0*SIGRCAA6)/CUTAA**2
      RCONSTAB=(6.0D0*SIGRCAB12-3.0D0*SIGRCAB6)/CUTAB**2
      RCONSTBB=(6.0D0*SIGRCBB12-3.0D0*SIGRCBB6)/CUTBB**2

      SIGAA6=1.0D0
      EPSAA=1.0D0
!
!  Work out cutoff for potential. Two particles interact if r<c, but
!  we will use the equivalent condition 1/r^2 > 1/c^2.
!
!  Deal with any atoms that have left the box.
!
      IF (PERIODIC) THEN
         DO J1=1,N
            J2 = 3*(J1-1)
            X(J2+1)=X(J2+1) - BOXLX*ANINT(X(J2+1)/BOXLX)
            X(J2+2)=X(J2+2) - BOXLY*ANINT(X(J2+2)/BOXLY)
            X(J2+3)=X(J2+3) - BOXLZ*ANINT(X(J2+3)/BOXLZ)
         ENDDO
      ENDIF
!
!  Calculate interatomic vectors using the minimum image convention.
!  VEC(i,j,alpha) is the alpha (x, y or z) component of the vector between
!  atoms i and j.
!

      POTEL=0.0D0

      IF ( GTEST .OR. STEST ) THEN

        DO J1=1,3*NATOMS
          V(J1) = 0
        END DO

        DO J1=1,NTYPEA
          DO J2=J1+1,NTYPEA
            CALL LJPSHIFT_UPDATE_PAIRG(X, V, EPSAA, SIGAA6, RCONSTAA, &
       CONSTAA, IRCUT2AA, POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, NATOMS,&
       PERIODIC )
          ENDDO
        ENDDO
        DO J1=1,NTYPEA
          DO J2=NTYPEA+1,N
            CALL LJPSHIFT_UPDATE_PAIRG(X, V, EPSAB, SIGAB6, RCONSTAB,&
       CONSTAB, IRCUT2AB, POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, NATOMS,&
       PERIODIC )
          ENDDO
        ENDDO
        DO J1=NTYPEA+1,N
          DO J2=J1+1,N
            CALL LJPSHIFT_UPDATE_PAIRG(X, V, EPSBB, SIGBB6, RCONSTBB, &
       CONSTBB, IRCUT2BB, POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, NATOMS,&
       PERIODIC )
          ENDDO
        ENDDO

      ELSE
         !only update POTEL
         DO J1=1,NTYPEA
            DO J2=J1+1,NTYPEA
               CALL LJPSHIFT_UPDATE_PAIR(X, EPSAA, SIGAA6, RCONSTAA, &
      CONSTAA, IRCUT2AA, POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, NATOMS,&
       PERIODIC )
            ENDDO
         ENDDO
         DO J1=1,NTYPEA
            DO J2=NTYPEA+1, N
               CALL LJPSHIFT_UPDATE_PAIR(X, EPSAB, SIGAB6, RCONSTAB, &
      CONSTAB, IRCUT2AB, POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, NATOMS,&
       PERIODIC )
            ENDDO
         ENDDO
         DO J1=NTYPEA+1, N
            DO J2=J1+1, N
               CALL LJPSHIFT_UPDATE_PAIR(X, EPSBB, SIGBB6, RCONSTBB, &
      CONSTBB, IRCUT2BB, POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, NATOMS,&
       PERIODIC )
            ENDDO
         ENDDO
      ENDIF
      !POTEL = POTEL * 4.D0

      RETURN
      END

!*******************************************************************

!*******************************************************************

      SUBROUTINE LJPSHIFT_UPDATE_PAIRG(X, V, EPSG, SIGG6, &
       RCONSTG, CONSTG, &
       IRCUT2G, POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, NATOMS, PERIODIC )
        !Calculate the potential energy and gradient between atoms j1 and j2.
        !Update V, and POTEL.
        !This routine should be equivalent to LJPSHIFT_UPDATE_PAIR except for
        !the update of V


        IMPLICIT NONE

        LOGICAL, INTENT(IN) :: PERIODIC
        INTEGER, INTENT(IN) :: J1, J2, NATOMS
        DOUBLE PRECISION, INTENT(IN) :: EPSG, SIGG6, RCONSTG, CONSTG
        DOUBLE PRECISION, INTENT(IN) :: BOXLX, BOXLY, BOXLZ, IRCUT2G
        DOUBLE PRECISION :: XVEC(3)
        DOUBLE PRECISION :: G
        DOUBLE PRECISION, INTENT(INOUT) :: V(3*NATOMS)
        DOUBLE PRECISION, INTENT(IN) :: X(3*NATOMS)
        DOUBLE PRECISION :: R2
        DOUBLE PRECISION :: R8
        DOUBLE PRECISION :: R14
        DOUBLE PRECISION, INTENT(INOUT) :: POTEL

        DOUBLE PRECISION :: R6
        INTEGER :: J3, J4, J5

        J3=3*(J1-1)
        J4=3*(J2-1)

        !calculate atom separation
        XVEC(1)=X(J3+1)-X(J4+1)
        XVEC(2)=X(J3+2)-X(J4+2)
        XVEC(3)=X(J3+3)-X(J4+3)
        IF (PERIODIC) THEN
           XVEC(1)=XVEC(1)-BOXLX*NINT(XVEC(1)/BOXLX)
           XVEC(2)=XVEC(2)-BOXLY*NINT(XVEC(2)/BOXLX)
           XVEC(3)=XVEC(3)-BOXLZ*NINT(XVEC(3)/BOXLX)
        ENDIF


        R2=1.0D0/(XVEC(1)**2+XVEC(2)**2+XVEC(3)**2)
        IF (R2.GT.IRCUT2G) THEN

          R6=R2**3
          R8=R6*R2
          R14=R8*R6

          !update potential energy
          POTEL = POTEL + 4.D0*EPSG*(SIGG6*R6*(SIGG6*R6-1.0D0) + &
       RCONSTG/R2 + CONSTG )

          !update derivative of potential energy
          G=-8.0D0*EPSG*(3.0D0*(2.0D0*R14*(SIGG6*SIGG6)&
       -R8*SIGG6)-RCONSTG)
          DO J5=1,3
            V(J3+J5)=V(J3+J5)+G*XVEC(J5)
            V(J4+J5)=V(J4+J5)-G*XVEC(J5)
          END DO

        ENDIF

      END SUBROUTINE LJPSHIFT_UPDATE_PAIRG

!*******************************************************************


!*******************************************************************

      SUBROUTINE LJPSHIFT_UPDATE_PAIR(X, EPSG, SIGG6, RCONSTG, &
       CONSTG, IRCUT2G, &
         POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, NATOMS, PERIODIC )
        !calculate the potential energy between atoms j1 and j2 
        !X are the atom positions, must apply BC's as needed


        IMPLICIT NONE

        LOGICAL, INTENT(IN) :: PERIODIC
        INTEGER, INTENT(IN) :: J1, J2, NATOMS
        DOUBLE PRECISION, INTENT(IN) :: EPSG, SIGG6, RCONSTG, CONSTG
        DOUBLE PRECISION, INTENT(IN) :: boxlx, boxly, boxlz, IRCUT2G
        DOUBLE PRECISION, INTENT(IN) :: X(3*NATOMS)
        DOUBLE PRECISION, INTENT(INOUT) :: POTEL

        DOUBLE PRECISION :: R2DUM, R6, VEC1, VEC2, VEC3
        INTEGER :: J3, J4

!       POTELDUM = 0.D0;

        J3=3*(J1-1)
        J4=3*(J2-1)
        !find particle separation with suitable boundary conditions
        VEC1=X(J3+1)-X(J4+1)
        VEC2=X(J3+2)-X(J4+2)
        VEC3=X(J3+3)-X(J4+3)
        IF (PERIODIC) THEN
           VEC1=VEC1-BOXLX*NINT(VEC1/BOXLZ)
           VEC2=VEC2-BOXLY*NINT(VEC2/BOXLZ)
           VEC3=VEC3-BOXLZ*NINT(VEC3/BOXLZ)
        ENDIF
        !calculate the potential
        R2DUM=VEC1**2+VEC2**2+VEC3**2
        R2DUM=1.0D0/R2DUM
        R6=R2DUM**3
        IF (R2DUM.GT.IRCUT2G) THEN
          POTEL = POTEL + 4.D0*EPSG*(SIGG6*R6*(SIGG6*R6-1.0D0) &
       + RCONSTG/R2DUM + CONSTG)  ! AB
        ENDIF

      END SUBROUTINE LJPSHIFT_UPDATE_PAIR

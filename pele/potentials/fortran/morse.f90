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
!  Subroutine MORSE calculates the cartesian gradient analytically for 
!  the Morse potential.
!
!*************************************************************************
!
      SUBROUTINE MORSE(X,V,EMORSE,GTEST, natoms, rho, R0, A)
      ! R0 is the position of the bottom of the well
      ! rho is the width of the well and has units of inverse length
      ! A is the energy scale
!      USE commons
      IMPLICIT NONE 
      LOGICAL GTEST
      INTEGER J1, J2, J3, J4
      integer, intent(IN) :: NATOMS
      DOUBLE PRECISION, intent(IN) :: X(3*NATOMS), rho, R0, A
      DOUBLE PRECISION, intent(OUT) :: V(3*NATOMS), EMORSE
      DOUBLE PRECISION DIST, R, DUMMY, &
                       RR(NATOMS,NATOMS), DUMMYX,  &
                       DUMMYY, DUMMYZ, XMUL2
!     LOGICAL EVAP, evapreject
!     COMMON /EV/ EVAP, evapreject

!     EVAP=.FALSE.
      V(:) = 0.D0
      EMORSE=0.0D0
      IF (GTEST) THEN
         DO J1=1,NATOMS
            J3=3*J1
            !DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
!           IF (DIST.GT.RADIUS) THEN
!              EVAP=.TRUE.
!              EMORSE=EMORSE+(DIST-RADIUS)**2
!           ENDIF
            RR(J1,J1)=0.0D0
            DO J2=J1+1,NATOMS
               J4=3*J2
               DIST=MAX(DSQRT((X(J3-2)-X(J4-2))**2 + (X(J3-1)-X(J4-1))**2 &
                        + (X(J3)-X(J4))**2),1.0D-5)
!              DIST=DSQRT((X(J3-2)-X(J4-2))**2 + (X(J3-1)-X(J4-1))**2 
!    1                  + (X(J3)-X(J4))**2)
               R=DEXP(RHO*R0-RHO*DIST)
               RR(J2,J1)=2.0D0*R*(R-1.0D0)/DIST
               RR(J1,J2)=RR(J2,J1)
               DUMMY=R*(R-2.0D0)
               EMORSE=EMORSE+DUMMY
            ENDDO
         ENDDO
      ELSE
         DO J1=1,NATOMS
            J3=3*J1
            DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
!           IF (DIST.GT.RADIUS) THEN
!              EVAP=.TRUE.
!              EMORSE=EMORSE+(DIST-RADIUS)**2
!           ENDIF
            DO J2=J1+1,NATOMS
               J4=3*J2
               DIST=DSQRT((X(J3-2)-X(J4-2))**2 + (X(J3-1)-X(J4-1))**2 &
                        + (X(J3)-X(J4))**2)
               R=DEXP(RHO*R0-RHO*DIST)
               DUMMY=R*(R-2.0D0)
               EMORSE=EMORSE+DUMMY
            ENDDO
         ENDDO
      ENDIF
      EMORSE = EMORSE * A

      IF (.NOT.GTEST) RETURN
 
      DO J1=1,NATOMS
         J3=3*J1
         DUMMYX=0.0D0
         DUMMYY=0.0D0
         DUMMYZ=0.0D0
         DO J4=1,NATOMS
            J2=3*J4
            XMUL2=RR(J4,J1)
            DUMMYX=DUMMYX-(X(J3-2)-X(J2-2))*XMUL2
            DUMMYY=DUMMYY-(X(J3-1)-X(J2-1))*XMUL2
            DUMMYZ=DUMMYZ-(X(J3)-X(J2))*XMUL2
         ENDDO
         DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
         V(J3-2)=DUMMYX
         V(J3-1)=DUMMYY
         V(J3)=DUMMYZ
      ENDDO

      V(:) = V(:) * A

      RETURN
      END

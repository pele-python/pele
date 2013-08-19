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
      SUBROUTINE MORSE(X,V,EMORSE,GTEST, natoms, rho, R0, A, periodic, boxvec)
      ! R0 is the position of the bottom of the well
      ! rho is the width of the well and has units of inverse length
      ! A is the energy scale
!      USE commons
      IMPLICIT NONE 
      LOGICAL, intent(IN) :: GTEST, periodic
      integer, intent(IN) :: NATOMS
      DOUBLE PRECISION, intent(IN) :: X(3*NATOMS), rho, R0, A, boxvec(3)
      DOUBLE PRECISION, intent(OUT) :: V(3*NATOMS), EMORSE
      INTEGER J1, J2, J3, J4
      DOUBLE PRECISION DIST, R, DUMMY, &
                       RR(NATOMS,NATOMS), &
                       XMUL2, iboxvec(3), dx(3)
!     LOGICAL EVAP, evapreject
!     COMMON /EV/ EVAP, evapreject
      if (periodic) iboxvec(:) = 1.d0 / boxvec(:)

!     EVAP=.FALSE.
      V(:) = 0.D0
      EMORSE=0.0D0
      DO J1=1,NATOMS
         J3=3*J1
         RR(J1,J1)=0.0D0
         DO J2=J1+1,NATOMS
            J4=3*J2
            dx(:) = X(J3-2:j3)-X(J4-2:j4)
            if (periodic) then
               dx = dx - boxvec * nint(dx * iboxvec)
            endif
            dist = max(sqrt(sum(dx**2)), 1.d-5)
            R=DEXP(RHO*R0-RHO*DIST)
            DUMMY=R*(R-2.0D0)
            EMORSE=EMORSE+DUMMY

            if (gtest) then
               xmul2 = 2.0D0*R*(R-1.0D0)/DIST * A
               V(J3-2:j2) = V(j3-2:j2) - xmul2 * dx
               V(J4-2:j2) = V(j4-2:j2) + xmul2 * dx
            endif
         ENDDO
      ENDDO
      EMORSE = EMORSE * A

      RETURN
      END

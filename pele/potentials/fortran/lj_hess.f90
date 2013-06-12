!   OPTIM: A program for optimizing geometries and calculating reaction pathways
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of OPTIM.
!
!   OPTIM is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   OPTIM is distributed in the hope that it will be useful,
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
!  Subroutine LJDIFF calculates the cartesian gradient and second
!  derivative matrix analytically. Reduced units.
!
!  This method stores the distances in large arrays, so it is not optimal for
!  large systems
!
!*************************************************************************
!
      SUBROUTINE LJDIFF(N, X, V, ENERGY, GTEST, STEST, HESS)
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  N
      LOGICAL, INTENT(IN) ::  GTEST, STEST
      DOUBLE PRECISION, INTENT(IN) :: X(3*N)
      DOUBLE PRECISION, INTENT(OUT) :: ENERGY, V(3*N), HESS(3*N,3*N)
      INTEGER J1, J2, J3, J4
      DOUBLE PRECISION R6,&
                       R2(N,N), R2T,&
                       R8(N,N), G(N,N), XG(N,N),&
                       R14(N,N), F(N,N), DUMMY, DUMMYX, DUMMYY, DUMMYZ, DIST, XMUL2
! 
!  Store distance matrices.
!
      ENERGY=0.0D0
      IF (GTEST.AND.(.NOT.STEST)) THEN
         DO J1=1,N
            J3=3*J1
            XG(J1,J1)=0.0D0
            DO J2=J1+1,N
               J4=3*J2
               DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
               DIST=1.0D0/DIST
               R6=DIST**3
               DUMMY=R6*(R6-1.0D0)
               ENERGY=ENERGY+DUMMY
               DIST=DIST*R6
               XG(J2,J1)=-24.0D0*(2.0D0*R6-1.0D0)*DIST
               XG(J1,J2)=XG(J2,J1)
            ENDDO
         ENDDO
      ELSEIF (GTEST) THEN
         DO J1=1,N
            R2(J1,J1)=0.0D0
            R8(J1,J1)=0.0D0
            R14(J1,J1)=0.0D0
            DO J2=J1+1,N
               R2(J2,J1)=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2&
                        +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2&
                        +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
               R2(J2,J1)=1.0D0/R2(J2,J1)
               R6=R2(J2,J1)**3
               ENERGY=ENERGY+R6*(R6-1.0D0)
               R8(J2,J1)=R2(J2,J1)**4
               R14(J2,J1)=R8(J2,J1)*R8(J2,J1)/R2(J2,J1)
               R2(J1,J2)=R2(J2,J1)
               XG(J2,J1)=-24.0D0*(2.0D0*R6-1.0D0)*R2(J1,J2)*R6
               XG(J1,J2)=XG(J2,J1)
            ENDDO
         ENDDO 
      ELSE
         DO J1=1,N
            J3=3*(J1-1)
            DO J2=J1+1,N
               J4=3*(J2-1)
               R2T=(X(J3+1)-X(J4+1))**2+(X(J3+2)-X(J4+2))**2+(X(J3+3)-X(J4+3))**2
               R2T=1.0D0/R2T
               R6=R2T**3
               ENERGY=ENERGY+R6*(R6-1.0D0)
            ENDDO
         ENDDO

      ENDIF
      ENERGY=4.0D0*ENERGY

      IF (.NOT.GTEST) RETURN
!     CALL LJG(G,R14,R8,V,X,N)
      DO J1=1,N
         J3=3*J1
         DUMMYX=0.0D0
         DUMMYY=0.0D0
         DUMMYZ=0.0D0
         DO J4=1,N
            J2=3*J4
            XMUL2=XG(J4,J1)
            DUMMYX=DUMMYX+XMUL2*(X(J3-2)-X(J2-2))
            DUMMYY=DUMMYY+XMUL2*(X(J3-1)-X(J2-1))
            DUMMYZ=DUMMYZ+XMUL2*(X(J3)  -X(J2))
         ENDDO
         V(J3-2)=DUMMYX
         V(J3-1)=DUMMYY
         V(J3)=DUMMYZ
      ENDDO
      
      IF (.NOT.STEST) RETURN
      CALL LJS(G,F,R2,R14,R8,X,N,HESS)

      RETURN
      END

!*****************************************************************************

      SUBROUTINE LJS(G,F,R2,R14,R8,X,N,HESS)
      !USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6
      DOUBLE PRECISION, INTENT(OUT) :: HESS(3*N,3*N)
      DOUBLE PRECISION G(N,N), R14(N,N), R8(N,N),&
                       R2(N,N), F(N,N), &
                       X(3*N),DUMMY

!
!  Calculate the g tensor.
!
      DO J1=1,N
         G(J1,J1)=0.0D0
         DO J2=J1+1,N
            G(J2,J1)=-24.0D0*(2.0D0*R14(J2,J1)-R8(J2,J1))
            G(J1,J2)=G(J2,J1)
         ENDDO
      ENDDO

      DO J1=1,N
         F(J1,J1)=0.0D0
         DO J2=J1+1,N 
            F(J2,J1)=672.0D0*R14(J2,J1)-192.0D0*R8(J2,J1)
            F(J1,J2)=F(J2,J1)
         ENDDO
      ENDDO
!
!  Now do the hessian. First are the entirely diagonal terms.
!
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DUMMY=0.0D0
            DO J4=1,N
               DUMMY=DUMMY+F(J4,J1)*R2(J4,J1)*&
                       (X(J3)-X(3*(J4-1)+J2))**2 + G(J4,J1)   
            ENDDO
            HESS(J3,J3)=DUMMY
         ENDDO
      ENDDO
!
!  Next are the terms where x_i and x_j are on the same atom
!  but are different, e.g. y and z.
!
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J2+1,3
               DUMMY=0.0D0
               DO J5=1,N
                  DUMMY=DUMMY + F(J5,J1)*R2(J5,J1)* &
                 (X(J3)-X(3*(J5-1)+J2))*(X(3*(J1-1)+J4)-X(3*(J5-1)+J4)) 
               ENDDO
               HESS(3*(J1-1)+J4,J3)=DUMMY
            ENDDO
         ENDDO
      ENDDO
!
!  Case III, different atoms, same cartesian coordinate.
!
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,N
               HESS(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*&
                                 (X(J3)-X(3*(J4-1)+J2))**2-G(J4,J1) 
            ENDDO
         ENDDO
      ENDDO
!
!  Case IV: different atoms and different cartesian coordinates.
!
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,N
               DO J5=1,J2-1
                  J6=3*(J4-1)+J5
                  HESS(J6,J3)=-F(J4,J1)*R2(J4,J1)&
                          *(X(J3)-X(3*(J4-1)+J2))&
                          *(X(3*(J1-1)+J5)-X(J6))
                  HESS(3*(J4-1)+J2,3*(J1-1)+J5)=HESS(J6,J3)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
!  Symmetrise Hessian
!
      DO J1=1,3*N
         DO J2=J1+1,3*N
            HESS(J1,J2)=HESS(J2,J1)
         ENDDO
      ENDDO
      RETURN
      END


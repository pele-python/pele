!
! GPL License Info 
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
!***********************************************************************

!
!  Normalize vector VEC1
!
      SUBROUTINE VECNORM(VEC1,NOPT)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NOPT
      DOUBLE PRECISION, INTENT(INOUT) :: VEC1(NOPT)
      INTEGER J2
      DOUBLE PRECISION DUMMY2

      DUMMY2=0.0D0
      DO J2=1,NOPT
         DUMMY2=DUMMY2+VEC1(J2)**2
      ENDDO
      IF (DUMMY2.GT.0.0D0) THEN
         DUMMY2=1.0D0/DSQRT(DUMMY2)
         DO J2=1,NOPT
            VEC1(J2)=VEC1(J2)*DUMMY2
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE VECNORM

      SUBROUTINE ORTHOGOPT(VEC1,COORDS,OTEST, natoms, translation_only) 
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NATOMS
      DOUBLE PRECISION, INTENT(IN) :: COORDS(3*NATOMS)
      DOUBLE PRECISION, INTENT(INOUT) :: VEC1(3*NATOMS)
      LOGICAL, INTENT(IN) :: OTEST, translation_only
      INTEGER J2, J3, NCHECK, NOPT
      DOUBLE PRECISION DUMMY1, DUMMY2, ROOTN, VDOT, &
                      CMX, CMY, CMZ, AMASS(NATOMS), TMASS, RMASS(NATOMS)
      NOPT = 3*NATOMS
      !MASST = .FALSE.

      ROOTN=SQRT(1.0D0*NATOMS)
      NCHECK=0
      CMX=0.0D0
      CMY=0.0D0
      CMZ=0.0D0
      TMASS=0.0D0
      DO J2=1,NATOMS
         AMASS(J2)=1.0D0
         RMASS(J2)=1.0D0
         !IF (MASST) AMASS(J2)=ATMASS(J2)
         !IF (MASST) RMASS(J2)=SQRT(ATMASS(J2))
         TMASS=TMASS+AMASS(J2)
      ENDDO
!
!  If MASST then the coordinates already have a square root of the mass in them
!
      DO J2=1,NATOMS
         CMX=CMX+COORDS(3*(J2-1)+1)*RMASS(J2)
         CMY=CMY+COORDS(3*(J2-1)+2)*RMASS(J2)
         CMZ=CMZ+COORDS(3*(J2-1)+3)*RMASS(J2)
      ENDDO
      CMX=CMX/TMASS
      CMY=CMY/TMASS
      CMZ=CMZ/TMASS
!     WRITE(*,'(A,3F20.10)') 'centre at ',CMX,CMY,CMZ
1     VDOT=0.0D0
      NCHECK=NCHECK+1
!
      !IF ((VDOT.GT.1.0D-6).AND.(SHIFTED).AND.(NCHECK.LT.100)) GOTO 1
      IF (NCHECK.GE.100) THEN
         write(*,*) '*** WARNING, cannot orthogonalise to known eigenvectors in ORTHOGOPT'
      ENDIF

      !
      !orthogonalize to translations in the x direction
      !
      DUMMY1=0.0D0
      DO J2=1,NOPT,3
         DUMMY1=DUMMY1+VEC1(J2)
      ENDDO
      DUMMY1=DUMMY1/ROOTN
      VDOT=MAX(VDOT,ABS(DUMMY1))
!     WRITE(*,'(A,F20.10)') 'X dot:',DUMMY1
      DO J2=1,NOPT,3
         VEC1(J2)=VEC1(J2)-DUMMY1/ROOTN
      ENDDO
      IF (OTEST) CALL VECNORM(VEC1,NOPT)

      !
      !orthogonalize to translations in the y direction
      !
      DUMMY1=0.0D0
      DO J2=2,NOPT,3
         DUMMY1=DUMMY1+VEC1(J2)
      ENDDO
      DUMMY1=DUMMY1/ROOTN
      VDOT=MAX(VDOT,ABS(DUMMY1))
!     WRITE(*,'(A,F20.10)') 'Y dot:',DUMMY1
      DO J2=2,NOPT,3
         VEC1(J2)=VEC1(J2)-DUMMY1/ROOTN
      ENDDO
      IF (OTEST) CALL VECNORM(VEC1,NOPT)

      !
      !orthogonalize to translations in the z direction
      !
30    DUMMY1=0.0D0
      DO J2=3,NOPT,3
         DUMMY1=DUMMY1+VEC1(J2)
      ENDDO
      DUMMY1=DUMMY1/ROOTN
      VDOT=MAX(VDOT,ABS(DUMMY1))
!     WRITE(*,'(A,F20.10)') 'Z dot:',DUMMY1
      DO J2=3,NOPT,3
         VEC1(J2)=VEC1(J2)-DUMMY1/ROOTN
      ENDDO
      IF (OTEST) CALL VECNORM(VEC1,NOPT)

31    CONTINUE
      if (translation_only) then
         IF ((VDOT.GT.1.0D-6).AND.(NCHECK.LT.100)) GOTO 1
         return
      endif
      !IF ((VDOT.GT.1.0D-6).AND.BULKT.AND.(NCHECK.LT.100)) GOTO 1
!     PRINT*,'after next part VDOT=',VDOT
      IF (NCHECK.GE.100) THEN
         PRINT*,'*** WARNING, cannot orthogonalise to known eigenvectors in ORTHOGOPT'
      ENDIF
      !IF (BULKT.AND.TWOD) GOTO 20
      !IF (BULKT) RETURN

10    DUMMY1=0.0D0
      DUMMY2=0.0D0
      DO J2=1,NATOMS
         J3=3*J2
         DUMMY1=DUMMY1+RMASS(J2)*(VEC1(J3-1)*(COORDS(J3)-CMZ)-VEC1(J3)*(COORDS(J3-1)-CMY))
         DUMMY2=DUMMY2+AMASS(J2)*((COORDS(J3)-CMZ)**2+(COORDS(J3-1)-CMY)**2)
      ENDDO
      IF (DUMMY2.GT.0.0D0) THEN
!        WRITE(*,'(A,F20.10)') 'XRdot:',DUMMY1/SQRT(DUMMY2)
         VDOT=MAX(VDOT,ABS(DUMMY1)/SQRT(DUMMY2))
         DUMMY2=DUMMY1/DUMMY2
         DO J2=1,NATOMS
            J3=3*J2
            VEC1(J3-1)=VEC1(J3-1)-DUMMY2*RMASS(J2)*(COORDS(J3)-CMZ)
            VEC1(J3)=VEC1(J3)+DUMMY2*RMASS(J2)*(COORDS(J3-1)-CMY)
         ENDDO
      ENDIF
      IF (OTEST) CALL VECNORM(VEC1,NOPT)
!     PRINT*,'rot x component=',DUMMY2

      DUMMY1=0.0D0
      DUMMY2=0.0D0
      DO J2=1,NATOMS
         J3=3*J2
         DUMMY1=DUMMY1-VEC1(J3-2)*RMASS(J2)*(COORDS(J3)-CMZ)+VEC1(J3)*RMASS(J2)*(COORDS(J3-2)-CMX)
         DUMMY2=DUMMY2+AMASS(J2)*((COORDS(J3)-CMZ)**2+(COORDS(J3-2)-CMX)**2)
      ENDDO
      IF (DUMMY2.GT.0.0D0) THEN
!        WRITE(*,'(A,F20.10)') 'YRdot:',DUMMY1/SQRT(DUMMY2)
         VDOT=MAX(VDOT,ABS(DUMMY1)/SQRT(DUMMY2))
         DUMMY2=DUMMY1/DUMMY2
         DO J2=1,NATOMS
            J3=3*J2
            VEC1(J3-2)=VEC1(J3-2)+DUMMY2*RMASS(J2)*(COORDS(J3)-CMZ)
            VEC1(J3)=VEC1(J3)-DUMMY2*RMASS(J2)*(COORDS(J3-2)-CMX)
         ENDDO
      ENDIF
      IF (OTEST) CALL VECNORM(VEC1,NOPT)
!     PRINT*,'rot y component=',DUMMY2

20    DUMMY1=0.0D0
      DUMMY2=0.0D0
      DO J2=1,NATOMS
         J3=3*J2
         DUMMY1=DUMMY1+VEC1(J3-2)*RMASS(J2)*(COORDS(J3-1)-CMY)-VEC1(J3-1)*RMASS(J2)*(COORDS(J3-2)-CMX)
         DUMMY2=DUMMY2+AMASS(J2)*((COORDS(J3-1)-CMY)**2+(COORDS(J3-2)-CMX)**2)
      ENDDO
      IF (DUMMY2.GT.0.0D0) THEN
!        WRITE(*,'(A,F20.10)') 'ZRdot:',DUMMY1/SQRT(DUMMY2)
         VDOT=MAX(VDOT,ABS(DUMMY1)/SQRT(DUMMY2))
         DUMMY2=DUMMY1/DUMMY2
         DO J2=1,NATOMS
            J3=3*J2
            VEC1(J3-2)=VEC1(J3-2)-DUMMY2*RMASS(J2)*(COORDS(J3-1)-CMY)
            VEC1(J3-1)=VEC1(J3-1)+DUMMY2*RMASS(J2)*(COORDS(J3-2)-CMX)
         ENDDO
      ENDIF
      IF (OTEST) CALL VECNORM(VEC1,NOPT)
!     PRINT*,'rot z component=',DUMMY2

!     WRITE(*,'(A,F20.10)') 'Largest remaining component in ORTHOGOPT=',VDOT
      IF ((VDOT.GT.1.0D-6).AND.(NCHECK.LT.100)) GOTO 1
!     PRINT*,'after next part VDOT=',VDOT
      IF (NCHECK.GE.100) THEN
         PRINT*,'*** WARNING, cannot orthogonalise to known eigenvectors in ORTHOGOPT'
      ENDIF
      RETURN
      END SUBROUTINE ORTHOGOPT

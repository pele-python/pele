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
!********************************************************************
!
! compute the energy and gradient of the sutton chen potential.  Add a cutoff
! which goes smoothly to zero.  The cutoff is implemented by replacing (1/r)^k
! with
!
! (1./r)**k - (1./rc)**k + k*(r-rc)*(1./rc)**(k+1) - 0.5*k*(k+1.)*(r-rc)**2*(1/rc)**(k+2)
!
!********************************************************************
!
      SUBROUTINE SCDIFF_periodic(N,X,V,EPS,C,SIG,NN,MM,PSC, boxvec, rcut)
      !USE MODHESS
      !use commons, only : param6, param7, param8
      IMPLICIT NONE
      integer, intent(IN) :: N, NN, MM
      double precision, intent(IN) :: x(3*N)
      double precision, intent(IN) :: eps, c, sig, rcut, boxvec(3)
      double precision, intent(out) :: psc ! energy
      double precision, intent(out) :: v(3*N) !gradient
      INTEGER J1, J2, J3, J4, I, J
      DOUBLE PRECISION RHO(N), RM(N,N), &
                       R1(N,N), &
                       RN(N,N),RN2(N,N), &
                       RM2(N,N), &
                       DIST,POTA,POTB, R
      double precision iboxvec(3), dx(3), dr
      double precision rho1, ircm, ircn, sigmm, signn
      IBOXVEC = 1.D0 / BOXVEC

      ! compute the parameters for the cutoff smoothing function
      IRCM = 1./RCUT**MM
      IRCN = 1./RCUT**NN
      SIGMM = SIG**MM
      sigNN = sig**NN
!
!
!  Store distance matrices.
!
      DO J1=1,N
         RM(J1,J1)=0.0D0
         RN(J1,J1)=0.0D0
         RM2(J1,J1)=0.0D0
         RN2(J1,J1)=0.0D0
         DO J2=1,J1-1
            dx(1) = (X(3*(J1-1)+1)-X(3*(J2-1)+1))
            dx(2) = (X(3*(J1-1)+2)-X(3*(J2-1)+2))
            dx(3) = (X(3*(J1-1)+3)-X(3*(J2-1)+3))
            dx = dx - boxvec * nint(dx * iboxvec)
            DIST = sum(dx**2)
            R1(J2,J1)=sqrt(DIST)
            RM(J2,J1)=DIST**(-FLOAT(MM)/2.0D0)
            RN(J2,J1)=DIST**(-FLOAT(NN)/2.0D0)
            RM2(J2,J1)=DIST**(-(FLOAT(MM)+2.0D0)/2.0D0)
            RN2(J2,J1)=DIST**(-(FLOAT(NN)+2.0D0)/2.0D0)
            R1(J1,J2)=R1(J2,J1)
            RM(J1,J2)=RM(J2,J1)
            RN(J1,J2)=RN(J2,J1)
            RM2(J1,J2)=RM2(J2,J1)
            RN2(J1,J2)=RN2(J2,J1)
         end do
      end do
!
! Store density matrix: In the case of the perfect fcc lattice,
! the infinitely extended crystal implies that every RHO(J) is
! equal to RHO(1).
!
      DO I=1,N
         RHO(I)=0.0D00
         DO J=1,N
            if (i.ne.j) then
               R = R1(I,J)
               IF (R .LT. RCUT) THEN
                  RHO1 = SIGMM*(RM(I,J) + &
                     IRCM*(-1.D0 + MM*(R-RCUT)/RCUT - 0.5*MM*(MM+1)*(R-RCUT)**2/RCUT**2))
                  RHO(I)=RHO(I) + rho1
               ENDIF
            endif
         end do
      end do
!
! First calculate the potential energy:
!
      POTA=0.0D0
      POTB=0.0D0
      DO I=1,N
         DO J=1,N
            if (i.ne.j) then
               IF (R1(I,J) .LT. RCUT) THEN
                  POTA=POTA + 0.50D00*EPS*SIGNN*(RN(I,J) + &
                     IRCN*(-1.D0 + NN*(R-RCUT)/RCUT - 0.5*NN*(NN+1)*(R-RCUT)**2/RCUT**2))
               endif
            endif
         end do
         POTB=POTB + EPS*DSQRT(RHO(I))*C
      end do
      PSC=POTA - POTB
      !PRINT*, ' Sutton-Chen n=', NN,' m=', MM
      !PRINT*,'Two-body contribution=',POTA,' eV'
      !PRINT*,'Many-body contribution=',-POTB,' eV'
!     PRINT*, 'Total Energy for last step=', PSC,' eV'
!
! Now calculate the gradient analytically.
!
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            V(J3)=0.0D0
            DO J4=1,N
               if (j1.ne.j4) then
                  IF (R1(j1,j4) .LT. RCUT) THEN
                     dr = X(J3)-X(3*(J4-1)+J2)
                     dr = dr - boxvec(j2) * nint(dr * iboxvec(j2))
                     V(J3)=V(J3)+(-NN*EPS*SIGNN &
                                  *(RN2(J1,J4) - IRCN*(NN/RCUT - NN*(NN+1)*(R-RCUT)/RCUT**2 ) ) &
                                  + &
                                  FLOAT(MM)/2.d0*EPS*C*SIGMM &
                                  *(RM2(J1,J4) - IRCM*(MM/RCUT - MM*(MM+1)*(R-RCUT)/RCUT**2)) &
                                  *(1.d0/DSQRT(RHO(J1)) + 1.d0/DSQRT(RHO(J4))) &
                                 )*dr
                  endif
              endif
            end do
!           PRINT*,'V(',J3,')=',V(J3)
         end do
      end do
      RETURN
      END SUBROUTINE SCDIFF_periodic

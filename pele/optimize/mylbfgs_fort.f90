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
!        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION
!                          JORGE NOCEDAL
!                        *** July 1990 ***
!
!        Line search removed plus small modifications, DJW 2001
!
      SUBROUTINE MYLBFGS(N, M, XCOORDS, EPS, MFLAG, ENERGY, &
         ITMAX, ITDONE, NPOTCALLS, MAXBFGS, MAXERISE, GRAD)

      IMPLICIT NONE
      INTEGER, INTENT(IN) ::             N
      INTEGER, INTENT(IN) ::             M
      DOUBLE PRECISION                       XCOORDS(N)
!f2py intent(in,out) XCOORDS(N)
      DOUBLE PRECISION, INTENT(IN) ::    EPS
      LOGICAL, INTENT(OUT) ::            MFLAG
      DOUBLE PRECISION, INTENT(OUT) ::   ENERGY
      INTEGER, INTENT(IN) ::             ITMAX
      INTEGER, INTENT(OUT) ::            ITDONE
      INTEGER, INTENT(OUT) ::            NPOTCALLS
!f2py INTENT(CALLBACK) POTENTIAL
      EXTERNAL POTENTIAL

      DOUBLE PRECISION, INTENT(IN) ::    MAXBFGS !max step size
      DOUBLE PRECISION             ::    MAXEFALL
      DOUBLE PRECISION, INTENT(IN) ::    MAXERISE
      DOUBLE PRECISION, INTENT(OUT) ::   GRAD(N)
!
!     variables that need to be saved from step to step
!
      INTEGER                            ITER
      INTEGER                            POINT
      INTEGER                            ISPT
      INTEGER                            IYPT
      INTEGER                            NPT
      DOUBLE PRECISION                   W(N*(2*M+1)+2*M)
      LOGICAL                            DIAGCO
      DOUBLE PRECISION                   DIAG(N)
!
!     input from the potential
!
      DOUBLE PRECISION                   XSAVE(N)
      DOUBLE PRECISION                   RMS
      DOUBLE PRECISION                   DGUESS
      DOUBLE PRECISION                   GNEW(N)
      DOUBLE PRECISION                   ENEW
!
!
      !INTEGER, INTENT(INOUT) ::          BOUND
      INTEGER     BOUND
      DOUBLE PRECISION                   STP
      INTEGER                            NFAIL
      INTEGER                            NDECREASE
!
!
!     finally the temporary variables
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCcc
      !LOGICAL FIXIMAGE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER J1 !definitely local variables
      DOUBLE PRECISION SLENGTH, DUMMY, WTEMP(N)
      DOUBLE PRECISION OVERLAP !definitely local
      DOUBLE PRECISION GNORM, YS, YY, SQ, YR, BETA !definitely local
      DOUBLE PRECISION DOT1, DOT2
      INTEGER CP, INMC, IYCN, ISCN !definintly local
      !
      !function definitions
      !
      DOUBLE PRECISION DDOT
      !SAVE W, DIAG, ITER, POINT, ISPT, IYPT, NPT

240         FORMAT(' IMPROPER INPUT PARAMETERS (N OR M ARE NOT POSITIVE)')
235               FORMAT(' THE',I5,'-TH DIAGONAL ELEMENT OF THE INVERSE HESSIAN APPROXIMATION IS NOT POSITIVE')

      NPOTCALLS = 0
      DIAGCO  = .FALSE.
      !write(*,*) "starting minimizer", n, m
      !write(*,*) "fort n=", n
      !write(*,*) "fort m=", m
      !CALL FLUSH()

      DGUESS=1.0d0
      !MAXBFGS=0.5d0 !input
      MAXEFALL=-1.D100
      RMS=1.0d0



      NFAIL=0
      !IF (RESET) ITER=0
      ITER = 0
      ITDONE=0
      !FIXIMAGE=.FALSE.

      !write(*,*) "calling potential"
      !CALL FLUSH()
      CALL POTENTIAL(XCOORDS,GRAD,ENERGY, N)
      !write(*,*) "done calling potential", energy
      !write(*,*) grad
      !CALL FLUSH()
      NPOTCALLS = NPOTCALLS + 1
      RMS=DOT_PRODUCT(GRAD,GRAD)

!
!  Termination test. 
!
! 10    CALL FLUSH(*)
10    MFLAG=.FALSE.
      IF (RMS.LE.EPS) THEN ! vr274> added cell change test for DMACRYS
        MFLAG=.TRUE.
        RETURN
      ENDIF
      

      IF (ITDONE.EQ.ITMAX) THEN
         RETURN
      ENDIF


      IF (ITER.EQ.0) THEN
         IF (N.LE.0.OR.M.LE.0) THEN
            WRITE(*,240)
            STOP
         ENDIF
         POINT=0
         MFLAG=.FALSE.
         IF (DIAGCO) THEN
            WRITE(*,'(A)') 'using estimate of the inverse diagonal elements'
            DO J1=1,N
               IF (DIAG(J1).LE.0.0D0) THEN
                  WRITE(*,235) J1
                  STOP
               ENDIF
            ENDDO
         ELSE
            DO J1=1,N
               DIAG(J1)=DGUESS
            ENDDO
         ENDIF

!
!     THE WORK VECTOR W IS DIVIDED AS FOLLOWS:
!     ---------------------------------------
!     THE FIRST N LOCATIONS ARE USED TO STORE THE GRADIENT AND
!         OTHER TEMPORARY INFORMATION.
!     LOCATIONS (N+1)...(N+M) STORE THE SCALARS RHO.
!     LOCATIONS (N+M+1)...(N+2M) STORE THE NUMBERS ALPHA USED
!         IN THE FORMULA THAT COMPUTES H*G.
!     LOCATIONS (N+2M+1)...(N+2M+NM) STORE THE LAST M SEARCH
!         STEPS.
!     LOCATIONS (N+2M+NM+1)...(N+2M+2NM) STORE THE LAST M
!         GRADIENT DIFFERENCES.
!
!     THE SEARCH STEPS AND GRADIENT DIFFERENCES ARE STORED IN A
!     CIRCULAR ORDER CONTROLLED BY THE PARAMETER POINT.
!
         ISPT= N+2*M    ! index for storage of search steps
         IYPT= ISPT+N*M ! index for storage of gradient differences
!
!  NR step for diagonal inverse Hessian
!
         DO J1=1,N
            DUMMY=-GRAD(J1)*DIAG(J1)
            W(ISPT+J1)=DUMMY
            W(J1)=DUMMY
         ENDDO
         GNORM=DSQRT(DDOT(N,GRAD,1,GRAD,1))
!
!  Make the first guess for the step length cautious.
!
         STP=MIN(1.0D0/GNORM,GNORM)
      ELSE 
         BOUND=ITER
         IF (ITER.GT.M) BOUND=M
         YS= DDOT(N,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)
!
!  Update estimate of diagonal inverse Hessian elements
!
         IF (.NOT.DIAGCO) THEN
            YY= DDOT(N,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)
            IF (YY.EQ.0.0D0) THEN
               WRITE(*,'(A)') 'WARNING, resetting YY to one in mylbfgs'
               YY=1.0D0
            ENDIF
            IF (YS.EQ.0.0D0) THEN
               WRITE(*,'(A)') 'WARNING, resetting YS to one in mylbfgs'
               YS=1.0D0
            ENDIF
!           WRITE(*,'(A,2F20.10)') 'YS/YY,STP=',YS/YY,STP
            DO J1=1,N
!              DIAG(J1)= ABS(YS/YY) ! messes up after step reversals!
               DIAG(J1)= YS/YY
            ENDDO
         ELSE
            WRITE(*,'(A)') 'using estimate of the inverse diagonal elements'
            DO J1=1,N
               IF (DIAG(J1).LE.0.0D0) THEN
                  WRITE(*,235) J1
                  STOP
               ENDIF
            ENDDO
         ENDIF
!
!     COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
!     "Updating quasi-Newton matrices with limited storage",
!     Mathematics of Computation, Vol.24, No.151, pp. 773-782.
!     ---------------------------------------------------------
!

         CP= POINT
         IF (POINT.EQ.0) CP=M
         W(N+CP)= 1.0D0/YS
         DO J1=1,N
            W(J1)= -GRAD(J1)
         ENDDO
         CP= POINT
         DO J1= 1,BOUND
            CP=CP-1
            IF (CP.EQ.-1) CP=M-1
            SQ=DDOT(N,W(ISPT+CP*N+1),1,W,1)
            INMC=N+M+CP+1
            IYCN=IYPT+CP*N
            W(INMC)=W(N+CP+1)*SQ
            CALL DAXPY(N,-W(INMC),W(IYCN+1),1,W,1)
         ENDDO
        
         DO J1=1,N
            W(J1)=DIAG(J1)*W(J1)
         ENDDO

         DO J1=1,BOUND
            YR= DDOT(N,W(IYPT+CP*N+1),1,W,1)
            BETA= W(N+CP+1)*YR
            INMC=N+M+CP+1
            BETA= W(INMC)-BETA
            ISCN=ISPT+CP*N
            CALL DAXPY(N,BETA,W(ISCN+1),1,W,1)
            CP=CP+1
            IF (CP.EQ.M) CP=0
         ENDDO
         STP=1.0D0  
      ENDIF
!
!  Store the new search direction
!
      IF (ITER.GT.0) THEN
         DO J1=1,N
            W(ISPT+POINT*N+J1)= W(J1)
         ENDDO
      ENDIF

      DOT1=SQRT(DDOT(N,GRAD,1,GRAD,1))
!
!  Overflow has occasionally occurred here.
!  We only need the sign of the overlap, so use a temporary array with
!  reduced elements.
!
      DUMMY=1.0D0
      DO J1=1,N
         IF (ABS(W(J1)).GT.DUMMY) DUMMY=ABS(W(J1))
      ENDDO
      DO J1=1,N
         WTEMP(J1)=W(J1)/DUMMY
      ENDDO
      DOT2=SQRT(DDOT(N,WTEMP,1,WTEMP,1))
      OVERLAP=0.0D0
      IF (DOT1*DOT2.NE.0.0D0) THEN
        OVERLAP=DDOT(N,GRAD,1,WTEMP,1)/(DOT1*DOT2)
      ENDIF
      IF (OVERLAP.GT.0.0D0) THEN
         DO J1=1,N
            W(ISPT+POINT*N+J1)= -W(J1)  !!! DJW, reverses step
         ENDDO
      ENDIF

      DO J1=1,N
         W(J1)=GRAD(J1)
      ENDDO
      SLENGTH=0.0D0
      DO J1=1,N
         SLENGTH=SLENGTH+W(ISPT+POINT*N+J1)**2
      ENDDO
      SLENGTH=SQRT(SLENGTH)
      IF (STP*SLENGTH.GT.MAXBFGS) STP=MAXBFGS/SLENGTH


!
!  We now have the proposed step.
!!
! Save XCOORDS here so that we can undo the step reliably including the
! non-linear projection for Thomson for the angular coordinates.
!
      XSAVE(1:N)=XCOORDS(1:N) 
      DO J1=1,N
         XCOORDS(J1)=XCOORDS(J1)+STP*W(ISPT+POINT*N+J1)
      ENDDO 
      NDECREASE=0
20    CALL POTENTIAL(XCOORDS,GNEW,ENEW, N)
      NPOTCALLS = NPOTCALLS + 1
      !js850> must do something with ENEW / ENERGY here
      RMS=DOT_PRODUCT(GRAD,GRAD)
!
      IF (((ENEW-ENERGY.LE.MAXERISE)).AND.(ENEW-ENERGY.GT.MAXEFALL)) THEN
         ITER=ITER+1
         ITDONE=ITDONE+1
         ENERGY=ENEW
         DO J1=1,N
            GRAD(J1)=GNEW(J1)
         ENDDO
!
!  Step finished so can reset OLDQ to new XINT, OLDCART to new LCART,
!  as well as the Cartesian and internal gradients.
!
      ELSEIF (ENEW-ENERGY.LE.MAXEFALL) THEN
!
!  Energy decreased too much - try again with a smaller step size
!
         IF (NDECREASE.GT.5) THEN
            NFAIL=NFAIL+1
            WRITE(*,'(A,G20.10)') ' in mylbfgs LBFGS step cannot find an energy in the required range, NFAIL=',NFAIL
!
! Resetting to XSAVE should be the same as subtracting the step. 
! If we have tried PROJI with Thomson then the projection is non-linear
! and we need to reset to XSAVE. This should always be reliable!
!
            XCOORDS(1:N)=XSAVE(1:N)
            GRAD(1:N)=GNEW(1:N) ! GRAD contains the gradient at the lowest energy point

            ITER=0   !  try resetting
            IF (NFAIL.GT.20) THEN
               WRITE(*,'(A)') ' Too many failures - giving up '
               !FIXIMAGE=.FALSE.
!              STOP
               RETURN
            ENDIF
            GOTO 30
         ENDIF
!
! Resetting to XSAVE and adding half the step should be the same as subtracting 
! half the step. 
! If we have tried PROJI with Thomson then the projection is non-linear
! and we need to reset to XSAVE. This should always be reliable!
!
         XCOORDS(1:N)=XSAVE(1:N)
         DO J1=1,N
            XCOORDS(J1)=XCOORDS(J1)+0.5*STP*W(ISPT+POINT*N+J1)
         ENDDO 
         STP=STP/2.0D0
         NDECREASE=NDECREASE+1      
         !FIXIMAGE=.TRUE.
         GOTO 20
      ELSE
!
!  Energy increased - try again with a smaller step size
!
         IF (NDECREASE.GT.10) THEN ! DJW
            NFAIL=NFAIL+1
            WRITE(*,'(A,G20.10)') ' in mylbfgs LBFGS step cannot find a lower energy, NFAIL=',NFAIL

! Resetting to XSAVE should be the same as subtracting the step. 
! If we have tried PROJI with Thomson then the projection is non-linear
! and we need to reset to XSAVE. This should always be reliable!
            XCOORDS(1:N)=XSAVE(1:N)
            GRAD(1:N)=GNEW(1:N) ! GRAD contains the gradient at the lowest energy point
            ITER=0   !  try resetting
!            IF (NFAIL.GT.20) THEN
! bs360: smaller NFAIL 
             IF (NFAIL.GT.5) THEN         
               WRITE(*,'(A)') ' Too many failures - giving up '
               !FIXIMAGE=.FALSE.
!              STOP
               RETURN
            ENDIF
            GOTO 30
         ENDIF
!
! Resetting to XSAVE and adding 0.1 of the step should be the same as subtracting 
! 0.9 of the step. 
! If we have tried PROJI with Thomson then the projection is non-linear
! and we need to reset to XSAVE. This should always be reliable!
!
         XCOORDS(1:N)=XSAVE(1:N)
         DO J1=1,N
            XCOORDS(J1)=XCOORDS(J1)+0.1D0*STP*W(ISPT+POINT*N+J1)
         ENDDO 
         STP=STP/1.0D1
         NDECREASE=NDECREASE+1
         !FIXIMAGE=.TRUE.
         GOTO 20
      ENDIF
!
!     Compute the new step and gradient change
!
30    NPT=POINT*N

      DO J1=1,N
         W(ISPT+NPT+J1)= STP*W(ISPT+NPT+J1) ! save the step taken
         W(IYPT+NPT+J1)= GRAD(J1)-W(J1)     ! save gradient difference: W(1:N) contains the old gradient
      ENDDO
      POINT=POINT+1
      IF (POINT.EQ.M) POINT=0
      !FIXIMAGE=.FALSE.
  
      GOTO 10

      RETURN
      END


      DOUBLE PRECISION FUNCTION DDOT (N, DX, INCX, DY, INCY)
         !js850> my simple implimentation
!     --Input--
!        N  number of elements in input vector(s)
!       DX  double precision vector with N elements
! ignored     INCX  storage spacing between elements of DX
!       DY  double precision vector with N elements
! ignored     INCY  storage spacing between elements of DY
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: N, INCX, INCY
         DOUBLE PRECISION, INTENT(IN) :: DX(N), DY(N)
         INTEGER J1
         DDOT = SUM( DX(:) * DY(:) )
         return
         DDOT = 0.0D0
         IF (N .LE. 0) RETURN
         DO J1=1,N
            DDOT = DDOT + DX(J1) * DY(J1)
         ENDDO
         RETURN
      END FUNCTION DDOT

      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
         IMPLICIT NONE
         INTEGER, INTENT(IN) ::  INCX, INCY, N
         DOUBLE PRECISION, INTENT(IN):: DA
         DOUBLE PRECISION, INTENT(IN) :: DX(N)
         DOUBLE PRECISION, INTENT(INOUT)::  DY(N)
         DY(:) = DY(:) + DA*DX(:)
      END SUBROUTINE DAXPY

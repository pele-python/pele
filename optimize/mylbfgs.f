C   GMIN: A program for finding global minima
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of GMIN.
C
C   GMIN is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   GMIN is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
C
C        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION
C                          JORGE NOCEDAL
C                        *** July 1990 ***
C
C        Line search removed plus small modifications, DJW 2001
C
      SUBROUTINE MYLBFGS(POTENTIAL, N,M,XCOORDS,DIAGCO,EPS,MFLAG,ENERGY,ITMAX,ITDONE,RESET,NP)
      IMPLICIT NONE
CCCCCCCCCCCCCCCCCCCCCCCCCCCcc
      DOUBLE PRECISION DGUESS, MAXBFGS, MAXEFALL, RMS, MAXERISE, SMINKCURRENT
      LOGICAL FIXIMAGE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER N,M,J1,ITMAX,ITDONE,NP,J2,J3,NFAIL,NDECREASE,NGUESS,NDUMMY
      DOUBLE PRECISION XCOORDS(N),GRAD(N),SLENGTH,DDOT,EPLUS,EMINUS,DIFF,DUMMY,WTEMP(N)
      DOUBLE PRECISION TMPANG(N), TMPCOORDS(N)
      DOUBLE PRECISION EPS,ENERGY,ENEW,GNEW(N),OVERLAP,OLDX(N),OLDOLDX(N),VGUESS(3),
     1                 X1, Y1, Z1, X2, Y2, Z2, TRY(N), D1, D2, RBCOORDS(18), DUMMY2, DIST, DIST1
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DIAG, W
      LOGICAL DIAGCO, YESNO, RESET, NOTCALLED, CTEST, MFLAG
      DOUBLE PRECISION GNORM,STP,YS,YY,SQ,YR,BETA,POTEL,QSTART,QFINISH
      DOUBLE PRECISION OLDCART(N), DELTAQ(N),DELTACART(N),LEPSILON,DOT1,DOT2
      DOUBLE PRECISION LCART(N),OLDQ(N),NEWQ(N),XSAVE(N),SMINKCURRENTP
      DOUBLE PRECISION, ALLOCATABLE :: FRAMES(:,:), PE(:), MODGRAD(:)
      LOGICAL NOCOOR, FAILED, COREDONE
      INTEGER ITER,POINT,ISPT,IYPT,BOUND,NPT,CP,INMC,IYCN,ISCN
      INTEGER KD, NNZ, GETUNIT
      LOGICAL EVAP, GUIDECHANGET, GUIDET, EVAPREJECT, SMINKCHANGET, CSMDOGUIDET
      SAVE W, DIAG, ITER, POINT, ISPT, IYPT, NPT


      DGUESS=1.0d0
      MAXBFGS=0.5d0
      MAXEFALL=1.0d0
      RMS=1.0d0

      IF (.NOT.ALLOCATED(DIAG)) ALLOCATE(DIAG(N))       ! SAVE doesn't work otherwise for Sun
      IF (.NOT.ALLOCATED(W)) ALLOCATE(W(N*(2*M+1)+2*M)) ! SAVE doesn't work otherwise for Sun

      IF (SIZE(W,1).NE.N*(2*M+1)+2*M) THEN ! mustn't call mylbfgs with changing number of variables!!!
         WRITE(*, '(A,I10,A,I10,A)') 'ERROR, dimension of W=',SIZE(W,1),' but N*(2*M+1)+2*M=',N*(2*M+1)+2*M,' in mylbfgs'
         call exit(10)
      ENDIF
      COREDONE=.FALSE.
      SMINKCHANGET=.FALSE.
      SMINKCURRENT=0.0D0
      SMINKCURRENTP=0.0D0


      NFAIL=0
      IF (GUIDECHANGET) ITER=0
      IF (RESET) ITER=0
      ITDONE=0
      FIXIMAGE=.FALSE.

      CALL POTENTIAL(XCOORDS,GRAD,ENERGY,.TRUE.,.FALSE.)
      POTEL=ENERGY
      RMS=DOT_PRODUCT(GRAD,GRAD)

C
C  Termination test. 
C
C 10    CALL FLUSH(*)
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
 240        FORMAT(' IMPROPER INPUT PARAMETERS (N OR M ARE NOT POSITIVE)')
            STOP
         ENDIF
         POINT=0
         MFLAG=.FALSE.
         IF (DIAGCO) THEN
            WRITE(*,'(A)') 'using estimate of the inverse diagonal elements'
            DO J1=1,N
               IF (DIAG(J1).LE.0.0D0) THEN
                  WRITE(*,235) J1
 235              FORMAT(' THE',I5,'-TH DIAGONAL ELEMENT OF THE',/,
     1                   ' INVERSE HESSIAN APPROXIMATION IS NOT POSITIVE')
                  STOP
               ENDIF
            ENDDO
         ELSE
            DO J1=1,N
               DIAG(J1)=DGUESS
            ENDDO
         ENDIF

C
C     THE WORK VECTOR W IS DIVIDED AS FOLLOWS:
C     ---------------------------------------
C     THE FIRST N LOCATIONS ARE USED TO STORE THE GRADIENT AND
C         OTHER TEMPORARY INFORMATION.
C     LOCATIONS (N+1)...(N+M) STORE THE SCALARS RHO.
C     LOCATIONS (N+M+1)...(N+2M) STORE THE NUMBERS ALPHA USED
C         IN THE FORMULA THAT COMPUTES H*G.
C     LOCATIONS (N+2M+1)...(N+2M+NM) STORE THE LAST M SEARCH
C         STEPS.
C     LOCATIONS (N+2M+NM+1)...(N+2M+2NM) STORE THE LAST M
C         GRADIENT DIFFERENCES.
C
C     THE SEARCH STEPS AND GRADIENT DIFFERENCES ARE STORED IN A
C     CIRCULAR ORDER CONTROLLED BY THE PARAMETER POINT.
C
         ISPT= N+2*M    ! index for storage of search steps
         IYPT= ISPT+N*M ! index for storage of gradient differences
C
C  NR step for diagonal inverse Hessian
C
         DO J1=1,N
            DUMMY=-GRAD(J1)*DIAG(J1)
            W(ISPT+J1)=DUMMY
            W(J1)=DUMMY
         ENDDO
         GNORM=DSQRT(DDOT(N,GRAD,1,GRAD,1))
C
C  Make the first guess for the step length cautious.
C
         STP=MIN(1.0D0/GNORM,GNORM)
      ELSE 
         BOUND=ITER
         IF (ITER.GT.M) BOUND=M
         YS= DDOT(N,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)
C
C  Update estimate of diagonal inverse Hessian elements
C
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
C           WRITE(*,'(A,2F20.10)') 'YS/YY,STP=',YS/YY,STP
            DO J1=1,N
C              DIAG(J1)= ABS(YS/YY) ! messes up after step reversals!
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
C
C     COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
C     "Updating quasi-Newton matrices with limited storage",
C     Mathematics of Computation, Vol.24, No.151, pp. 773-782.
C     ---------------------------------------------------------
C

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
C
C  Store the new search direction
C
      IF (ITER.GT.0) THEN
         DO J1=1,N
            W(ISPT+POINT*N+J1)= W(J1)
         ENDDO
      ENDIF

      DOT1=SQRT(DDOT(N,GRAD,1,GRAD,1))
C
C  Overflow has occasionally occurred here.
C  We only need the sign of the overlap, so use a temporary array with
C  reduced elements.
C
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


C
C  We now have the proposed step.
C!
! Save XCOORDS here so that we can undo the step reliably including the
! non-linear projection for Thomson for the angular coordinates.
!
      XSAVE(1:N)=XCOORDS(1:N) 
      DO J1=1,N
         XCOORDS(J1)=XCOORDS(J1)+STP*W(ISPT+POINT*N+J1)
      ENDDO 
20    CALL POTENTIAL(XCOORDS,GNEW,ENEW,.TRUE.,.FALSE.)
      RMS=DOT_PRODUCT(GRAD,GRAD)
C
      IF (((ENEW-ENERGY.LE.MAXERISE).OR.EVAP.OR.GUIDECHANGET.OR.SMINKCHANGET).AND.(ENEW-ENERGY.GT.MAXEFALL)) THEN
         ITER=ITER+1
         ITDONE=ITDONE+1
         ENERGY=ENEW
         DO J1=1,N
            GRAD(J1)=GNEW(J1)
         ENDDO
C
C  Step finished so can reset OLDQ to new XINT, OLDCART to new LCART,
C  as well as the Cartesian and internal gradients.
C
      ELSEIF (ENEW-ENERGY.LE.MAXEFALL) THEN
C
C  Energy decreased too much - try again with a smaller step size
C
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
               FIXIMAGE=.FALSE.
!              STOP
!              IF (QUENCHDOS) DEALLOCATE(FRAMES, PE, MODGRAD)
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
C
C  Energy increased - try again with a smaller step size
C
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
               FIXIMAGE=.FALSE.
!              STOP
!              IF (QUENCHDOS) DEALLOCATE(FRAMES, PE, MODGRAD)
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
C
C     Compute the new step and gradient change
C
30    NPT=POINT*N

      DO J1=1,N
         W(ISPT+NPT+J1)= STP*W(ISPT+NPT+J1) ! save the step taken
         W(IYPT+NPT+J1)= GRAD(J1)-W(J1)     ! save gradient difference: W(1:N) contains the old gradient
      ENDDO
      POINT=POINT+1
      IF (POINT.EQ.M) POINT=0
      FIXIMAGE=.FALSE.
  
      SMINKCURRENTP=SMINKCURRENT
      GOTO 10

      RETURN
      END

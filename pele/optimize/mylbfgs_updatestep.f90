!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!start a new subroutine here
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE MYLBFGS_UPDATESTEP(ITER, N, M, GRAD, W, DIAG, POINT, &
         STPOUT)
      IMPLICIT NONE
      INTEGER, INTENT(IN)             :: N
      INTEGER, INTENT(IN)             :: M
      INTEGER, INTENT(IN)             :: ITER
      INTEGER, INTENT(INOUT)          :: POINT
      !DOUBLE PRECISION, INTENT(IN)    :: X(N)
      DOUBLE PRECISION, INTENT(IN)    :: GRAD(N)
      DOUBLE PRECISION, INTENT(INOUT) :: W(N*(2*M+1)+2*M)
      DOUBLE PRECISION, INTENT(INOUT) :: DIAG(N)
      DOUBLE PRECISION                :: STP  !not passed
      DOUBLE PRECISION, INTENT(OUT)   :: STPOUT(N)
      !INTEGER, INTENT(OUT)            :: error ! if not zero something went wrong

      INTEGER J1
      INTEGER NPT  !js850>  I think this is (POINT-1)*N
      INTEGER ISPT, IYPT, CP, INMC, ISCN, IYCN, BOUND
      DOUBLE PRECISION DGUESS, DUMMY, WTEMP(N), DOT1, DOT2, GNORM
      DOUBLE PRECISION OVERLAP, SQ, BETA
      DOUBLE PRECISION YR, YS, YY

      LOGICAL DIAGCO !this could be passed. It would be used to specify what the diagonal components are
      LOGICAL DIAGCO_STEP0 !use DIAG input only on initial step
      DOUBLE PRECISION DDOT ! A FUNCTION
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
      !error = 0
      ! error = 1 # Improper input parameters (N or M are not positive)
      ! error = 2 # A diagonal element of the inverse hessian approximation is not positive

      ISPT= N+2*M    ! index for storage of search steps
      IYPT= ISPT+N*M ! index for storage of gradient differences

      DGUESS=1.0d0
      DIAGCO = .FALSE.
      DIAGCO_STEP0 = .TRUE.

 240        FORMAT(' IMPROPER INPUT PARAMETERS (N OR M ARE NOT POSITIVE)')
 235              FORMAT(' THE',I5,'-TH DIAGONAL ELEMENT OF THE',/, &
                        ' INVERSE HESSIAN APPROXIMATION IS NOT POSITIVE')

      IF (ITER.EQ.0) THEN

         !a bug check
         IF (N.LE.0.OR.M.LE.0) THEN
            WRITE(*,240)
            STOP
         ENDIF

         POINT=0

         !set the inverse hessian
         IF (DIAGCO .OR. DIAGCO_STEP0) THEN
            !WRITE(*,*) 'using estimate of the inverse diagonal elements', diag(1)
            DO J1=1,N
               IF (DIAG(J1).LE.0.0D0) THEN
                  WRITE(*,*) "The diagonal inverse hessian approximation is negative, using", dguess, "instead"
                  DIAG(:) = DGUESS
                  exit ! exit loop
                  !STOP
               ENDIF
            ENDDO
         ELSE
            DO J1=1,N
               DIAG(J1)=DGUESS
            ENDDO
         ENDIF

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

      ELSE !iter > 0
         !IF (POINT == 0) THEN
            !NPT = (POINT-1)*N !js850> I'm not sure about this
         !ELSE
            !ENDIF
         NPT = N*MOD(POINT + M - 1, M)

         BOUND=ITER
         IF (ITER.GT.M) BOUND=M
         !CALL FLUSH()
         YS= DDOT(N,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)
!
!  Update estimate of diagonal inverse Hessian elements
!
         IF (.NOT.DIAGCO) THEN
            YY= DDOT(N,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)
            !YY = SUM( W( IYPT+NPT+1 : IYPT+NPT+N )**2 )
            !write(*,*) "YS YY fo",  YS, YY, ISPT+NPT+1
            IF (YY.EQ.0.0D0) THEN
               WRITE(*,'(A)') 'WARNING, resetting YY to one in mylbfgs'
               YY=1.0D0
            ENDIF
            IF (YS.EQ.0.0D0) THEN
               WRITE(*,'(A)') 'WARNING, resetting YS to one in mylbfgs'
               YS=1.0D0
            ENDIF
            !write(*,*) W(IYPT+NPT+1:IYPT+NPT+N)
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

      DO J1=1,N
         STPOUT(J1) = STP*W(ISPT+POINT*N+J1)
      ENDDO


   END SUBROUTINE MYLBFGS_UPDATESTEP

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

subroutine lbfgs_get_step( G, s, y, rho, k, H0, N, M, step)
implicit none
integer, intent(IN) :: N, M, k
double precision, intent(IN) :: G(N)
double precision, intent(IN) :: s(N,M)
double precision, intent(IN) :: y(N,M)
double precision, intent(IN) :: rho(M)
double precision, intent(IN) :: H0
double precision, intent(OUT) :: step(N)

double precision :: a(N), sq, yz, beta, gnorm
integer i, j1, jstart, jstop

step(:) = G(:)

if (k.eq.0) then
   gnorm = sqrt(sum(g(:) * g(:)))
   gnorm = min(gnorm, 1.d0 / gnorm)
   step(:) = -gnorm * H0 * step(:)
   return
endif

jstart = max(1, k - M + 1)
jstop = k

! loop through the history, most recent first
do j1 = jstop, jstart, -1
   i = mod(j1, M)
   if (i .eq. 0) i = M
   sq = sum(s(:,i) * step(:)) 
   a(i) = rho(i) * sq
   step(:) = step(:) - a(i) * y(:,i)
enddo

! include our estimate for diagonal component of the inverse hessian
step(:) = step(:) * H0

! loop through the history, most recent last
   !write(*,*) " "
do j1 = jstart, jstop
   !write(*,*) "j1"
   i = mod(j1, M)
   if (i .eq. 0) i = M
   yz = sum(y(:,i) * step(:))
   beta = rho(i) * yz
   step(:) = step(:) + s(:,i) * (a(i) - beta)
enddo
   !write(*,*) " "

! step should point downhill
step(:) = -step(:)

! make first guess for the step length cautious
!if (k .eq. 0) then
   !gnorm = sqrt(sum(g(:) * g(:)))
   !gnorm = min(gnorm, 1.d0 / gnorm)
   !step(:) = step(:) * gnorm
!endif

end subroutine lbfgs_get_step

subroutine lbfgs_get_step_wrapper( G, s, y, rho, k, H0, N, M, step)
! wrap lbfgs_get_step so that only 1d arrays are passed
implicit none
integer, intent(IN) :: N, M, k
double precision, intent(IN) :: G(N)
double precision, intent(IN) :: s(N*M)
double precision, intent(IN) :: y(N*M)
double precision, intent(IN) :: rho(M)
double precision, intent(IN) :: H0
double precision, intent(OUT) :: step(N)
call lbfgs_get_step( G, s, y, rho, k, H0, N, M, step)
end subroutine lbfgs_get_step_wrapper

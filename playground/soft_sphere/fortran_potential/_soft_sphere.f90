      SUBROUTINE soft_sphere(N, X, gradient, hess, energy, BOXLX, BOXLY, BOXLZ, &
            radii, power, GTEST, STEST)
      IMPLICIT NONE
      INTEGER, intent(IN) :: N
      DOUBLE PRECISION, intent(IN) :: X(3*N), BOXLX, BOXLY, BOXLZ, radii(N), power
      DOUBLE PRECISION, intent(OUT) :: gradient(3*N), energy, hess(3*N,3*N)
      LOGICAL, INTENT(IN) :: GTEST, STEST
      INTEGER J1, J2, J3, J4, k, l 
      DOUBLE PRECISION dr(3), boxvec(3), eps, gij, hij, r, r0, drpow
      DOUBLE PRECISION hii_diag, hij_diag
      DOUBLE PRECISION hii_off, hij_off
      boxvec(1) = boxlx
      boxvec(2) = boxly
      boxvec(3) = boxlz

      eps = 1.d0
      energy=0.0D0
      gradient(:) = 0.d0
      if (stest) hess(:,:) = 0.d0
      DO J1=1,N
         J3=3*(J1-1)
         DO J2=J1+1,N
            J4=3*(J2-1)

            do k=1,3
               dr(k)=X(J3+k)-X(J4+k)
               dr(k)=dr(k) - BOXLX * NINT(dr(k) / BOXLX) 
            enddo


            r0 = radii(j1) + radii(j2)
            r = sqrt(sum(dr(:)**2))
            drpow = (1.d0 - r/r0)**power

            IF (r.lt.r0) THEN
               ! do energy
               energy=energy + drpow * eps / power
               
               ! do gradient
               gij = -eps * drpow  / (r * (r-r0))
               do k = 1,3
                  gradient(J3+k) = gradient(J3+k) - gij * dr(k)
                  gradient(J4+k) = gradient(J4+k) + gij * dr(k)
               enddo

               ! do hessian
               if (stest) then
                  hij = eps * (power-1) * drpow  / ( ((r-r0)) * ((r-r0)) )

                  do k=1,3
                     !diagonal block - diagonal terms
                     Hii_diag = (hij+gij)*dr(k)*dr(k)/(r*r) - gij
                     hess(j3+k, j3+k) = hess(j3+k, j3+k) + Hii_diag
                     hess(j4+k, j4+k) = hess(j4+k, j4+k) + Hii_diag
                     !off diagonal block - diagonal terms
                     Hij_diag = -Hii_diag
                     hess(j3+k, j4+k) = Hij_diag
                     hess(j4+k, j3+k) = Hij_diag
                     do l = k+1,3
                        !diagonal block - off diagonal terms
                        Hii_off = (hij+gij)*dr(k)*dr(l)/(r*r)
                        hess(j3+k, j3+l) = hess(j3+k, j3+l) + Hii_off
                        hess(j3+l, j3+k) = hess(j3+l, j3+k) + Hii_off
                        hess(j4+k, j4+l) = hess(j4+k, j4+l) + Hii_off
                        hess(j4+l, j4+k) = hess(j4+l, j4+k) + Hii_off
                        !off diagonal block - off diagonal terms
                        Hij_off = -Hii_off
                        hess(j3+k, j4+l) = Hij_off
                        hess(j3+l, j4+k) = Hij_off
                        hess(j4+k, j3+l) = Hij_off
                        hess(j4+l, j3+k) = Hij_off
                     enddo
                  enddo

               endif
            ENDIF
         ENDDO
      ENDDO
      END SUBROUTINE soft_sphere

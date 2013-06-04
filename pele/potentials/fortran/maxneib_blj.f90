function fermi(x, mu, T) result(ferm)
   implicit none
double precision, intent(in) :: x, mu, T
double precision ferm
ferm = 1.d0 / (exp(-(x - mu) / T) + 1.d0)
end function

function dfermi_dx(x, mu, T) result(ferm)
   implicit none
double precision, intent(in) :: x, mu, T
double precision ferm, myexp
myexp = exp(-(x - mu) / T)

ferm = myexp * (1.d0/T) / (myexp + 1.d0)**2
end function


subroutine maxneib_ljenergy( coords, natoms, e, ntypea, &
   epsa, siga, &
   epsb, sigb, &
   epsab, sigab, &
   periodic, boxl, &
   rneib, rneib_crossover, max_neibs, neib_crossover, epsneibs, &
   only_AB_neibs)
implicit none
integer, intent(in) :: natoms, ntypea
double precision, intent(in) :: coords(3*natoms), boxl
double precision, intent(in) :: siga, epsa
double precision, intent(in) :: sigb, epsb
double precision, intent(in) :: sigab, epsab
double precision, intent(out) :: e
double precision dr(3), r2, ir2, ir6, ir12, iboxl
logical, intent(in) :: periodic, only_AB_neibs
integer j1, j2
double precision sig, eps
double precision, intent(in) :: rneib, rneib_crossover, max_neibs
double precision, intent(in) :: neib_crossover, epsneibs
double precision nneibs(natoms)
double precision eneibs, eneibs1, areneibs
double precision fermi

nneibs(:) = 0.d0

!rneib = 1.4
!rneib_crossover = .08
!max_neibs = 5.
!neib_crossover = .4
!epsneibs = 1.


if (periodic) iboxl = 1.d0 / boxl

e = 0.d0
do j1 = 1,natoms
   do j2 = 1,j1-1
      if (j1 .le. ntypea .and. j2 .le. ntypea) then
         sig = siga
         eps = epsa
      elseif (j1 .gt. ntypea .and. j2 .gt. ntypea) then
         sig = sigb
         eps = epsb
      else
         sig = sigab
         eps = epsab
      endif
      dr(:) = coords(3*(j1-1)+1 : 3*(j1-1) + 3) - coords(3*(j2-1)+1 : 3*(j2-1) + 3)
      if (periodic)  dr(:) = dr(:) - nint( dr(:) * iboxl ) * boxl
      r2 = sum( dr(:)**2 )
      r2 = r2 / sig**2
      ir2 = 1.d0 / r2
      ir6 = ir2**3
      ir12 = ir6**2
      e = e - 4.d0 * eps * (ir6 - ir12)
      if (only_AB_neibs) then
         if (j1 .le. ntypea .and. j2 .le. ntypea) then
            areneibs = 0.d0
         elseif (j1 .gt. ntypea .and. j2 .gt. ntypea) then
            areneibs = 0.d0
         else
            areneibs = 1.d0 - fermi(sqrt(r2), rneib, rneib_crossover)
         endif
      else
         areneibs = 1.d0 - fermi(sqrt(r2), rneib, rneib_crossover)
      endif
      nneibs(j1) = nneibs(j1)  + areneibs
      nneibs(j2) = nneibs(j2)  + areneibs
   enddo
enddo

!calculate energy from neighbors
eneibs = 0.d0
do j1=1,natoms
   eneibs1 = fermi(nneibs(j1), max_neibs, neib_crossover)
   eneibs = eneibs + eneibs1
enddo
eneibs = eneibs * epsneibs

e = e + eneibs

end subroutine maxneib_ljenergy

subroutine maxneib_ljenergy_gradient( coords, natoms, e, grad, ntypea, &
   epsa, siga, &
   epsb, sigb, &
   epsab, sigab, &
   periodic, boxl, &
   rneib, rneib_crossover, max_neibs, neib_crossover, epsneibs, &
   only_AB_neibs)
implicit none
integer, intent(in) :: natoms, ntypea
double precision, intent(in) :: coords(3*natoms), boxl
double precision, intent(in) :: siga, epsa
double precision, intent(in) :: sigb, epsb
double precision, intent(in) :: sigab, epsab
double precision, intent(out) :: e, grad(3*natoms)
double precision dr(3), r2, ir2, ir6, ir12, iboxl
logical, intent(in) :: periodic, only_AB_neibs
integer j1, j2, i1, i2
double precision sig, eps
double precision, intent(in) :: rneib, rneib_crossover, max_neibs
double precision, intent(in) :: neib_crossover, epsneibs
double precision nneibs(natoms)
double precision eneibs, eneibs1, areneibs
double precision fermi, dfermi_dx
double precision r, dfr, dfn1, dfn2, g

nneibs(:) = 0.d0

!rneib = 1.4
!rneib_crossover = .08
!max_neibs = 5.
!neib_crossover = .4
!epsneibs = 1.


if (periodic) iboxl = 1.d0 / boxl


!get energy and gradient of lj.
! also calculate the number of neighbors of each atom
e = 0.d0
grad(:) = 0.d0
do j1 = 1,natoms
   i1 = 3*(j1-1)
   do j2 = 1,j1-1
      i2 = 3*(j2-1)
      if (j1 .le. ntypea .and. j2 .le. ntypea) then
         sig = siga
         eps = epsa
      elseif (j1 .gt. ntypea .and. j2 .gt. ntypea) then
         sig = sigb
         eps = epsb
      else
         sig = sigab
         eps = epsab
      endif
      dr(:) = coords(i1+1 : i1 + 3) - coords(i2+1 : i2 + 3)
      if (periodic)  dr(:) = dr(:) - nint( dr(:) * iboxl ) * boxl
      r2 = sum( dr(:)**2 )
      r2 = r2 / sig**2
      ir2 = 1.d0 / r2
      ir6 = ir2**3
      ir12 = ir6**2
      e = e - 4.d0 * eps * (ir6 - ir12)
      g = 4.d0 * eps * (12.d0 * ir12 -  6.d0 * ir6) * ir2 / sig**2;

      grad(i1+1 : i1+3) = grad(i1+1 : i1+3) - g * dr(:)
      grad(i2+1 : i2+3) = grad(i2+1 : i2+3) + g * dr(:)

      if (only_AB_neibs) then
         if (j1 .le. ntypea .and. j2 .le. ntypea) then
            areneibs = 0.d0
         elseif (j1 .gt. ntypea .and. j2 .gt. ntypea) then
            areneibs = 0.d0
         else
            areneibs = 1.d0 - fermi(sqrt(r2), rneib, rneib_crossover)
         endif
      else
         areneibs = 1.d0 - fermi(sqrt(r2), rneib, rneib_crossover)
      endif
      !write(*,*) "areneibs", areneibs
      nneibs(j1) = nneibs(j1)  + areneibs
      nneibs(j2) = nneibs(j2)  + areneibs
   enddo
enddo

!calculate gradient from neighbors
do j1 = 1,natoms
   i1 = 3*(j1-1)
   dfn1 = dfermi_dx(nneibs(j1), max_neibs, neib_crossover)
   !write(*,*) nneibs(j1), dfn1
   do j2 = 1,j1-1
      i2 = 3*(j2-1)
      if (j1 .le. ntypea .and. j2 .le. ntypea) then
         sig = siga
      elseif (j1 .gt. ntypea .and. j2 .gt. ntypea) then
         sig = sigb
      else
         sig = sigab
      endif
      if (only_AB_neibs) then
         if (j1 .le. ntypea .and. j2 .le. ntypea) then
            exit
         elseif (j1 .gt. ntypea .and. j2 .gt. ntypea) then
            exit
         endif
      endif
      dr(:) = coords(i1+1 : i1 + 3) - coords(i2+1 : i2 + 3)
      if (periodic)  dr(:) = dr(:) - nint( dr(:) * iboxl ) * boxl
      r2 = sum( dr(:)**2 )
      r2 = r2 / sig**2
      r = sqrt(r2)

      dfr = dfermi_dx(r, rneib, rneib_crossover)
      dfn2 = dfermi_dx(nneibs(j2), max_neibs, neib_crossover)
      g = epsneibs * dfr * (dfn1 + dfn2)
      !write(*,*) dfr, dfn1, dfn2
      !write(*,*) nneibs(j2), dfn2
      grad(i1+1 : i1+3) = grad(i1+1 : i1+3) - g * dr(:) / r / sig**2
      grad(i2+1 : i2+3) = grad(i2+1 : i2+3) + g * dr(:) / r / sig**2
   enddo
enddo


!calculate energy from neighbors
eneibs = 0.d0
do j1=1,natoms
   eneibs1 = fermi(nneibs(j1), max_neibs, neib_crossover)
   eneibs = eneibs + eneibs1
enddo
eneibs = eneibs * epsneibs

e = e + eneibs

end subroutine maxneib_ljenergy_gradient

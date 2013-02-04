function fermi(x, mu, T) result(ferm)
double precision, intent(in) :: x, mu, T
double precision ferm
ferm = 1.d0 / (exp(-(x - mu) / T) + 1.d0)
end function

subroutine maxneib_ljenergy( coords, natoms, e, eps, sig, periodic, boxl, &
   rneib, rneib_crossover, max_neibs, neib_crossover, epsneibs)
implicit none
integer, intent(in) :: natoms
double precision, intent(in) :: coords(3*natoms), sig, eps, boxl
double precision, intent(out) :: e
double precision dr(3), sig6, sig12, r2, ir2, ir6, ir12, iboxl
logical, intent(in) :: periodic
integer j1, j2
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

sig6 = sig**6
sig12 = sig6*sig6

e = 0.d0
do j1 = 1,natoms
   do j2 = 1,j1-1
      dr(:) = coords(3*(j1-1)+1 : 3*(j1-1) + 3) - coords(3*(j2-1)+1 : 3*(j2-1) + 3)
      if (periodic)  dr(:) = dr(:) - nint( dr(:) * iboxl ) * boxl
      r2 = sum( dr(:)**2 )
      ir2 = 1.d0/r2
      ir6 = ir2**3
      ir12 = ir6**2
      e = e - 4.d0 * eps * (sig6*ir6 - sig12*ir12)
      areneibs = 1.d0 - fermi(sqrt(r2), rneib, rneib_crossover)
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

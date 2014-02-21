subroutine mypot_energy_gradient( coords, natoms, e, grad, eps, sig )
implicit none
integer, intent(in) :: natoms
double precision, intent(in) :: coords(3*natoms), sig, eps
double precision, intent(out) :: e, grad(3*natoms)
double precision dr(3), sig12, sig24, r2, ir2, ir12, ir24, g
integer j1, j2, i1, i2


sig12 = sig**12
sig24 = sig12*sig12

e = 0.d0
grad(:) = 0.d0
do j1 = 1,natoms
   i1 = 3*(j1-1)
   do j2 = 1,j1-1
      i2 = 3*(j2-1)
      dr(:) = coords(i1+1 : i1 + 3) - coords(i2+1 : i2 + 3)
      r2 = sum( dr(:)**2 )
      ir2 = 1.d0/r2
      ir12 = ir2**6
      ir24 = ir12**2
      e = e - 4.d0 * eps * (sig12*ir12 - sig24*ir24)

      g = 4.d0 * eps * (24.d0 * sig24 * ir24 -  12.d0 * sig12 * ir12) * ir2;
      grad(i1+1 : i1+3) = grad(i1+1 : i1+3) - g * dr(:)
      grad(i2+1 : i2+3) = grad(i2+1 : i2+3) + g * dr(:)
   enddo
enddo
end subroutine mypot_energy_gradient

subroutine ljenergy( coords, natoms, e, eps, sig )
implicit none
integer, intent(in) :: natoms
double precision, intent(in) :: coords(3*natoms), sig, eps
double precision, intent(out) :: e
double precision dr(3), sig6, sig12, r2, ir2, ir6, ir12
integer j1, j2

sig6 = sig**6
sig12 = sig6*sig6

e = 0.d0
do j1 = 1,natoms
   do j2 = 1,j1-1
      dr(:) = coords(3*(j1-1)+1 : 3*(j1-1) + 3) - coords(3*(j2-1)+1 : 3*(j2-1) + 3)
      r2 = sum( dr(:)**2 )
      ir2 = 1.d0/r2
      ir6 = ir2**3
      ir12 = ir6**2
      e = e - 4.d0 * eps * (sig6*ir6 - sig12*ir12)
   enddo
enddo

end subroutine ljenergy

subroutine ljenergy_gradient( coords, natoms, e, grad, eps, sig )
implicit none
integer, intent(in) :: natoms
double precision, intent(in) :: coords(3*natoms), sig, eps
double precision, intent(out) :: e, grad(3*natoms)
double precision dr(3), sig6, sig12, r2, ir2, ir6, ir12, g
integer j1, j2, i1, i2

sig6 = sig**6
sig12 = sig6*sig6

e = 0.d0
grad(:) = 0.d0
do j1 = 1,natoms
   i1 = 3*(j1-1)
   do j2 = 1,j1-1
      i2 = 3*(j2-1)
      dr(:) = coords(i1+1 : i1 + 3) - coords(i2+1 : i2 + 3)
      r2 = sum( dr(:)**2 )
      ir2 = 1.d0/r2
      ir6 = ir2**3
      ir12 = ir6**2
      e = e - 4.d0 * eps * (sig6*ir6 - sig12*ir12)

      g = 4.d0 * eps * (12.d0 * sig12 * ir12 -  6.d0 * sig6 * ir6) * ir2;
      grad(i1+1 : i1+3) = grad(i1+1 : i1+3) - g * dr(:)
      grad(i2+1 : i2+3) = grad(i2+1 : i2+3) + g * dr(:)
   enddo
enddo
end subroutine ljenergy_gradient

!!!
! lj with a cutoff that is continuous and smooth
!!!

!
!note: it happens quite often that integers are defined differently in python
!and fortran.  (e.g. 32 bit in fortran and 64 bit in python).  This can be a
!problem because if a large array of integers needs to be passed, it must be
!coppied first.  The solution I've found is to define integers to be kind=8,
!then make sure they're defined in python with np.array( ... , np.int64)
!

subroutine ljenergy( coords, natoms, e, eps, sig, periodic, boxl, rcut )
implicit none
integer, intent(in) :: natoms
double precision, intent(in) :: coords(3*natoms), sig, eps, boxl, rcut
double precision, intent(out) :: e
double precision dr(3), sig6, sig12, r2, ir2, ir6, ir12, iboxl
logical, intent(in) :: periodic
integer j1, j2
double precision rcut2, rcut6, A1, B1

if (periodic) iboxl = 1.d0 / boxl

sig6 = sig**6
sig12 = sig6*sig6
rcut2 = rcut**2
rcut6 = rcut**6
A1 = 4.0D0*(sig6/rcut6) - 7.0D0*(sig12/rcut6**2)
B1 = (-3.0D0*(sig6/rcut6) + 6.0D0*(sig12/rcut6**2)) * (1.d0/rcut)**2

e = 0.d0
do j1 = 1,natoms
   do j2 = 1,j1-1
      dr(:) = coords(3*(j1-1)+1 : 3*(j1-1) + 3) - coords(3*(j2-1)+1 : 3*(j2-1) + 3)
      if (periodic)  dr(:) = dr(:) - nint( dr(:) * iboxl ) * boxl
      r2 = sum( dr(:)**2 )
      if (r2 .le. rcut2) then
         ir2 = 1.d0/r2
         ir6 = ir2**3
         ir12 = ir6**2
         e = e - 4.d0 * eps * (sig6*ir6 - sig12*ir12 - A1 - B1*r2)
      endif
   enddo
enddo

end subroutine ljenergy

subroutine ljenergy_gradient( coords, natoms, e, grad, eps, sig, periodic, boxl, rcut )
implicit none
integer, intent(in) :: natoms
double precision, intent(in) :: coords(3*natoms), sig, eps, boxl, rcut
double precision, intent(out) :: e, grad(3*natoms)
logical, intent(in) :: periodic
double precision dr(3), sig6, sig12, r2, ir2, ir6, ir12, g, iboxl
integer j1, j2, i1, i2
double precision rcut2, rcut6, A1, B1

if (periodic) iboxl = 1.d0 / boxl

sig6 = sig**6
sig12 = sig6*sig6
rcut2 = rcut**2
rcut6 = rcut**6
A1 = 4.0D0*(sig6/rcut6) - 7.0D0*(sig12/rcut6**2)
B1 = (-3.0D0*(sig6/rcut6) + 6.0D0*(sig12/rcut6**2)) * (1.d0/rcut)**2

e = 0.d0
grad(:) = 0.d0
do j1 = 1,natoms
   i1 = 3*(j1-1)
   do j2 = 1,j1-1
      i2 = 3*(j2-1)
      dr(:) = coords(i1+1 : i1 + 3) - coords(i2+1 : i2 + 3)
      if (periodic)  dr(:) = dr(:) - nint( dr(:) * iboxl ) * boxl
      r2 = sum( dr(:)**2 )
      if (r2 .le. rcut2) then
         ir2 = 1.d0/r2
         ir6 = ir2**3
         ir12 = ir6**2
         e = e - 4.d0 * eps * (sig6*ir6 - sig12*ir12 - A1 - B1*r2)

         g = 4.d0 * eps * ((12.d0 * sig12 * ir12 -  6.d0 * sig6 * ir6) * ir2 - 2.d0*B1)
         grad(i1+1 : i1+3) = grad(i1+1 : i1+3) - g * dr(:)
         grad(i2+1 : i2+3) = grad(i2+1 : i2+3) + g * dr(:)
      endif
   enddo
enddo
end subroutine ljenergy_gradient

subroutine energy_ilist( coords, natoms, e, eps, sig, ilist, nlist, periodic, boxl, rcut )
implicit none
integer(kind=8), intent(in) :: natoms, nlist, ilist(nlist* 2)
double precision, intent(in) :: coords(3*natoms), sig, eps, boxl, rcut
double precision, intent(out) :: e
logical, intent(in) :: periodic
double precision dr(3), sig6, sig12, r2, ir2, ir6, ir12, iboxl
integer(kind=8) j1, j2, i1, i2, n
double precision rcut2, rcut6, A1, B1

if (periodic) iboxl = 1.d0 / boxl

sig6 = sig**6
sig12 = sig6*sig6
rcut2 = rcut**2
rcut6 = rcut**6
A1 = 4.0D0*(sig6/rcut6) - 7.0D0*(sig12/rcut6**2)
B1 = (-3.0D0*(sig6/rcut6) + 6.0D0*(sig12/rcut6**2)) * (1.d0/rcut2)

e = 0.d0

do n = 1,nlist
   j1 = ilist(2*(n-1)+1) + 1 !convert to fortran indexing
   j2 = ilist(2*(n-1)+2) + 1 !convert to fortran indexing

   i1 = 3*(j1-1)
   i2 = 3*(j2-1)

   dr(:) = coords(i1+1 : i1 + 3) - coords(i2+1 : i2 + 3)
   if (periodic)  dr(:) = dr(:) - nint( dr(:) * iboxl ) * boxl
   r2 = sum( dr(:)**2 )
   if (r2 .le. rcut2) then
      ir2 = 1.d0/r2
      ir6 = ir2**3
      ir12 = ir6**2
      e = e - 4.d0 * eps * (sig6*ir6 - sig12*ir12 - A1 - B1*r2)
   endif

enddo
end subroutine energy_ilist

subroutine energy_gradient_ilist( coords, natoms, e, grad, eps, sig, ilist, nlist, periodic, boxl, rcut )
implicit none
integer(kind=8), intent(in) :: natoms, nlist, ilist(nlist* 2)
double precision, intent(in) :: coords(3*natoms), sig, eps, boxl, rcut
double precision, intent(out) :: e, grad(3*natoms)
logical, intent(in) :: periodic
double precision dr(3), sig6, sig12, r2, ir2, ir6, ir12, g, iboxl
integer(kind=8) j1, j2, i1, i2, n
double precision rcut2, rcut6, A1, B1

if (periodic) iboxl = 1.d0 / boxl

sig6 = sig**6
sig12 = sig6*sig6
rcut2 = rcut**2
rcut6 = rcut2**3
A1 = 4.0D0*(sig6/rcut6) - 7.0D0*(sig12/rcut6**2)
B1 = (-3.0D0*(sig6/rcut6) + 6.0D0*(sig12/rcut6**2)) * (1.d0/rcut)**2

e = 0.d0
grad(:) = 0.d0

do n = 1,nlist
   j1 = ilist(2*(n-1)+1) + 1 !convert to fortran indexing
   j2 = ilist(2*(n-1)+2) + 1 !convert to fortran indexing

   i1 = 3*(j1-1)
   i2 = 3*(j2-1)

   dr(:) = coords(i1+1 : i1 + 3) - coords(i2+1 : i2 + 3)
   if (periodic)  dr(:) = dr(:) - nint( dr(:) * iboxl ) * boxl
   r2 = sum( dr(:)**2 )
   if (r2 .le. rcut2) then
      ir2 = 1.d0/r2
      ir6 = ir2**3
      ir12 = ir6**2
      e = e - 4.d0 * eps * (sig6*ir6 - sig12*ir12 - A1 - B1*r2)

      g = 4.d0 * eps * ((12.d0 * sig12 * ir12 -  6.d0 * sig6 * ir6) * ir2 - 2.d0*B1)
      grad(i1+1 : i1+3) = grad(i1+1 : i1+3) - g * dr(:)
      grad(i2+1 : i2+3) = grad(i2+1 : i2+3) + g * dr(:)
   endif
enddo
end subroutine energy_gradient_ilist

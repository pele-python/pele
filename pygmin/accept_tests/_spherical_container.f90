subroutine check_sphereical_container(coords, N, radius, res)
   implicit none
   integer, intent(in) :: N
   double precision, intent(in) :: coords(3*N), radius
   logical, intent(out) :: res
   integer i, j
   double precision radius2, r2
   radius2 = radius**2
   res = .true.
   do i=1,N
      j = i*3 - 2
      r2 = sum(coords(j:j+2)**2)
      if (r2 .ge. radius2) then
         res = .false.
         return
      endif
   enddo
end subroutine check_sphereical_container

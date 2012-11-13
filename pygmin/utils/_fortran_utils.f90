
!
!note: it happens quite often that integers are defined differently in python
!and fortran.  (e.g. 32 bit in fortran and 64 bit in python).  This can be a
!problem because if a large array of integers needs to be passed, it must be
!coppied first.  The solution I've found is to define integers to be kind=8,
!then make sure they're defined in python with np.array( ... , np.int64)
!

subroutine build_neighbor_list1(coords, natoms, atomlist, natomlist, list, nlistmax, nlist, rlist2)
implicit none
integer(kind=8), intent(in) :: natoms, natomlist, atomlist(natomlist), nlistmax
double precision, intent(in) :: coords(3*natoms), rlist2
integer(kind=8), intent(out) :: list(nlistmax)
integer(kind=8), intent(out) :: nlist
integer(kind=8) k1, k2, j1, j2
double precision r2
!note that atomlist is indexed as in python, e.g. 0,1,2,
!write(*,*) "in build_neighbor_list1"
nlist = 0
do k1=1,natomlist
    j1 = atomlist(k1) 
    do k2 = 1,k1-1
        j2 = atomlist(k2)
        r2 = sum( (coords(3*(j1)+1 : 3*(j1)+3) - coords(3*(j2)+1 : 3*(j2)+3))**2 )
        if (r2 .le. rlist2) then
            list(nlist*2+1) = j1
            list(nlist*2+2) = j2
            nlist = nlist + 1
        endif
    enddo
enddo
end subroutine build_neighbor_list1

subroutine build_neighbor_list2(coords, natoms, Alist, nAlist, Blist, nBlist, list, nlistmax, nlist, rlist2)
implicit none
integer(kind=8), intent(in) :: natoms, nAlist, Alist(nAlist), nBlist, Blist(nBlist), nlistmax
double precision, intent(in) :: coords(3*natoms), rlist2
integer(kind=8), intent(out) :: list(nlistmax)
integer(kind=8), intent(out) :: nlist
integer(kind=8) k1, k2, j1, j2
double precision r2
!note that Alist is indexed as in python, e.g. 0,1,2,
!write(*,*) "in build_neighbor_list2"
nlist = 0
do k1=1,nAlist
   j1 = Alist(k1) 
   do k2=1,nBlist
        j2 = Blist(k2) 
        r2 = sum( (coords(3*(j1)+1 : 3*(j1)+3) - coords(3*(j2)+1 : 3*(j2)+3))**2 )
        if (r2 .le. rlist2) then
            list(nlist*2+1) = j1
            list(nlist*2+2) = j2
            nlist = nlist + 1
        endif
    enddo
enddo
end subroutine build_neighbor_list2

subroutine build_neighbor_list1_periodic(coords, natoms, atomlist, natomlist, list, nlistmax, nlist, rlist2, boxl)
implicit none
integer(kind=8), intent(in) :: natoms, natomlist, atomlist(natomlist), nlistmax
double precision, intent(in) :: coords(3*natoms), rlist2, boxl
integer(kind=8), intent(out) :: list(nlistmax)
integer(kind=8), intent(out) :: nlist
integer(kind=8) k1, k2, j1, j2
double precision r2, dr(3), iboxl
!write(*,*) "in build_neighbor_list1_periodic"
iboxl = 1.d0/boxl
!note that atomlist is indexed as in python, e.g. 0,1,2,
nlist = 0
do k1=1,natomlist
    j1 = atomlist(k1) 
    do k2 = 1,k1-1
        j2 = atomlist(k2)
        dr = (coords(3*(j1)+1 : 3*(j1)+3) - coords(3*(j2)+1 : 3*(j2)+3))
        dr = dr - boxl * nint( dr * iboxl )
        r2 = sum( dr**2 )
        if (r2 .le. rlist2) then
            list(nlist*2+1) = j1
            list(nlist*2+2) = j2
            nlist = nlist + 1
        endif
    enddo
enddo
end subroutine build_neighbor_list1_periodic

subroutine build_neighbor_list2_periodic(coords, natoms, Alist, nAlist, Blist, nBlist, list, nlistmax, nlist, rlist2, boxl)
implicit none
integer(kind=8), intent(in) :: natoms, nAlist, Alist(nAlist), nBlist, Blist(nBlist), nlistmax
double precision, intent(in) :: coords(3*natoms), rlist2, boxl
integer(kind=8), intent(out) :: list(nlistmax)
integer(kind=8), intent(out) :: nlist
integer(kind=8) k1, k2, j1, j2
double precision r2, dr(3), iboxl
!write(*,*) "in build_neighbor_list2_periodic"
iboxl = 1.d0/boxl
!note that Alist is indexed as in python, e.g. 0,1,2,
nlist = 0
do k1=1,nAlist
   j1 = Alist(k1) 
   do k2=1,nBlist
        j2 = Blist(k2) 
        dr = (coords(3*(j1)+1 : 3*(j1)+3) - coords(3*(j2)+1 : 3*(j2)+3))
        dr = dr - boxl * nint( dr * iboxl )
        r2 = sum( dr**2 )
        if (r2 .le. rlist2) then
            list(nlist*2+1) = j1
            list(nlist*2+2) = j2
            nlist = nlist + 1
        endif
    enddo
enddo
end subroutine build_neighbor_list2_periodic

subroutine check_neighbor_lists(coordsold, coords, natoms, atomlist, natomlist, drmax, rebuild, periodic, boxl)
   integer(kind=8), intent(in) :: natoms, natomlist, atomlist(natomlist)
   double precision, intent(in) :: coordsold(3*natoms), coords(3*natoms), drmax, boxl
   logical, intent(in) :: periodic
   logical, intent(out) :: rebuild
   double precision dr2max, dr2, dr(3)
   integer(kind=8) j1, j2, i1
   rebuild = .False.
   dr2max = drmax*drmax
   do j1=1,natomlist
      j2 = atomlist(j1) + 1 !+1 for fortran indices vs python indices
      i1 = 3*(j2-1) + 1
      dr2 = sum((coordsold(i1:i1+2) - coords(i1:i1+2))**2)
      if (dr2 .gt. dr2max) then
         if (periodic) then
            !try recalculating dr2 with periodic boundary conditions
            dr(:) = (coordsold(i1:i1+2) - coords(i1:i1+2))
            dr(:) = dr(:) - boxl*nint(dr(:) / boxl)
            dr2 = sum(dr**2)
            if (dr2 .le. dr2max) then
               !write(*,*) "    periodic. ok"
               cycle
            endif
         endif
         rebuild = .True.
         exit !exit do loop
      endif
   enddo
end subroutine check_neighbor_lists

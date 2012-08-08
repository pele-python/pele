subroutine build_neighbor_list1(coords, natoms, atomlist, natomlist, list, nlistmax, nlist, rlist2)
implicit none
integer, intent(in) :: natoms, natomlist, atomlist(natomlist), nlistmax
double precision, intent(in) :: coords(3*natoms), rlist2
integer, intent(out) :: list(nlistmax)
integer, intent(out) :: nlist
integer k1, k2, j1, j2
double precision r2
!note that atomlist is indexed as in python, e.g. 0,1,2,
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
integer, intent(in) :: natoms, nAlist, Alist(nAlist), nBlist, Blist(nBlist), nlistmax
double precision, intent(in) :: coords(3*natoms), rlist2
integer, intent(out) :: list(nlistmax)
integer, intent(out) :: nlist
integer k1, k2, j1, j2
double precision r2
!note that Alist is indexed as in python, e.g. 0,1,2,
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
integer, intent(in) :: natoms, natomlist, atomlist(natomlist), nlistmax
double precision, intent(in) :: coords(3*natoms), rlist2, boxl
integer, intent(out) :: list(nlistmax)
integer, intent(out) :: nlist
integer k1, k2, j1, j2
double precision r2, dr(3), iboxl
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
integer, intent(in) :: natoms, nAlist, Alist(nAlist), nBlist, Blist(nBlist), nlistmax
double precision, intent(in) :: coords(3*natoms), rlist2, boxl
integer, intent(out) :: list(nlistmax)
integer, intent(out) :: nlist
integer k1, k2, j1, j2
double precision r2, dr(3), iboxl
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

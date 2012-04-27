function overlap( xa, xb, l2, n, atomlist, nlist ) result(q)
implicit none
integer, intent(in) :: n, nlist, atomlist(nlist)
double precision, intent(in) :: xa(3*n), xb(3*n), l2
double precision dr2, q
integer j1, j2, i1, i2
q = 0.d0
do j1=1,nlist
    i1 = (atomlist(j1))*3
    do j2=1,nlist
        i2 = (atomlist(j2))*3
        dr2 = sum( (xa(i1+1:i1+3) - xb(i2+1:i2+3))**2 )
        q = q + exp( -dr2 / l2 )
    enddo
enddo
end function overlap

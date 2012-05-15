function overlap( xa, xb, L2, n, atomlist, nlist ) result(q)
implicit none
integer, intent(in) :: n, nlist, atomlist(nlist)
double precision, intent(in) :: xa(3*n), xb(3*n), L2
double precision dr2, q
integer j1, j2, i1, i2
q = 0.d0
do j1=1,nlist
    i1 = (atomlist(j1))*3
    do j2=1,nlist
        i2 = (atomlist(j2))*3
        dr2 = sum( (xa(i1+1:i1+3) - xb(i2+1:i2+3))**2 )
        q = q + exp( -dr2 / L2 )
    enddo
enddo
end function overlap

subroutine overlap_gradient( xa, xb, L2, n, atomlist, nlist, q, gradient )
implicit none
integer, intent(in) :: n, nlist, atomlist(nlist)
double precision, intent(in) :: xa(3*n), xb(3*n), L2
double precision, intent(out) :: q, gradient(3*n)
double precision dr2, dr(3), q0, iL2
integer j1, j2, i1, i2
iL2 = 1.d0/L2
q = 0.d0
gradient(:) = 0.d0
do j1=1,nlist
    i1 = (atomlist(j1))*3
    do j2=1,nlist
        i2 = (atomlist(j2))*3
        dr = xa(i1+1:i1+3) - xb(i2+1:i2+3)
        dr2 = sum( dr**2 )
        q0 = exp( -dr2 * iL2 )
        q = q + q0
        !gradient(i1+1 : i1+3) = gradient(i1+1 : i1+3) + (q0 * 2.d0 / L2) * dr(:)
        gradient(i2+1 : i2+3) = gradient(i2+1 : i2+3) + (q0 * 2.d0 * iL2) * dr(:)
    enddo
enddo
end  subroutine overlap_gradient

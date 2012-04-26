function OVERLAP( XA, XB, L2, N ) result(q)
IMPLICIT NONE
integer, intent(in) :: N
DOUBLE PRECISION, INTENT(IN) :: XA(3*N), xb(3*N), L2
double precision dr2, q
integer j1, j2, i1, i2
q = 0.D0
do j1=1,N
    i1 = (j1-1)*3
    do j2=1,N
        i2 = (j2-1)*3
        dr2 = sum( (XA(i1+1:i1+3) - XB(i2+1:i2+3))**2 )
        q = q + exp( -dr2 / l2 )
    enddo
enddo
end function overlap

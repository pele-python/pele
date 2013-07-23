function dprand() result(r)
!impliment a crappy random number generator.  To be changed later
implicit none
double precision r
r = dble(rand())
end function dprand

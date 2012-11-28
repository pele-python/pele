function dot(x1, x2, N) result(d)
   implicit none
   integer, intent(in) :: N
   double precision, intent(IN) :: x1(N)
   double precision, intent(IN) :: x2(N)
   double precision d
   integer j1
   d = 0.d0
   do j1=1,N
      d = d + x1(j1) * x2(j1)
   enddo
end function dot



SUBROUTINE NEB_FORCE(t, greal,  gspring, k, N, dneb, E, g_tot)
   implicit none
   integer, intent(in) :: N
   double precision, intent(IN) :: t(N)
   double precision, intent(IN) :: greal(N)
   double precision, intent(IN) :: gspring(N)
   double precision, intent(IN) :: k
   logical, intent(IN) :: dneb
   double precision, intent(OUT) :: g_tot(N), E
   double precision gperp(N)
   double precision gs_par(N)
   double precision gs_perp(N)
   double precision dot
        ! project out parallel part
        gperp = greal - dot(greal, t, N) * t
        ! the parallel part
        gs_par = dot(gspring, t, N) * t
        ! perpendicular part
        gs_perp = gspring - gs_par
                                
        g_tot = gperp + gs_par

        if (dneb) then
            ! double nudging
            g_tot = g_tot + gs_perp - dot(gs_perp, gperp, N) * gperp / dot(gperp, gperp, N)
        endif

        E = 0.5d0 * dot(gspring, gspring, N) / k
        
END SUBROUTINE NEB_FORCE

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
        ! project out parallel part
        gperp = greal - dot_product(greal, t) * t
        ! the parallel part
        gs_par = dot_product(gspring, t) * t
                                
        g_tot = gperp + gs_par

        if (dneb) then
            ! perpendicular part
            gs_perp = gspring - gs_par
            ! double nudging
            g_tot = g_tot + gs_perp - dot_product(gs_perp, gperp) * gperp / dot_product(gperp, gperp)
        endif

        E = 0.5d0 * dot_product(gspring, gspring) / k
        
END SUBROUTINE NEB_FORCE

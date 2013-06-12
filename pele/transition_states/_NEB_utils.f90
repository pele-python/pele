subroutine neb_force(t, greal, d_left, g_left, d_right, g_right, k, N, dneb, E, g_tot)
    implicit none
    integer, intent(in) :: N
    double precision, intent(IN) :: t(N)
    double precision, intent(IN) :: greal(N)
    double precision, intent(IN) :: g_left(N), g_right(N)
    double precision, intent(IN) :: d_left, d_right, k
    logical, intent(IN) :: dneb
    double precision, intent(OUT) :: g_tot(N), E
    double precision :: g_spring(N)
    double precision gperp(N)
    double precision gs_par(N)
    double precision gs_perp(N)

    ! project out parallel part
    gperp = greal - dot_product(greal, t) * t

    ! the parallel part
    gs_par = k*(d_left - d_right)*t

    g_tot = gperp + gs_par

    if (dneb) then
        g_spring = k*(g_left + g_right)
        ! perpendicular part of spring
        gs_perp = g_spring - dot_product(g_spring,t)*t
        ! double nudging
        g_tot = g_tot + gs_perp - dot_product(gs_perp,gperp)*gperp/dot_product(gperp,gperp)
    endif

    E = 0.5d0 * dot_product(g_spring, g_spring) / k

end subroutine neb_force

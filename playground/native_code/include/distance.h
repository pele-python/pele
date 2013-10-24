struct cartesian_distance {
    inline void get_rij(double * __restrict__ r_ij, 
                 double const * __restrict__ const r1, 
                 double const * __restrict__ const r2) 
    {
        r_ij[0] = r1[0] - r2[0];
        r_ij[1] = r1[1] - r2[1];
        r_ij[2] = r1[2] - r2[2];
    } 
};


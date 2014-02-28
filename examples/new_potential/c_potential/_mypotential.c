#include <math.h>
#include <stdio.h>

// TODO: make a lennard jones class! for now just do a function for testing


double mypotential(double *x, int N, double *grad, double eps, double sig)
{
  int i, j, k;
  double energy = 0;
  double sig12, sig24;
  double g, r, r2, ir2, ir12, ir24;
  double dr[3];
  sig12 = sig * sig * sig;
  sig12 *= sig12;
  sig12 *= sig12;
  sig24 = sig12 * sig12;


  for(i=0; i<N; i+=3) 
    for(j=i+3; j<N; j+=3)
      r2 = 0;
      for(k=0; k<3; ++k) {
        r = x[i+k] - x[j+k];
        r2 += r*r;
      }
      ir2 = 1. / r2;
      ir12 = ir2 * ir2 * ir2;
      ir12 *= ir12;
      ir24 = ir12 * ir12;
      energy += 4. * eps * (sig24 * ir24 - sig12 * ir12);
      g = 4. * eps * (-24. * sig24 * ir24 + 12. * sig12 * ir12);
      for(k=0; k<3; ++k) {
        grad[i+k] += g * ir2 * dr[k];
        grad[j+k] -= g * ir2 * dr[k];
      }

  return energy;
}

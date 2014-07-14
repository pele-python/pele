#ifndef _MYPOTENTIAL_H_
#define _MYPOTENTIAL_H_
#include <math.h>
#include <stdio.h>

double mypotential(double *x, int N, double *grad, double eps, double sig)
{
  int i, j, k;
  double energy = 0;
  double sig12, sig24;
  double g, r, r2, ir2, ir12, ir24;
  double dr[3];
  //printf("entering mypotential N %d\n", N);
  //for(i=0; i<N; i+=1) {
    //printf("coords %f\n", x[i]);
  //}

  sig12 = sig * sig * sig;
  sig12 *= sig12;
  sig12 *= sig12;
  sig24 = sig12 * sig12;

  for(i=0; i<N; i+=3){
    for(j=i+3; j<N; j+=3){
      r2 = 0;
      for(k=0; k<3; ++k) {
        dr[k] = x[i+k] - x[j+k];
        r2 += dr[k] * dr[k];
      }
      ir2 = 1. / r2;
      ir12 = ir2 * ir2 * ir2;
      ir12 *= ir12;
      ir24 = ir12 * ir12;
      energy += 4. * eps * (sig24 * ir24 - sig12 * ir12);
      //printf("i j %d %d %g %g\n", i, j, ir24, ir12);
      g = 4. * eps * (-24. * sig24 * ir24 + 12. * sig12 * ir12);
      for(k=0; k<3; ++k) {
        grad[i+k] += g * ir2 * dr[k];
        grad[j+k] -= g * ir2 * dr[k];
      }
    }
  }

  //printf("exiting c mypotential %g\n", energy);
  return energy;
}
#endif

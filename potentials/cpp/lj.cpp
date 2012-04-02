#include "nparray.h"
#include <math.h>

// TODO: make a lennard jones class! for now just do a function for testing
// also the whole calculation is crappy and no convenience wrappers for arrays

int python_array_pointer(boost::python::numeric::array& p, double **data);

double lj(double *x1, double *x2)
{
  double d[3], r2=0;
  
  for(int i=0; i<3; ++i) {
    d[i]=x2[i] - x1[i];
    r2+=d[i]*d[i];
  }
  
  return 4.*(pow(1./r2,6) - pow(1./r2, 3));
}

double ljg(double * __restrict__ x1, double * __restrict__ x2, double * __restrict__ g)
{
  double r2=0;
  
  for(int i=0; i<3; ++i) {
    g[i]=x2[i] - x1[i];
    r2+=g[i]*g[i];
  }
  //r2+=g[0]*g[0] + g[1]*g[1] + g[2]*g[2];
  
  double r6 = 1./(r2*r2*r2);
  double r12 = r6*r6;
  
  for(int i=0; i<3; ++i)
    g[i]*=4.0*(12.*r12 -  6.*r6)/r2;
  
  return 4.*(r12 - r6);
}

double energy(boost::python::numeric::array& px)
{
  int N;
  double r[3];
  
  NPArray<1> x(px);
  N = x.size(0);
 
  double energy = 0;
  for(int i=0; i<N; i+=3) 
    for(int j=i+3; j<N; j+=3)
      energy+=lj(&x[i], &x[j]);
  return energy;
}

double gradient(boost::python::numeric::array& px, boost::python::numeric::array& pgrad)
{
  int N;
  double gtmp[3];
  
  NPArray<1> x(px);
  NPArray<1> g(pgrad);
 
  N=x.size(0);
  for(int i=0; i<N; ++i)
    g[i]=0.0;
  
  double energy = 0;
  for(int i=0; i<N; i+=3)
    for(int j=i+3; j<N; j+=3) {
      energy+=ljg(&x[i], &x[j], gtmp);
      for(int k=0; k<3; ++k) {
        g[i+k]+=gtmp[k];
        g[j+k]-=gtmp[k];	
      }
    }
  return energy;
}

BOOST_PYTHON_MODULE(ljcpp_)
{
  using namespace boost::python;
  import_array();
  boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
  def("energy", energy);
  def("gradient", gradient);
}

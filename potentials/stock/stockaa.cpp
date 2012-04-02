#include "nparray.h"
#include <math.h>

extern "C" {
void stockaa_(int *n, double *x, double *g, double *e);
}
double gradient(boost::python::numeric::array& px, boost::python::numeric::array& pgrad)
{
  int N;
  double energy; 
  NPArray<1> x(px);
  NPArray<1> g(pgrad);
 
  N=x.size(0)/3;
  stockaa_(&N,&x[0],&g[0],&energy);  
  return energy;
}

BOOST_PYTHON_MODULE(stockaa_)
{
  using namespace boost::python;
  import_array();
  boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
  def("gradient", gradient);
}

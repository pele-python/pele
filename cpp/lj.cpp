#include <boost/python/numeric.hpp>
#include <boost/python.hpp>
#include <numpy/arrayobject.h>
#include <math.h>

// TODO: make a lennard jones class! for now just do a function for testing
// also the whole calculation is crappy and no convenience wrappers for arrays

int python_array_pointer(boost::python::numeric::array& p, double **data);

// not optimized!!
double lj(double *x1, double *x2)
{
  double d[3], r=0;
  
  for(int i=0; i<3; ++i) {
    d[i]=x2[i] - x1[i];
    r+=d[i]*d[i];
  }
  r=sqrt(r);
  
  return 4.*(pow(1./r,12) - pow(1./r, 6));
}

double ljg(double *x1, double *x2, double *g)
{
  double d[3], r=0;
  
  for(int i=0; i<3; ++i) {
    d[i]=x2[i] - x1[i];
    r+=d[i]*d[i];
  }
  r=sqrt(r);
  for(int i=0; i<3; ++i)
    d[i]/=r;
  
  for(int i=0; i<3; ++i)
    g[i]=4.0*d[i]*(12.*pow(1./r,13) -  6.*pow(1./r,7));
  
  return 4.*(pow(1./r,12) - pow(1./r, 6));
}

double energy(boost::python::numeric::array& px)
{
  int N;
  double *x, r[3];
  
  N = python_array_pointer(px, &x);
   
  double energy = 0;
  for(int i=0; i<N; i+=3) 
    for(int j=i+3; j<N; j+=3)
      energy+=lj(&x[i], &x[j]);
  return energy;
}

double gradient(boost::python::numeric::array& px, boost::python::numeric::array& pgrad)
{
  int N;
  double *x, *g, gtmp[3];
  
  N = python_array_pointer(px, &x);
  python_array_pointer(pgrad, &g);
  
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

int python_array_pointer(boost::python::numeric::array& p, double **data)
{
  const int ndims = 1;
  // Get pointer to np array
  PyArrayObject* a = (PyArrayObject*)PyArray_FROM_O(p.ptr());
  if (a == NULL) {
    throw std::runtime_error("Could not get NP array.");
  }
  if (a->descr->elsize != sizeof(double)) {
    throw std::runtime_error("Must be double ndarray");
  }
  
  if (a->nd != ndims) {
    throw std::runtime_error("Wrong dimension on array.");
  }
  
  *data = (double*)a->data;
  return *(a->dimensions);
}

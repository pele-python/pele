#include <math.h>
#include <wales/nparray.h>

using namespace wales;

double g_eps11=1.0;
double g_eps22=1.0;
double g_eps12=1.0;
double g_sigma11=1.0;
double g_sigma22=1.0;
double g_sigma12=1.0;
double g_q1=1;
double g_q2=-1;

void get_mindist(double * __restrict__ x1, double * __restrict__ x2, double * __restrict__ d)
{
	for(int i=0; i<3; ++i)
		d[i] = x2[i] - x1[i];
}

class PairInteraction
{
public:
	PairInteraction(double eps, double sigma, double q12)
	: _eps(eps), _sigma(sigma), _q12(q12) {}

	PairInteraction(int offset1, int n1, double eps, double sigma, double q12)
	: _offset1(offset1), _offset2(offset1), _n1(n1), _n2(n1), _eps(eps), _sigma(sigma), _q12(q12) {}

	PairInteraction(int offset1, int n1, int offset2, int n2, double eps, double sigma, double q12)
	: _offset1(offset1), _offset2(offset2), _n1(n1), _n2(n2), _eps(eps), _sigma(sigma), _q12(q12) {}


	void setSpecies1(int offset1, int n1) { _offset1 = offset1; _n1 = n1; }
	void setSpecies2(int offset2, int n2) { _offset2 = offset2; _n2 = n2; }


	double CalcEnergy(double *x) {
		double energy=0;
		if(_offset1 == _offset2) {
			for(int i=_offset1; i<_offset1 + _n1; ++i)
				for(int j=i+1; j<_offset1 + _n1; ++j) {
					energy
					+= pair_energy(&x[3*i], &x[3*j]);
				}
		} else {
			for(int i=_offset1; i<_offset1 + _n1; ++i)
				for(int j=_offset2; j<_offset2 + _n2; ++j) {
					energy += pair_energy(&x[3*i], &x[3*j]);
				}
		}
		return energy;
	}

	double CalcGradient(double *x, double *g) {
		double energy=0, gtmp[3];
		if(_offset1 == _offset2) {
			for(int i=_offset1; i<_offset1 + _n1; ++i)
				for(int j=i+1; j<_offset1 + _n1; ++j) {
					energy += eval_pair_gradient(x, g, i, j);
				}
		} else {
			for(int i=_offset1; i<_offset1 + _n1; ++i)
				for(int j=_offset2; j<_offset2 + _n2; ++j) {
					energy += eval_pair_gradient(x, g, i, j);
				}
		}
		return energy;
	}

	double eval_pair_gradient(double *x, double *g, int i, int j) {
		double gtmp[3];
		double energy=pair_gradient(&x[3*i], &x[3*j], gtmp);
		for(int k=0; k<3; ++k) {
			g[3*i+k]+=gtmp[k];
			g[3*j+k]-=gtmp[k];
		}
		return energy;
	}

	double pair_energy(double * __restrict__ x1, double * __restrict__ x2)
	{
		double d[3], r2=0;
        get_mindist(x1, x2, d);
		for(int i=0; i<3; ++i)
			r2+=d[i]*d[i];
		return 4.*_eps*(pow(_sigma/r2,6) - pow(_sigma/r2, 3)) + _q12/sqrt(r2);
	}

	double pair_gradient(double * __restrict__ x1, double * __restrict__ x2, double * __restrict__ g)
	{
		double r2=0;

		get_mindist(x1, x2, g);

		for(int i=0; i<3; ++i)
			r2+=g[i]*g[i];

		double r6 = pow(_sigma,6)/(r2*r2*r2);
		double r12 = r6*r6;

		for(int i=0; i<3; ++i)
			g[i]*=4.0*_eps*(12.*r12 -  6.*r6)/r2 + _q12/r2;

		return 4.*_eps*(r12 - r6) + _q12/sqrt(r2);
	}

private:
	double _eps, _sigma;
	double _q12;

	int _offset1, _n1;
	int _offset2, _n2;
};

double energy(boost::python::numeric::array& px)
{
	int N;
	double r[3];

	NPArray<1> x(px);
	N = x.size(0)/3;

	PairInteraction aa(0, N/2, g_eps11, g_sigma11, g_q1*g_q1);
	PairInteraction bb(N/2, N/2, g_eps22, g_sigma22, g_q2*g_q2);
	PairInteraction ab(0, N/2, N/2, N/2, g_eps12, g_sigma12, g_q1*g_q2);

	double energy=0;
	energy += aa.CalcEnergy(&x[0]);
	energy += ab.CalcEnergy(&x[0]);
	energy += bb.CalcEnergy(&x[0]);

	return energy;
}

double gradient(boost::python::numeric::array& px, boost::python::numeric::array& pgrad)
{
	int N;
	double gtmp[3];

	NPArray<1> x(px);
	NPArray<1> g(pgrad);

	N=x.size(0)/3;
	for(int i=0; i<3*N; ++i)
		g[i]=0.0;

	PairInteraction aa(0, N/2, g_eps11, g_sigma11, g_q1*g_q1);
	PairInteraction bb(N/2, N/2, g_eps22, g_sigma22, g_q2*g_q2);
	PairInteraction ab(0, N/2, N/2, N/2, g_eps12, g_sigma12, g_q1*g_q2);

	double energy=0;
	energy += aa.CalcGradient(&x[0], &g[0]);
	energy += ab.CalcGradient(&x[0], &g[0]);
	energy += bb.CalcGradient(&x[0], &g[0]);

	return energy;
}

BOOST_PYTHON_MODULE(salt_)
{
	using namespace boost::python;
	import_array();
	boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
	def("energy", energy);
	def("gradient", gradient);
}

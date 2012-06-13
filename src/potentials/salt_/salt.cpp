#include <math.h>
#include <wales/nparray.h>
#include "vec.h"
#include "matrix.h"
#include <vector>
#include "ewald.h"

#define USE_EWALD

using namespace wales;

double g_eps11=0.;//16.8*123.984/1000./96.;
double g_eps22=0.;//14.9*123.984/1000./96.;
double g_eps12=sqrt(g_eps11*g_eps22);
double g_sigma11=1.98;
double g_sigma22=5.24;
double g_sigma12=0.5*(g_sigma11+g_sigma22);
double g_q1=1.18; //;*sqrt(17.5);
double g_q2=-1.18;
double _cutoff=1.0;

using namespace votca::tools;

double entier(double x)
{
	double e = int(x);
	if (x-e < 0.) e -= 1.;
	if (x-e >= 1.) e += 1.;
	return e;
}

class System {
public:
	System(int n) {
		_pos.resize(n);
		_grad.resize(n);
		_charges.resize(n);
		for(int i=0; i<n/2; ++i) {
			_charges[i] = g_q1;
			_charges[i+n/2] = g_q2;
		}

		//_ewald.setTolerance(1.0, 1e-5);
		_ewald.setAlpha(1./0.320163);
	}

	void setLattice(double lattice[]) {
		_box=0.;
		for(int i=0; i<3; ++i)
			_box[i][i] = lattice[i];
		_box[1][0]=lattice[3];
		_box[2][0]=lattice[4];
		_box[2][1]=lattice[5];
		_boxinv = _box;
		_boxinv.Invert();
		buildImages();
		_ewald.setBox(getA(), getB(), getC());
		//printf("I have %d ghosts\n", _images.size());
		return;
		vec a = _box.getCol(0); vec b = _box.getCol(1); vec c = _box.getCol(2);

		std::cout << "The lattice matrix is: "
				  << _box << std::endl;
		std::cout << "a: " << a << std::endl;
		std::cout << "b: " << b << std::endl;
		std::cout << "c: " << c << std::endl;
		std::cout << "z: " << _box*vec(0, 0, 1) << std::endl;
	}

	vec getA() {
		return _box.getCol(0);
	}
	vec getB() {
		return _box.getCol(1);
	}
	vec getC() {
		return _box.getCol(2);
	}

	void setCoordinates(double *coords) {
		setLattice(&coords[3*_pos.size()]);
		for(int i=0; i<_pos.size(); ++i) {
			pos(i) = vec(&coords[3*i]);
			pos(i).x() -= entier(pos(i).getX());
			pos(i).y() -= entier(pos(i).getY());
			pos(i).z() -= entier(pos(i).getZ());
			//printf("a %f %f %f\n", pos(i).x(),pos(i).y(),pos(i).z());
			pos(i) = _box*pos(i);
			//printf("%f %f %f\n", pos(i).x(),pos(i).y(),pos(i).z());
		}
	}

	void getGradient(double *g) {
		matrix m=_box;
		m.Transpose();
		for(int i=0; i<_pos.size(); ++i) {
			vec tmp = m*grad(i);
			g[3*i+0] = tmp.getX();
			g[3*i+1] = tmp.getY();
			g[3*i+2] = tmp.getZ();
			//printf("gr %f %f %f\n", tmp.getX(),tmp.getY(),tmp.getZ());
		}

		for(int i=0; i<6; ++i)
			g[3*_pos.size()+i] = _glatt[i];
	}

	double addPressure() {
		return 0.;
		double volume = (getA()^getB())*getC();
		double p = 0.1;
		double gp[6];
		gp[0] = (getB()^getC()).getX();
		gp[1] = -(getA()^getC()).getY();
		gp[2] = (getA()^getB()).getZ();
		gp[3] = -(getB()^getC()).getY();
		gp[4] = (getB()^getC()).getZ();
		gp[5] = -(getA()^getC()).getZ();
		for(int i=0; i<6; ++i) {
			_glatt[i] += p*gp[i];
		}

		return p*volume;
	}
	void getReal(double *x) {
		for(int i=0; i<_pos.size(); ++i) {
			x[3*i+0] = pos(i).getX();
			x[3*i+1] = pos(i).getY();
			x[3*i+2] = pos(i).getZ();
		}
	}

	void resetGrad() {
		for(int i=0; i<_grad.size(); ++i)
			grad(i)=vec(0., 0, 0);
		for(int i=0; i<6; ++i)
			_glatt[i]=0.0;
	}

	vec mindist(const vec &r_i, const vec &r_j) const
	{
		//return r_j - r_i;
	    vec r_tp, r_dp, r_sp, r_ij;
	    vec a = _box.getCol(0); vec b = _box.getCol(1); vec c = _box.getCol(2);
	    r_tp = r_j - r_i;
	    r_dp = r_tp - c*round(r_tp.getZ()/c.getZ());
	    r_sp = r_dp - b*round(r_dp.getY()/b.getY());
	    r_ij = r_sp - a*round(r_sp.getX()/a.getX());
	    return r_ij;
	}

	vec &pos(int i) { return _pos[i]; }
	vec &grad(int i) { return _grad[i]; }

	std::vector<vec> &pos() { return _pos; }
	std::vector<double> &charges() { return _charges; }

	void addGLatt(const vec &g, const vec &r_ij)
	{
		vec u = -(_boxinv*r_ij);
		_glatt[0] += g.getX()*u.getX();
		_glatt[1] += g.getY()*u.getY();
		_glatt[2] += g.getZ()*u.getZ();
		_glatt[3] += g.getY()*u.getX();
		_glatt[4] += g.getZ()*u.getX();
		_glatt[5] += g.getZ()*u.getY();
		return;
	}

	void buildImages() {
		vec a = getA();
		vec b = getB();
		vec c = getC();
		//_images.clear();
		//_images.push_back(vec(0,0,0));
		//return;
		int n=1;
		for(int i=-n; i<=n; ++i)
			for(int j=-n; j<=n; ++j)
				for(int k=-n; k<=n; ++k) {
					vec v = i*a + j*b + k*c;
//	if(abs(v) < 2.*_cutoff) {
						_images.push_back(v);
					//}
				}
	}

	Ewald &ewald() { return _ewald; }

	std::vector<vec> &images() { return _images; }
private:
	matrix _box, _boxinv;
	std::vector<vec> _pos;
	std::vector<double> _charges;
	std::vector<vec> _grad;
	std::vector<vec> _images;
	double _glatt[6];
	Ewald _ewald;
};

class PairInteraction
{
public:
	PairInteraction(System &sys, double eps, double sigma, double q12)
	: _sys(sys), _eps(eps), _sigma(sigma), _q12(q12) {}

	PairInteraction(System &sys, int offset1, int n1, double eps, double sigma, double q12)
	: _sys(sys), _offset1(offset1), _offset2(offset1), _n1(n1), _n2(n1), _eps(eps), _sigma(sigma), _q12(q12) {}

	PairInteraction(System &sys, int offset1, int n1, int offset2, int n2, double eps, double sigma, double q12)
	: _sys(sys), _offset1(offset1), _offset2(offset2), _n1(n1), _n2(n2), _eps(eps), _sigma(sigma), _q12(q12) {}


	void setSpecies1(int offset1, int n1) { _offset1 = offset1; _n1 = n1; }
	void setSpecies2(int offset2, int n2) { _offset2 = offset2; _n2 = n2; }


	double CalcEnergy() {
		double energy=0;
		_cut_shift = 4.*_eps*(pow(_sigma/_cutoff,12) - pow(_sigma/_cutoff, 6))
#ifdef USE_EWALD
				+ 0; //_sys.ewald().PairEnergy(_cutoff, _q12);
#else
				+ _q12/_cutoff;
#endif
		_cut_shift = 0.;
		if(_offset1 == _offset2) {
			for(int i=_offset1; i<_offset1 + _n1; ++i)
				for(int j=i; j<_offset1 + _n1; ++j) {
					energy
					+= eval_pair_energy(i, j);
				}
		} else {
			for(int i=_offset1; i<_offset1 + _n1; ++i)
				for(int j=_offset2; j<_offset2 + _n2; ++j) {
					energy += eval_pair_energy(i, j);
				}
		}

		return energy;
	}

	double CalcGradient() {
		double energy=0;
		_cut_shift = 4.*_eps*(pow(_sigma/_cutoff,12) - pow(_sigma/_cutoff, 6)) + _sys.ewald()._f*_q12/_cutoff;
		_cut_shift = 0.;
		if(_offset1 == _offset2) {
			for(int i=_offset1; i<_offset1 + _n1; ++i)
				for(int j=i; j<_offset1 + _n1; ++j) {
					energy += eval_pair_gradient(i, j);
				}
		} else {
			for(int i=_offset1; i<_offset1 + _n1; ++i)
				for(int j=_offset2; j<_offset2 + _n2; ++j) {
					energy += eval_pair_gradient(i, j);
				}
		}
		return energy;
	}

	double eval_pair_gradient(int i, int j) {
		vec gtmp;
		vec d = _sys.pos(i) - _sys.pos(j); //_sys.mindist(_sys.pos(i), _sys.pos(j));
		double energy = 0;
		for(int l=0; l<_sys.images().size(); ++l) {
		energy+=pair_gradient(d+_sys.images()[l], gtmp);
			vec v =_sys.images()[l];
			//printf("%d %d %f %f %f\n", i, j, d.getX(), d.getY(), d.getZ());
			//if(i==j && abs(gtmp) > 1e-8) cout << i << " " << v << " " << gtmp << endl;
			_sys.grad(i)+=gtmp;
			_sys.grad(j)-=gtmp;
			//if(i==j) gtmp=2.*gtmp;
			_sys.addGLatt(gtmp, d+v);
		}
		return energy;
	}

	double eval_pair_energy(int i, int j) {
		vec d = _sys.mindist(_sys.pos(i), _sys.pos(j));
		double energy=0;
		for(int l=0; l<_sys.images().size(); ++l) {
			energy+=pair_energy(d+_sys.images()[l]);
		}
		return energy;
	}

	double pair_energy(const vec &d)
	{
		double r2;
        r2=d*d;
		if(r2>_cutoff*_cutoff || r2 < 1e-8) {
			return 0.;
		}
		return 4.*_eps*(pow(_sigma*_sigma/r2,6) - pow(_sigma*_sigma/r2, 3)) - _cut_shift
#ifdef USE_EWALD
				+ _sys.ewald().PairEnergy(sqrt(r2), _q12);
#else
				+ _q12/sqrt(r2);
#endif
	}

	double pair_gradient(const vec &r_ij, vec &g)
	{
		g = r_ij;
		double r2 = g*g;
		double r=sqrt(r2);
//cout << r << endl;
		if(r>_cutoff || r2 < 1e-8 ) {
			g=vec(0., 0., 0.);
			return 0.;
		}
		//printf("%f \n", r2);
		double r6 = pow(_sigma,6)/(r2*r2*r2);
		double r12 = r6*r6;

		g*=4.0*_eps*(12.*r12 -  6.*r6)/r2 + _q12/(r2*r);

	//	printf("%f %f %f\n", r, _q12/r , _cut_shift);
		return 4.*_eps*(r12 - r6) - _cut_shift
//#ifdef USE_EWALD
//				+ _sys.ewald().PairEnergy(r, _q12);
//#else
				+ _sys.ewald()._f*_q12/r;
//#endif
	}

private:
	double _eps, _sigma;
	double _q12;

	int _offset1, _n1;
	int _offset2, _n2;
	double _cut_shift;

	System &_sys;
};

double energy(boost::python::numeric::array& px)
{
	int N;
	double r[3];

	NPArray<1> x(px);
	N = x.size(0)/3-2;

	System sys(N);
	sys.setCoordinates(&x[0]);

	PairInteraction aa(sys, 0, N/2, g_eps11, g_sigma11, g_q1*g_q1);
	PairInteraction bb(sys, N/2, N/2, g_eps22, g_sigma22, g_q2*g_q2);
	PairInteraction ab(sys, 0, N/2, N/2, N/2, g_eps12, g_sigma12, g_q1*g_q2);

	double energy=0;
	energy += aa.CalcEnergy();
	energy += ab.CalcEnergy();
	energy += bb.CalcEnergy();
	//energy += sys.addPressure();
#ifdef USE_EWALD
	Ewald::zip_vectors adapter(sys.pos(), sys.charges());
	int ii=0;
	double dd = energy;
	for(Ewald::zip_vectors::iterator i= adapter.begin();
			i != adapter.end(); ++i) {
		//cout << "bla" << (*i)->getPos() << " " << (*i)->getQ() << endl
		//		<< sys.pos(ii) << " " << sys.charges()[ii] << endl;
		//++ii;
	}
	energy += sys.ewald().EnergyKSpace(adapter);
	cout << "Pair " << dd << endl;

#endif

	return energy;
}

double gradient(boost::python::numeric::array& px, boost::python::numeric::array& pgrad)
{
	int N;
	double gtmp[3];

	NPArray<1> x(px);
	NPArray<1> g(pgrad);

	N=x.size(0)/3-2;

	System sys(N);
	sys.setCoordinates(&x[0]);
	sys.resetGrad();
	PairInteraction aa(sys, 0, N/2, g_eps11, g_sigma11, g_q1*g_q1);
	PairInteraction bb(sys, N/2, N/2, g_eps22, g_sigma22, g_q2*g_q2);
	PairInteraction ab(sys, 0, N/2, N/2, N/2, g_eps12, g_sigma12, g_q1*g_q2);

	double energy=0;
	energy += aa.CalcGradient();
	cout << energy << endl;
	energy += ab.CalcGradient();
	cout << energy << endl;
	energy += bb.CalcGradient();
	cout << energy << endl;
	//energy += sys.addPressure();
	sys.getGradient(&g[0]);
	return energy;
}

void toReal(boost::python::numeric::array& px)
{
	int N;

	NPArray<1> x(px);

	N=x.size(0)/3-2;

	System sys(N);
	sys.setCoordinates(&x[0]);
	sys.getReal(&x[0]);
}

BOOST_PYTHON_MODULE(salt_)
{
	using namespace boost::python;
	import_array();
	boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
	def("energy", energy);
	def("gradient", gradient);
	def("toReal", toReal);
}

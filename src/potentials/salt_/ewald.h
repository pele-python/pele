#ifndef MINGMIN__EWALD_H
#define MINGMIN__EWALD_H

#include <complex>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/tuple/tuple.hpp>
#include "vec.h"

using namespace votca::tools;

/**
 * \brief Class for ewald electrostatics
 *
 * The class calculates the k-space part of the electrostatic
 * interactions using ewald summation. Iterating over pairs for
 * the real space part has to be done by user, summing over
 * PairEnergy for all pairs.
 *
 *
 */
class Ewald
{
public:
	Ewald();

	virtual ~Ewald() {};

	/**
	 * Calculate the K-space contribution of the energy
	 *
	 * \f[
	 *   E^{(k)} = \frac{1}{2V} \sum\limits_{k \neq 0} \frac{4\pi}{k^2}e^{-k^2/4\alpha^2} | S({\bf k}) |^2
	 * \f]
	 *
	 * \param charges container/adapter for charge iteration
	 */
	template<typename container>
	double EnergyKSpace(container &charges) const;

	/**
	 * Calculate real space energy for a pair
	 */
	double PairEnergy(double r, double q) const;

	/**
	 *
	 * Calculate reals space gradient for a pair
	 */
	//double pair_gradient(const vec r_ij, vec &g) const;

	/**
	 * set the box size
	 *
	 * \param a box vector a
	 * \param b box vector b
	 * \param c box vector c
	 */
	void setBox(vec a, vec b, vec c);
	void setAlpha(double alpha) { _alpha = alpha; }

	/***
	 *	\brief set the permitivity in MD units.
	 *
	 *	Sets the permitivity in MD units. The default
	 *	is in kJ nm / mol e^2 which means distances are
	 *	in nm, charges in e and energies in kJ/mol
	 *
	 *	\parameter eps0 permitivity in MD units
	 *
	 */void setPermitivity(double eps0=5.727656e-4) { _f = 0.25 / (M_PI*eps0); }
	//void setTolerance(double rcut, double rtol) { _alpha = erfcinv(rtol)/rcut; }

	class zip_vectors {
	public:
		zip_vectors(std::vector<vec> &positions, std::vector<double> &charges)
			: _positions(positions), _charges(charges) {}

		struct charge_t {
			charge_t(vec &pos, double &q) : _pos(pos), _q(q) {}
			const vec &getPos() { return _pos; }
			double getQ() { return _q; }
			charge_t *operator->() { return this; }
		private:
			double &_q;
			vec &_pos;
		};

		struct iterator {
			iterator();
			iterator(std::vector<vec>::iterator ipos, std::vector<double>::iterator iq)
				: _ipos(ipos), _iq(iq) {}
			iterator(const iterator &i) : _ipos(i._ipos), _iq(i._iq) {}

			charge_t operator*() {
				return charge_t(*_ipos, *_iq);
			}

			iterator &operator++() {
				++_ipos; ++_iq; return *this;
			}

			iterator &operator=(const iterator &i) {
				_ipos = i._ipos; _iq = i._iq; return *this;
			}

			bool operator!=(const iterator &i) { return i._ipos != _ipos; }
		private:
			std::vector<vec>::iterator _ipos;
			std::vector<double>::iterator _iq;
		};

		iterator begin() { return iterator(_positions.begin(), _charges.begin()); }
		iterator end() { return iterator(_positions.end(), _charges.end()); }

	private:
		std::vector<vec> &_positions;
		std::vector<double> &_charges;
	};
	double _f;

protected:
	double _alpha;

	// reciprocal lattice vectors
	vec _ar, _br, _cr;
	double _lmax, _mmax, _nmax;
	// box volume
	double _V;

	/**
	 * Calculate structure factor
	 *
	 * \f[
	 *   S({\bf k}) = \sum\limits_{j=1}^N q_j e^{i{\bf k} \dot {\bf r_j}}
	 * \f]
	 */
	template<typename container>
	std::complex<double> S(const vec &k, container &charges) const;

	template<typename container>
	double Dipolar(container &charges) const;

	/**
	 * Calculate self energy
	 *
	 * \f[
	 *   E^{(s)} = - \frac{\alpha}{\sqrt{\pi}} \sum\limits_i q_i^2
	 * \f]
	 */
	template<typename container>
	double SelfEnergy(container &charges) const;
};

inline Ewald::Ewald()
	: _alpha(1.), _V(0.0)
{
	setPermitivity();
}

inline void Ewald::setBox(vec a, vec b, vec c)
{
	_V = abs((a^b)*c);
	_ar = 2.*M_PI*b^c / _V;
	_br = 2.*M_PI*c^a / _V;
	_cr = 2.*M_PI*a^b / _V;

	_lmax = _mmax = _nmax = 20;
}


inline double Ewald::PairEnergy(double r, double q) const
{
	double z = _alpha * r;
	return _f*q/r*erfc(z);
}

template<typename container>
inline double Ewald::EnergyKSpace(container &charges) const
{
	double E = 0;

	for(int l=-_lmax; l<_lmax; ++l) {
		for(int m=-_mmax; m<_mmax; ++m) {
			for(int n=-_nmax; n<_nmax; ++n) {
				if(l == 0 && m == 0 && n == 0) continue;
				vec k = double(l)*_ar + double(m)*_br + double(n)*_cr;
				double k2 = k*k;
				std::complex<double> s=S(k, charges);
				double S2 = s.real()*s.real() + s.imag()*s.imag();
				E+= exp(-0.25*k2/(_alpha*_alpha)) / k2 * S2;
				//std::cout << "E " << E << " " << k2 << k << std::endl;
			}
		}
	}
	E = _f * 2.*M_PI/_V*E ;
	printf("\nKSPace %f\n", E);
	printf("Self %f\n", SelfEnergy(charges));
	printf("Recip %f\n", E - SelfEnergy(charges));
//	printf("Dipolar %f\n", Dipolar(charges));
	return E - SelfEnergy(charges); //+ Dipolar(charges);
}


template<typename container>
inline std::complex<double> Ewald::S(const vec &k, container &charges) const
{
	std::complex<double> S(0.,0.);
	for(typename container::iterator crg=charges.begin();
			crg!=charges.end();
			++crg) {
		S+=(*crg)->getQ()*exp(complex<double>(0.,-k*(*crg)->getPos()));
	}
	//std::cout << S << std::endl;
	return S;
}

template<typename container>
inline double Ewald::SelfEnergy(container &charges) const
{
	double E=0;
	for(typename container::iterator crg=charges.begin();
		crg!=charges.end(); ++crg) {
		double q = (*crg)->getQ();
		E+=q*q;
//std::cout << (*crg)->getPos() << " " << q << endl;
	}
	return _f*E*_alpha/sqrt(M_PI);
}

template<typename container>
double Ewald::Dipolar(container &charges) const
{
	double E=0;
	vec d(0.,0.,0.);
	for(typename container::iterator crg=charges.begin();
		crg!=charges.end(); ++crg) {
		d += (*crg)->getQ()*(*crg)->getPos();
//std::cout << (*crg)->getPos() << " " << q << endl;
	}
	return 2.*M_PI/(3. * _V)*d*d;
}

#endif

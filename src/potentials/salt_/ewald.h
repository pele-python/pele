#ifndef MINGMIN__EWALD_H
#define MINGMIN__EWALD_H

#include <complex>
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

protected:
	double _alpha;
	double _sigma;

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
	: _alpha(1.), _V(0.0) {}

inline void Ewald::setBox(vec a, vec b, vec c)
{
	_V = abs((a^b)*c);
	_ar = b^c / _V;
	_br = c^a / _V;
	_cr = a^b / _V;
}


inline double Ewald::PairEnergy(double r, double q) const
{
	double z = _alpha * r;
	return q/r*erfc(z);
}

template<typename container>
inline double Ewald::EnergyKSpace(container &charges) const
{
	vec a,b,c;
	double E = 0;
	double V;

	for(int l=0; l<_lmax; ++l) {
		for(int m=0; m<_mmax; ++m) {
			for(int n=0; n<_nmax; ++n) {
				vec k = double(l)*a + double(m)*b + double(n)*c;
				double k2 = k*k;
				std::complex<double> S=S(k, charges);
				double S2 = S.real()*S.real() + S.imag()*S.imag();
				E+= exp(-0.25*k2/(_alpha*_alpha)) / k2 * S2;
			}
		}
	}
	E = 2.*M_PI/V*E - SelfEnergy(charges);
}


template<typename container>
inline std::complex<double> Ewald::S(const vec &k, container &charges) const
{
	std::complex<double> S(0.,0.);
	for(typename container::iterator crg=charges.begin();
			crg!=charges.end();
			++crg) {
		S.real()+=crg.q()*exp(complex<double>(0.,-k*crg.r()));
	}
	return S;
}

template<typename container>
inline double Ewald::SelfEnergy(container &charges) const
{
	double E=0;
	for(typename container::iterator crg=charges.begin();
		crg!=charges.end(); ++crg) {
		E+=crg.q()*crg.q();
	}
	return E*_alpha/sqrt(M_PI);
}

#endif

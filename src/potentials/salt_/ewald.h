#ifndef MINGMIN__EWALD_H
#define MINGMIN__EWALD_H

#include <complex>
#include "vec.h"


/**
 * \brief Class for ewald electrostatics
 *
 * The class calculates the k-space part of the electrostatic
 * interactions using ewald summation. Iterating over pairs for
 * the real space part has to be done by user, summing over
 * PairEnergy for all pairs.
 */
class Ewald
{
public:
	Ewald();

	virtual ~Ewald() {};

	/**
	 * Calculate the K-space contribution of the energy
	 *
	 * \param charges container/adapter for charge iteration
	 */
	template<container>
	double EnergyKSpace(container &charges);

	/**
	 * Calculate real space energy for a pair
	 */
	double pair_energy(const vec3 &r_ij) const;

	/**
	 * Calculate reals space gradient for a pair
	 */
	double pair_gradient(const vec3 r_ij, vec3 &g) const;

	/**
	 * set the box size
	 *
	 * \param a box vector a
	 * \param b box vector b
	 * \param c box vector c
	 */
	void setBox(vec3 a, vec3 b, vec3 c);

protected:
	double _eps0;
	double _sigma;

	// reciprocal lattice vectors
	vec _ar, _br, _cr;
	// box volume
	double _V;

	/**
	 * Calculate structure factor
	 */
	template<container>
	std::complex<double> S(const vec &k, container &charges) const;

	template<container>
	double SelfEnergy(container &charges);
};

inline Ewald::Ewald()
	: _sigma(0.), _eps0(1.0), _V(0.0) {}

inline void Ewal::setBox(vec3 a, vec3 b, vec3 c)
{
	_V = abs((a^b)*c);
	_ar = b^c / _V;
	_br = c^a / _V;
	_cr = a^b / _V;
}


inline double Ewald::PairEnergy(double r, double q)
{
	double z = r/(sqrt(2.)*_sigma);
	return 1./(4.*M_PI*_eps0)*q/r*erfc(z);
}

template<container>
inline double Ewald::EnergyKSpace(container &charges)
{
	vec a,b,c;
	double E = 0;
	double V;

	for(int l=0; l<lmax; ++l) {
		for(int m=0; m<mmax; ++m) {
			for(int n=0; n<nmax; ++n) {
				vec k = double(l)*a + double(m)*b + double(n)*c;
				double k2 = k*k;
				std::complex<double> S=S(k, charge);
				double S2 = S.real()*S.real() + S.imag()*S.imag();
				E+= exp(-0.5*_sigma*_sigma *k2) * S2 / k2;
			}
		}
	}
	E = E/(2.*V*_eps0) - SelfEnergy(charges);
}


template<container>
inline std::complex<double> Ewald::S(const vec &k, container &charges) const
{
	std::complex<double> S(0.,0.);
	for(container::iterator crg=charges.begin();
			crg!=charges.end();
			++crg) {
		S.real()+=crg.q()*exp(complex<double(0.,k*crg.r()));
	}
	return S;
}

template<container>
inline double SelfEnergy(container &charges)
{
	double E=0;
	for(container::iterator crg=charges.begin();
		crg!=charges.end(); ++crg) {
		E+=crg.q()*crg.q();
	}
	return E/(3.*sqrt(M_PI)*_eps0*_sigma)/sqrt(2.*M_PI);
}
#endif

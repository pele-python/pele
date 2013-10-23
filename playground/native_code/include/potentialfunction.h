#include "array.h"
#include "base_potential.h"

namespace pele {

/**
 * This class wraps a get_energy function and a get_energy_gradient function in
 * a class that derives from BasePotential.  This is necessary to be able to use
 * the functions in the pele c++ interface.  This is a backup method, the
 * preferred method is to define a class separately for each potential, which
 * eliminates the need to pass around the void * userdata parameter.
 */
class PotentialFunction : public BasePotential
{
public:
	typedef double EnergyCallback(Array, void *);
	typedef double EnergyGradientCallback(Array, Array, void *);

	PotentialFunction(EnergyCallback *get_energy, EnergyGradientCallback *get_energy_gradient, void *userdata)
		:	_get_energy(get_energy), _get_energy_gradient(get_energy_gradient), _userdata(userdata) {}

	virtual double get_energy(Array x) { return (*_get_energy)(x, _userdata); } ;
	virtual double get_energy_gradient(Array x, Array grad) {  return (*_get_energy_gradient)(x, grad, _userdata); }

private:
	EnergyCallback *_get_energy;
	EnergyGradientCallback *_get_energy_gradient;
	void *_userdata;
};


}

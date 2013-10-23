#include "include/array.h"
#include "include/base_potential.h"

namespace pele {

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

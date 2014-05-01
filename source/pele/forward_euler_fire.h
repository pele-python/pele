#ifndef _PELE_FORWARD_EULER_H__
#define _PELE_FORWARD_EULER_H__

#include "base_integrator.h"

namespace pele {

    /**
       * An implementation of the forward euler algorithm in c++ *specific*
       * to the modifiedFIRE algorithm.
       */

    class ForwardEuler: public BaseIntegrator
      {
      public:
          /*Constructor*/
        ForwardEuler(pele::BasePotential * potential, pele::Array<double> x, double dt, double maxstep,
                    pele::Array<double> v = pele::Array<double>(), pele::Array<double> g = pele::Array<double>(),
                    pele::Array<double> m = pele::Array<double>());

        inline void oneiteration()
              {
                   /* the minuses in the following expressions are due to the fact that
                   * the gradients rather than the forces appear in the expression
                   */
                  size_t i;
                  double normdx;

                  for(i =0; i < _x.size(); ++i) //this was after get_energy_gradient, moved for testing
                  {
                      _v[i] -= _dt * _g[i] / _m[i];     //update velocity
                      _dx[i] = _dt * _v[i];    //build displacement vector
                      _gold[i] = _g[i]; //save gradient as old g
                  }

                  normdx = norm(_dx);

                  if(normdx > _maxstep){
                    _dx *= (_maxstep / normdx); //resize displacement vector is greater than _maxstep
                    }

                  _x += _dx;

                  *_E = _potential->get_energy_gradient(_x, _g);    //update gradient
              }

        void run(size_t const N)
                {
                   for(size_t i =0; i < N; ++i)
                   {
                       oneiteration();
                   }
                }

      };

    ForwardEuler::ForwardEuler(pele::BasePotential * potential, pele::Array<double> x, double dt, double maxstep,
            pele::Array<double> v, pele::Array<double> g, pele::Array<double> m):
            BaseIntegrator(potential, x, dt, maxstep, v, g, m) //initialise base integrator from which this class is inherited
            {}

}
#endif

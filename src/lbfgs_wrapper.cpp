#include <stdio.h>
#include <lbfgs.h>
#include "potential.h"
#include <iostream>

namespace pygmin {

class LBFGS
{
protected:
	Potential *_pot;

public:
    LBFGS() : _pot(NULL) {}
    LBFGS(Potential *pot) : _pot(pot) { std::cout << "potential set" << std::endl;}

    ~LBFGS() {}

    int run(Array &x)
    {
    	double fx;
        int ret = lbfgs(x.size(), x.data(), &fx, _evaluate, _progress, this, NULL);
        return ret;
    }

protected:
    int progress(
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
    {
        /*printf("Iteration %d:\n", k);
        printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
        printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
        printf("\n");*/
        return 0;
    }

    lbfgsfloatval_t evaluate(
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
    		Array _x((double*)x, n), _g(g, n);
    		return _pot->get_energy_gradient(_x, _g);

    }

    static lbfgsfloatval_t _evaluate(
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
    	return reinterpret_cast<LBFGS*>(instance)->evaluate(x, g, n, step);
    }

    static int _progress(
        void *instance,
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
    {
        return reinterpret_cast<LBFGS*>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
    }
};
}


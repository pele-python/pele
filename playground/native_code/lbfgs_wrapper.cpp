#include <stdio.h>
#include <lbfgs.h>
#include "potential.h"
#include <iostream>
#include <string>

namespace pele {

class LBFGS
{
protected:
    Potential *_pot;
    lbfgs_parameter_t _params;
    int _error_code;

public:
    LBFGS() : _pot(NULL) { setup(); }
    LBFGS(Potential *pot) : _pot(pot) { setup(); }

    ~LBFGS() {}

    int run(Array &x)
    {
    	double fx;
        int ret = lbfgs(x.size(), x.data(), &fx, _evaluate, _progress, this, &_params);
        _error_code = ret;
        return ret;
    }

    int error_code() {
    	return _error_code;
    }

    const char *error_string() {
    	switch(_error_code) {
    	case LBFGS_SUCCESS: return "LBFGS_SUCCESS";
    	//case LBFGS_CONVERGENCE: return "LBFGS_CONVERGENCE";
    	case LBFGS_STOP: return "LBFGS_STOP";
    	case LBFGS_ALREADY_MINIMIZED: return "LBFGS_ALREADY_MINIMIZED";
    	case LBFGSERR_UNKNOWNERROR: return "LBFGSERR_UNKNOWNERROR";
    	case LBFGSERR_LOGICERROR: return "LBFGSERR_LOGICERROR";
    	case LBFGSERR_OUTOFMEMORY: return "LBFGSERR_OUTOFMEMORY";
    	case LBFGSERR_CANCELED: return "LBFGSERR_CANCELED";
    	case LBFGSERR_INVALID_N: return "LBFGSERR_INVALID_N";
    	case LBFGSERR_INVALID_N_SSE: return "LBFGSERR_INVALID_N_SSE";
    	case LBFGSERR_INVALID_X_SSE: return "LBFGSERR_INVALID_X_SSE";
    	case LBFGSERR_INVALID_EPSILON: return "LBFGSERR_INVALID_EPSILON";
    	case LBFGSERR_INVALID_TESTPERIOD: return "LBFGSERR_INVALID_TESTPERIOD";
    	case LBFGSERR_INVALID_DELTA: return "LBFGSERR_INVALID_DELTA";
    	case LBFGSERR_INVALID_LINESEARCH: return "LBFGSERR_INVALID_LINESEARCH";
    	case LBFGSERR_INVALID_MINSTEP: return "LBFGSERR_INVALID_MINSTEP";
    	case LBFGSERR_INVALID_MAXSTEP: return "LBFGSERR_INVALID_MAXSTEP";
    	case LBFGSERR_INVALID_FTOL: return "LBFGSERR_INVALID_FTOL";
    	case LBFGSERR_INVALID_WOLFE: return "LBFGSERR_INVALID_WOLFE";
    	case LBFGSERR_INVALID_GTOL: return "LBFGSERR_INVALID_GTOL";
    	case LBFGSERR_INVALID_XTOL: return "LBFGSERR_INVALID_XTOL";
    	case LBFGSERR_INVALID_MAXLINESEARCH: return "LBFGSERR_INVALID_MAXLINESEARCH";
    	case LBFGSERR_INVALID_ORTHANTWISE: return "LBFGSERR_INVALID_ORTHANTWISE";
    	case LBFGSERR_INVALID_ORTHANTWISE_START: return "LBFGSERR_INVALID_ORTHANTWISE_START";
    	case LBFGSERR_INVALID_ORTHANTWISE_END: return "LBFGSERR_INVALID_ORTHANTWISE_END";
    	case LBFGSERR_OUTOFINTERVAL: return "LBFGSERR_OUTOFINTERVAL";
    	case LBFGSERR_INCORRECT_TMINMAX: return "LBFGSERR_INCORRECT_TMINMAX";
    	case LBFGSERR_ROUNDING_ERROR: return "LBFGSERR_ROUNDING_ERROR";
    	case LBFGSERR_MINIMUMSTEP: return "LBFGSERR_MINIMUMSTEP";
    	case LBFGSERR_MAXIMUMSTEP: return "LBFGSERR_MAXIMUMSTEP";
    	case LBFGSERR_MAXIMUMLINESEARCH: return "LBFGSERR_MAXIMUMLINESEARCH";
    	case LBFGSERR_MAXIMUMITERATION: return "LBFGSERR_MAXIMUMITERATION";
    	case LBFGSERR_WIDTHTOOSMALL: return "LBFGSERR_WIDTHTOOSMALL";
    	case LBFGSERR_INVALIDPARAMETERS: return "LBFGSERR_INVALIDPARAMETERS";
    	case LBFGSERR_INCREASEGRADIENT: return "LBFGSERR_INCREASEGRADIENT";
    	}
    }

protected:
    void setup(void) {
        lbfgs_parameter_init(&_params);
        //_params.max_step=0.1;
        _params.min_step=1e-320;
    }

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


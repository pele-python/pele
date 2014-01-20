#ifndef _PELE_LBFGS_H__
#define _PELE_LBFGS_H__

#include <vector>
#include "base_potential.h"
#include "array.h"
#include "optimizer.h"

using std::vector;

namespace pele{
	/**
	 * An implementation of the LBFGS optimization algorithm in c++.  This
	 * Implementation uses a backtracking linesearch.
	 */
	class LBFGS : public GradientOptimizer{
	private:

		int M_; /**< The length of the LBFGS memory */
		double max_f_rise_; /**< The maximum the function is allowed to rise in a
							 * given step.  This is the criterion for the
							 * backtracking line search.
							 */

		// places to store the lbfgs memory
		std::vector<Array<double> > s_;
		std::vector<Array<double> > y_;
		Array<double> rho_;
		double H0_;
		int k_; /**< Counter for how many times the memory has been updated */

	public:
		/**
		 * Constructor
		 */
		LBFGS(
			pele::BasePotential * potential,
			const pele::Array<double> x0,
			double tol = 1e-4,
			int M = 4
			);

		/**
		 * Destructor
		 */
		~LBFGS() {}

		/**
		 * Do one iteration iteration of the optimization algorithm
		 */
		void one_iteration();

		// functions for setting the parameters
		void set_H0(double H0)
		{
			if (iter_number_ > 0){
				cout << "warning: setting H0 after the first iteration.\n";
			}
			H0_ = H0;
		}
		void set_max_f_rise(double max_f_rise) { max_f_rise_ = max_f_rise; }

		// functions for accessing the results
		double get_H0() { return H0_; }

	private:

		/**
		 * Add a step to the LBFGS Memory
		 * This updates s_, y_, rho_, H0_, and k_
		 */
		void update_memory(
			Array<double> xold,
			Array<double> gold,
			Array<double> xnew,
			Array<double> gnew);

		/**
		 * Compute the LBFGS step from the memory
		 */
		void compute_lbfgs_step(Array<double> step);

		/**
		 * Take the step and do a backtracking linesearch if necessary.
		 * Apply the maximum step size constraint and ensure that the function
		 * does not rise more than the allowed amount.
		 */
		double backtracking_linesearch(Array<double> step);

	};
}

#endif

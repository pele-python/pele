#ifndef PYGMIN_NEB_H
#define PYGMIN_NEB_H

#include <vector>
#include <list>
#include "potential.h"

namespace pele {

	/**
	 * Doubly-Nudged Elastic Band
	 */
	class NEB {
	public:
		class NEBDistance;

	protected:
		double _k;

		// number of coordinates per image
		int _N;
		// number of images
		int _nimages;

		// callback function to calculate distances
		NEBDistance *_neb_distance;

		// callback function to calculate distances
		Potential *_potential;

		// full neb coordinates
		Array _neb_coords;

		void initialize(void);

		// interpolate between 2 images
		void interpolate(Array &x1, Array &x2, Array &xout, double t);

	public:
		NEB(Potential *potential, NEBDistance *neb_distance) :
			_potential(potential), _neb_distance(neb_distance) { initialize(); };

		void set_nimages(int nimages);

		void set_path(Array &x1, Array &x2);
		void set_path(std::list< Array *> &path);

		/// run a full optimization of neb
		void run();

		/// run one neb cycle (step), returns True when done
		bool one_cycle();

		/***
		 *  interface for distance calculations
		 */
		class NEBDistance {
		public:
			virtual ~NEBDistance() {}
			virtual double dist(Array &x1, Array &x2) = 0;
		};
	};

	void test_array(Array a);
	void test_py(void);
}

#endif

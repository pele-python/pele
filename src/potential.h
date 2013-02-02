#ifndef PYGMIN_POTENTIAL_H
#define PYGMIN_POTENTIAL_H

#include <stdexcept>
#include <iostream>

namespace pygmin {
	/**
	 * Simple wrapper class for arrays
	 *
	 * The Array class should provide a simple array handling interface
	 * and allow to efficiently transfer data between C++/Fortran/Python
	 */
	class Array {
		double *_data;
		size_t _size;

		bool _owner;
	public:
		Array() : _data(NULL), _size(0), _owner(true) {}

		/// create array with specific size and allocate memory
		Array(size_t size) : _size(size) { _data = new double[size]; _owner=true; }

		/// allocate memory with existing data
		///
		/// if owner is true the data will be deleted in the destructor
		Array(double *data, size_t size, bool owner=false)
			: _data(data), _size(size), _owner(owner) {}

		Array(std::vector<double> &x) : _data(x.data()), _size(x.size()), _owner(false) {}

		~Array() { if(_owner) delete[] _data; _data = NULL; _size = 0; }

		/// return pointer to data
		double *data() { return _data; }
		/// return size of array
		size_t size() const { return _size; }

		/// resize the array
		void resize(size_t size) {
			if(!_owner)
				throw std::runtime_error("Array: cannot resize Arrays if not owner of data");
			delete [] _data;
			_size = size;
			_data = new double[_size];
		}

		/// access an element in the array
		double &operator[](size_t i) { return _data[i]; }

		/// read only access to element in array
		double operator()(size_t i) const { return _data[i]; }

		Array &operator=(double d) {
			for(int i=0; i<_size; ++i)
				_data[i] = d;
			return *this;
		}
	};

	// for array printing
	inline std::ostream &operator<<(std::ostream &out, const Array &a) {
		out << "[ ";
		for(int i=0; i<a.size();++i) {
			if(i>0) out << ", ";
			out << a(i);
		}
		out << " ]";
	}

	/***
	 * basic potential interface for native potentials
	 */
	class Potential {
	public:
		virtual ~Potential() {}

		virtual double get_energy(Array &x) = 0;
		virtual double get_energy_gradient(Array &x, Array &grad) = 0;
	};

}



#endif

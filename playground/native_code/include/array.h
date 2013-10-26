#ifndef PYGMIN_ARRAY_H
#define PYGMIN_ARRAY_H

#include <vector>
#include <stdexcept>
#include <iostream>

namespace pele {
	/**
	 * Simple wrapper class for arrays
	 *
	 * The Array class should provide a simple array handling interface
	 * and allow to efficiently transfer data between C++/Fortran/Python
	 */
  template<typename dtype>
	class Array {
		dtype *_data;
		size_t _size;

		bool _owner;
	public:
		Array() : _data(NULL), _size(0), _owner(true) {}

		/// create array with specific size and allocate memory
		Array(size_t size) : _size(size) { _data = new dtype[size]; _owner=true; }

		/// allocate memory with existing data
		///
		/// if owner is true the data will be deleted in the destructor
		Array(dtype *data, size_t size, bool owner=false)
			: _data(data), _size(size), _owner(owner) {}

        // wrap a vector
		Array(std::vector<dtype> &x) : _data(x.data()), _size(x.size()), _owner(false) { }

        // wrap another array
		//Array(Array<dtype> &x) : _data(x.data()), _size(x.size()), _owner(false) { }


        // destructor
		~Array() { if(_owner && _data != NULL) delete[] _data; _data = NULL; _size = 0; }

		/// return pointer to data
		dtype *data() { return _data; }
		dtype const *data() const { return _data; }

		/// return size of array
		size_t size() const { return _size; }

        // return iterators over data
        typedef dtype * iterator;
        typedef dtype const * const_iterator;
        iterator begin() { return &_data[0]; }
        iterator end() { return &_data[_size]; }
        const_iterator begin() const { return &_data[0]; }
        const_iterator end() const { return &_data[_size]; }


		/// resize the array
		void resize(size_t size) {
			if(!_owner)
				throw std::runtime_error("Array: cannot resize Arrays if not owner of data");
			delete [] _data;
			_size = size;
			_data = new dtype[_size];
		}

		/// access an element in the array
		dtype &operator[](size_t i) { return _data[i]; }
		dtype operator[](size_t i) const { return _data[i]; }

		/// read only access to element in array
		dtype operator()(size_t i) const { return _data[i]; }

		Array<dtype> &operator=(dtype d) {
			for(size_t i=0; i<_size; ++i)
				_data[i] = d;
			return *this;
		}
	};

	// for array printing
	inline std::ostream &operator<<(std::ostream &out, const Array<double> &a) {
		out << "[ ";
		for(size_t i=0; i<a.size();++i) {
			if(i>0) out << ", ";
			out << a(i);
		}
		out << " ]";
		return out;
	}
}

//Array newa(old);  Array newa; newa.view(old)
//Array newa = old; 
//Array newa(old.copy()); Array newa(old); 
//Array newa(old.view(0, 100)); Array newai; newa.view(old, 0, 100);
//
//newa.reference(old);
//newa.copy(old);

#endif

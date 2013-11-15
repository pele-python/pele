#ifndef PYGMIN_ARRAY_H
#define PYGMIN_ARRAY_H

#include <assert.h>
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
        dtype *_allocated_memory;
        size_t _size;

        long int *_reference_count;
    public:
        Array() : _data(NULL), _allocated_memory(NULL), _size(0), _reference_count(NULL) {}

        /** 
         * create array with specific size and allocate memory
         */
        Array(size_t size) : _size(size) 
        { 
            _allocated_memory = new dtype[size]; 
            _data = _allocated_memory;
            _reference_count = new long int;
            *_reference_count = 1;
        }

        /** 
         * wrap another array
         */
        Array(Array<dtype> const &x) : 
            _data(x._data), _allocated_memory(x._allocated_memory),
            _size(x._size), _reference_count(x._reference_count) 
        {
            if (_data == NULL){
                throw std::runtime_error("cannot wrap an array with no data");
            }
            assert((_reference_count==NULL) == (_allocated_memory==NULL)); //both null or both not null
            if (_reference_count != NULL){
                *_reference_count += 1;
            }
        }

        /**
         * wrap some data that is passed.  Do not take ownership of the data.
         */
        Array(dtype *data, size_t size)
            : _data(data), _allocated_memory(NULL), _size(size), _reference_count(NULL) {}

        /**
         * wrap a vector.  This memory should never be deleted.  
         */
        Array(std::vector<dtype> &x) : _data(x.data()), _allocated_memory(NULL), _size(x.size()), _reference_count(NULL) { }

        /** 
         * destructor
         */
        ~Array() 
        {
            free();
        }

        /**
         * free all the memory and resize to zero
         */
        void free()
        {
            assert((_reference_count==NULL) == (_allocated_memory==NULL)); //both null or both not null
            if (_allocated_memory != NULL){
                *_reference_count -= 1;
                if (*_reference_count < 0)
                    throw std::logic_error("reference_count cannot be less than zero.  Something went wrong");
                if (*_reference_count == 0){
                    delete[] _allocated_memory; 
                    delete _reference_count; 
                    _data = NULL; 
                    _allocated_memory = NULL; 
                    _reference_count = NULL;
                    _size = 0; 
                }
            }

        }

        /*
         * wrap another array
         */
        void wrap(Array<dtype> x)
        {
            free();
            _size = x._size;
            _reference_count = x._reference_count;
            _data = x._data;
            _allocated_memory = x._allocated_memory;
            if (_reference_count != NULL){
                *_reference_count += 1;
            }
        }
        
        /**
         * return a copy of the array.
         */
        Array<dtype> const copy()
        {
            Array<dtype> newarray(_size);
            for (size_t i=0; i<_size; ++i){
                newarray[i] = _data[i];
            }
            return newarray;
        }

        /**
         * return a view of the array.
         */
        Array<dtype> view(size_t ibegin, size_t iend)
        {
            if (iend <= ibegin) {
                throw std::invalid_argument("iend must larger than ibegin");
            }
            Array<dtype> newarray();
            newarray._data = &_data[ibegin];
            newarray._allocated_memory = _allocated_memory;
            newarray._size = iend - ibegin;
            newarray._reference_count = _reference_count;
            assert((_reference_count==NULL) == (_allocated_memory==NULL)); //both null or both not null
            if (_allocated_memory != NULL){
                *_reference_count += 1;
            }
            return newarray;
        }

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


        /**
         * Resize the array.  Only allowed if we are the sole owner of this
         * data or if the data has not been allocated yet.
         */
        void resize(size_t size) {
            // do sanity checks
            assert((_reference_count==NULL) == (_allocated_memory==NULL)); //both null or both not null
            if (_allocated_memory == NULL && _data != NULL){
                // this instance can occur if you wrap data that was not allocated by the Array class.
                // e.g. if this wraps data in a std::vector.
                throw std::runtime_error("Array: cannot resize Arrays if not sole owner of data");
            }
            if (_reference_count != NULL){
                if (*_reference_count != 1){
                    throw std::runtime_error("Array: cannot resize Arrays if not sole owner of data");
                }
            }
            // all seems ok. resize the array
            if (_allocated_memory != NULL) delete[] _allocated_memory;
            if (_reference_count != NULL) delete _reference_count;
            _size = size;
            _allocated_memory = new dtype[_size];
            _data = _allocated_memory;
            _reference_count = new long int;
            *_reference_count = 1;
        }

        /// access an element in the array
        dtype &operator[](size_t i) { return _data[i]; }
        dtype operator[](size_t i) const { return _data[i]; }

        /// read only access to element in array
//        dtype operator()(size_t i) const { return _data[i]; }

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
            out << a[i];
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

#ifndef PYGMIN_ARRAY_H
#define PYGMIN_ARRAY_H

#include <assert.h>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <math.h>
#include <memory>

namespace pele {

template<typename dtype>
class _ArrayMemory{
public:
    std::vector<dtype> _vector;
    dtype *_data;
    size_t _size;
    bool _owned;

    _ArrayMemory()
        : _vector(), _data(_vector.data()), _size(_vector.size()), _owned(true)
    {}

    _ArrayMemory(size_t size)
        : _vector(size), _data(_vector.data()), _size(_vector.size()), _owned(true)
    {}

    _ArrayMemory(size_t size, dtype const & val)
        : _vector(size, val), _data(_vector.data()), _size(_vector.size()), _owned(true)
    {}

    /**
     * wrap some data that is passed.  Do not take ownership of the data.
     */
    _ArrayMemory(dtype * data, size_t size)
        : _vector(), _data(data), _size(size), _owned(false)
    {}


    bool owned() const { return _owned; }

//    /**
//     * free all the memory and resize to zero
//     */
//    void free()
//    {
//        _vector.clear();
//        _data = _vector.data();
//        _size = 0;
//        _owned = true;
//    }



    void resize(size_t size) {
//        if (!_owned){
//            raise std::runtime_error("can't resize an array that is not owned")
//        }
        _vector.resize(size);
        _data = _vector.data();
        _size = size;
        _owned = true;
    }

    size_t size() const {
//        if (_owned)
//            return _vector.size();
//        else
        return _size;
    }

    /// return pointer to data
    dtype *data() { return _data; }
    dtype const *data() const { return _data; }



};

template<typename dtype>
class Array
{
    std::shared_ptr<_ArrayMemory<dtype> > _vector;
public:
    Array()
        : _vector(new _ArrayMemory<dtype>())
    {}
    Array(size_t size)
        : _vector(new _ArrayMemory<dtype>(size))
    {}
    Array(size_t size, dtype const & val)
    : _vector(new _ArrayMemory<dtype>(size, val))
    {}

    /**
     * wrap some data that is passed.  Do not take ownership of the data.
     */
    Array(dtype *data, size_t size)
    : _vector(new _ArrayMemory<dtype>(data, size))
    {}

    /**
     * wrap a vector.  This memory should never be deleted.
     */
    Array(std::vector<dtype> &x)
    : _vector(new _ArrayMemory<dtype>(x.data(), x.size()))
    {}

    /**
     * wrap another array is implemented by shared_ptr parent class
     *
     * Note, the input here is a const array, but this constructs a modifiable array.  This
     * is a loophole around the const declaration, but i'm not sure there is any way around it.
     * The compiler complains if this is not here.  e.g. for constructions like Array(x.copy())
     */
     Array(Array<dtype> const & x)
         : _vector(x._vector)
     {}

    /*
     * wrap another array
     */
    void wrap(Array<dtype> x)
    {
        _vector = x._vector;
    }

    /**
     * return pointer to data
     */
    inline dtype *data() { return _vector->data(); }
    inline dtype const *data() const { return _vector->data(); }

    /**
     * return the size of the array
     */
    inline size_t size() const { return _vector->size(); }

    /**
     * access an element in the array
     */
    inline dtype &operator[](const size_t i) { return data()[i]; }
    inline dtype const & operator[](const size_t i) const { return data()[i]; }

    /**
     * return iterators over data
     */
    typedef dtype * iterator;
    typedef dtype const * const_iterator;
    iterator begin() { return data(); }
    iterator end() { return data() + size(); }
    const_iterator begin() const { return data(); }
    const_iterator end() const { return data() + size(); }

    /**
     * test if they wrap the same data
     */
    inline bool operator==(Array<dtype> const rhs) const {
        return _vector.get() == rhs._vector.get();
    }
    inline bool operator!=(Array<dtype> const rhs) const {
        return !operator==(rhs);
    }


    /**
     * Assignment function: copy the data into the existing array
     *
     * arrays must be of same size
     */
    Array<dtype> &assign(const Array<dtype> & rhs) {
        if ((*this) != rhs) //check for self assignment
        {
            if (size() != rhs.size()){
                throw std::runtime_error("arrays must have the same size during assignment");
            }
            for (size_t i=0; i<size(); ++i)
                data()[i] = rhs[i];
        }
        return *this;
    }

    Array<dtype> &assign(dtype d) {
        for(size_t i=0; i<size(); ++i)
            data()[i] = d;
        return *this;
    }

    /**
     * return a copy of the array.
     */
    Array<dtype> copy() const
    {
        Array<dtype> newarray(size());
        newarray.assign(*this);
        return newarray;
    }

    /**
     * resize the array
     */
    void resize(size_t size){
        _vector->resize(size);
    }

    /**
     * free all the memory and set the size to zero
     */
    void free(){
        _vector = std::make_shared<_ArrayMemory<dtype> >();
    }

    /**
     * Returns whether the array is empty (i.e. whether its size is 0).
     */
    inline bool empty() const
    {
        return size() == 0;
    }

    long int reference_count() const
    {
        return _vector.use_count();
    }

    /*
     * Compound Assignment Operators += -= *=
     */
    Array<dtype> &operator+=(const Array<dtype> & rhs){
        if (size() != rhs.size()){
            throw std::runtime_error("operator+=: arrays must have the same size");
        }
        for (size_t i=0; i<size(); ++i)
            data()[i] += rhs[i];
        return *this;
    }

    Array<dtype> &operator+=(const dtype &rhs) {
        for (size_t i=0; i<size(); ++i)
            data()[i] += rhs;
        return *this;
   }

    Array<dtype> &operator-=(const Array<dtype> & rhs){
        if (size() != rhs.size()){
            throw std::runtime_error("operator-=: arrays must have the same size");
        }
        for (size_t i=0; i<size(); ++i)
            data()[i] -= rhs[i];
        return *this;
    }

    Array<dtype> &operator-=(const dtype &rhs) {
        for (size_t i=0; i<size(); ++i)
            data()[i] -= rhs;
        return *this;
   }

    Array<dtype> &operator*=(const Array<dtype> & rhs){
        if (size() != rhs.size()){
            throw std::runtime_error("operator*=: arrays must have the same size");
        }
        for (size_t i=0; i<size(); ++i)
            data()[i] *= rhs[i];
        return *this;
    }

    Array<dtype> &operator*=(const dtype &rhs) {
        for (size_t i=0; i<size(); ++i)
            data()[i] *= rhs;
        return *this;
    }


    Array<dtype> &operator/=(const Array<dtype> & rhs){
        if (size() != rhs.size()){
            throw std::runtime_error("operator/=: arrays must have the same size");
        }
        for (size_t i=0; i<size(); ++i)
            data()[i] /= rhs[i];
        return *this;
    }

    Array<dtype> &operator/=(const  dtype &rhs) {
        for (size_t i=0; i<size(); ++i)
            data()[i] /= rhs;
        return *this;
    }

/*SOME OTHER ARITHMETIC UTILITIES*/

    //returns the sum of all elements (reduces the array)
    const dtype sum() const {
        if (empty()){
            throw std::runtime_error("array::sum(): array is empty, can't sum array elements");
        }
        dtype sum_array=0;
        for (size_t i=0; i<size(); ++i)
            sum_array += data()[i] ;
        return sum_array;
    }

    //returns the product of all elements (reduces the array)
    const dtype prod() const {
        if (empty())
            throw std::runtime_error("array::prod(): array is empty, can't take product of array elements");

        dtype prod_array=1;
        for (size_t i=0; i<size(); ++i)
            prod_array *= data()[i] ;
        return prod_array;
    }


};


//    /**
//     * Simple wrapper class for arrays
//     *
//     * The Array class should provide a simple array handling interface
//     * and allow to efficiently transfer data between C++/Fortran/Python
//     */
//    template<typename dtype>
//    class Array {
//        dtype *_data;
//        dtype *_allocated_memory;
//        size_t _size;
//
//        long int *_reference_count;
//    public:
//        Array()
//            : _data(NULL), _allocated_memory(NULL), _size(0), _reference_count(NULL)
//        {}
//
//        /**
//         * create array with specific size and allocate memory
//         */
//        Array(size_t size)
//            : _data(NULL), _allocated_memory(NULL), _size(size), _reference_count(NULL)
//        {
//            if (size > 0){
//                _allocated_memory = new dtype[size];
//                _data = _allocated_memory;
//                _reference_count = new long int;
//                *_reference_count = 1;
//            }
//        }
//
//        /**
//         * create array with specific size and values and allocate memory
//         */
//
//        Array(size_t size, dtype val)
//            : _data(NULL), _allocated_memory(NULL), _size(size), _reference_count(NULL)
//        {
//            if (size > 0){
//                _allocated_memory = new dtype[size];
//                _data = _allocated_memory;
//                _reference_count = new long int;
//                *_reference_count = 1;
//                for (size_t i=0; i<_size; ++i)
//                {
//                    _data[i] = val;
//                }
//            }
//        }
//
//        /**
//         * wrap another array
//         */
//        Array(Array<dtype> const &x) :
//            _data(x._data), _allocated_memory(x._allocated_memory),
//            _size(x._size), _reference_count(x._reference_count)
//        {
//            //std::cout << "copy constructing Array\n";
//            if (_data == NULL){
//                throw std::runtime_error("cannot wrap an array with no data");
//            }
//            assert((_reference_count==NULL) == (_allocated_memory==NULL)); //both null or both not null
//            if (_reference_count != NULL){
//                *_reference_count += 1;
//            }
////            std::cout << "copy constructor: reference count " << _reference_count
////                    << " " << *_reference_count << "\n";
//        }
//
//        /**
//         * wrap some data that is passed.  Do not take ownership of the data.
//         */
//        Array(dtype *data, size_t size)
//            : _data(data), _allocated_memory(NULL), _size(size), _reference_count(NULL) {}
//
//        /**
//         * wrap a vector.  This memory should never be deleted.
//         */
//        Array(std::vector<dtype> &x) : _data(x.data()), _allocated_memory(NULL), _size(x.size()), _reference_count(NULL) { }
//
//        /**
//         * destructor
//         */
//        ~Array()
//        {
//            free();
//        }
//
//        /*Returns whether the array is empty (i.e. whether its size is 0).*/
//
//        bool empty() const
//        {
//            if (_data == NULL)
//                return true;
//            else
//                return false;
//        }
//
//        long int reference_count() const
//        {
//            if (_reference_count == NULL){
//                return 0;
//            } else {
//                return *_reference_count;
//            }
//        }
//
//        /**
//         * free all the memory and resize to zero
//         */
//        void free()
//        {
//            assert((_reference_count==NULL) == (_allocated_memory==NULL)); //both null or both not null
//            if (_allocated_memory != NULL){
//                *_reference_count -= 1;
//                if (*_reference_count < 0)
//                    throw std::logic_error("reference_count cannot be less than zero. Something went wrong");
//                if (*_reference_count == 0){
//                    delete[] _allocated_memory;
//                    delete _reference_count;
//                }
//            }
//            _allocated_memory = NULL;
//            _reference_count = NULL;
//            _data = NULL;
//            _size = 0;
//        }
//
//        /*
//         * wrap another array
//         */
//        void wrap(Array<dtype> x)
//        {
//            if (x._data == NULL){
//                throw std::runtime_error("cannot wrap an array with no data");
//            }
//            free();
//            _size = x._size;
//            _reference_count = x._reference_count;
//            _data = x._data;
//            _allocated_memory = x._allocated_memory;
//            if (_reference_count != NULL){
//                *_reference_count += 1;
//            }
//        }
//
//        /**
//         * return a copy of the array.
//         */
//        Array<dtype> const copy() const
//        {
//            Array<dtype> newarray(_size);
//            for (size_t i=0; i<_size; ++i){
//                newarray[i] = _data[i];
//            }
//            return newarray;
//        }
//
//        /**
//         * return a view of the array.
//         */
//        Array<dtype> view(size_t ibegin, size_t iend)
//        {
//            if (iend <= ibegin) {
//                throw std::invalid_argument("iend must larger than ibegin");
//            }
//            Array<dtype> newarray(*this);
//            newarray._data = &_data[ibegin];
//            newarray._size = iend - ibegin;
//            //newarray._allocated_memory = _allocated_memory;
//            //newarray._reference_count = _reference_count;
//            //assert((_reference_count==NULL) == (_allocated_memory==NULL)); //both null or both not null
//            //if (_allocated_memory != NULL){
//                //*_reference_count += 1;
//            //}
//            return newarray;
//        }
//
//        /**
//         * return pointer to data
//         */
//        dtype *data() { return _data; }
//        dtype const *data() const { return _data; }
//
//        /// return size of array
//        size_t size() const { return _size; }
//
//        /**
//         * Return true if this is the sole owner of the data
//         */
//        bool sole_owner() const
//        {
//            assert((_reference_count==NULL) == (_allocated_memory==NULL)); //both null or both not null
//            if (_reference_count == NULL){
//                if (_data != NULL){
//                    // the array wraps data which it doesn't own.
//                    return false;
//                }
//                assert(_size == 0);
//                // this array has no data
//                return true;
//            } else {
//                return (*_reference_count == 1);
//            }
//        }
//
//        // return iterators over data
//        typedef dtype * iterator;
//        typedef dtype const * const_iterator;
//        iterator begin() { return &_data[0]; }
//        iterator end() { return _data + _size; }
//        const_iterator begin() const { return _data; }
//        const_iterator end() const { return _data + _size; }
//
//
//        /**
//         * Resize the array.  Only allowed if we are the sole owner of this
//         * data or if the data has not been allocated yet.
//         */
//        void resize(size_t size) {
//            if (size == 0){
//                free();
//                return;
//            }
//            // do sanity checks
//            assert((_reference_count==NULL) == (_allocated_memory==NULL)); //both null or both not null
//            if (! sole_owner()){
//                // this instance can occur if you wrap data that was not allocated by the Array class.
//                // e.g. if this wraps data in a std::vector.
//                throw std::runtime_error("Array: cannot resize Arrays if not sole owner of data");
//            }
//            if (_size == size){
//                // no need to do anything.
//                // note, this is dangerous, as no new memory will be created
//                return;
//            }
//            // all seems ok. resize the array
//            free();
//            _size = size;
//            _allocated_memory = new dtype[_size];
//            _data = _allocated_memory;
//            _reference_count = new long int;
//            *_reference_count = 1;
//        }
//
//        /// access an element in the array
//        inline dtype &operator[](const size_t i) { return _data[i]; }
//        inline dtype const & operator[](const size_t i) const { return _data[i]; }
//
//        /**
//         * Assignment function: copy the data into the existing array
//         *
//         * arrays must be of same size
//         */
//        Array<dtype> &assign(const Array<dtype> & rhs) {
//            if (this != &rhs) //check for self assignement
//            {
//                if (_size != rhs.size()){
//                    throw std::runtime_error("arrays must have the same size during assignment");
//                }
//                for (size_t i=0; i<_size; ++i)
//                    _data[i] = rhs[i];
//            }
//            return *this;
//        }
//
//        Array<dtype> &assign(dtype d) {
//            for(size_t i=0; i<_size; ++i)
//                _data[i] = d;
//            return *this;
//        }
//
//
//        /*
//         * Assignment operator: WRAP the data
//         */
//        Array<dtype> &operator=(const Array<dtype> & rhs){
//            //if (_data != NULL) {
//                //std::cout << "operator=: cannot assign an array unless the array is unallocated\n";
//                //throw std::runtime_error("cannot assign an array unless the array is unallocated");
//            //}
//            if (this != &rhs) //check for self assignement
//            {
//                wrap(rhs);
//            }
//            return *this;
//        }
//
//        /*
//         * Compound Assignment Operators += -= *=
//         */
//        Array<dtype> &operator+=(const Array<dtype> & rhs){
//            if (_size != rhs.size()){
//                throw std::runtime_error("operator+=: arrays must have the same size");
//            }
//            for (size_t i=0; i<_size; ++i)
//                _data[i] += rhs[i];
//            return *this;
//        }
//
//        Array<dtype> &operator+=(const dtype &rhs) {
//            for (size_t i=0; i<_size; ++i)
//                _data[i] += rhs;
//            return *this;
//       }
//
//        Array<dtype> &operator-=(const Array<dtype> & rhs){
//            if (_size != rhs.size()){
//                throw std::runtime_error("operator-=: arrays must have the same size");
//            }
//            for (size_t i=0; i<_size; ++i)
//                _data[i] -= rhs[i];
//            return *this;
//        }
//
//        Array<dtype> &operator-=(const dtype &rhs) {
//            for (size_t i=0; i<_size; ++i)
//                _data[i] -= rhs;
//            return *this;
//       }
//
//        Array<dtype> &operator*=(const Array<dtype> & rhs){
//            if (_size != rhs.size()){
//                throw std::runtime_error("operator*=: arrays must have the same size");
//            }
//            for (size_t i=0; i<_size; ++i)
//                _data[i] *= rhs[i];
//            return *this;
//        }
//
//        Array<dtype> &operator*=(const dtype &rhs) {
//            for (size_t i=0; i<_size; ++i)
//                _data[i] *= rhs;
//            return *this;
//        }
//
//
//        Array<dtype> &operator/=(const Array<dtype> & rhs){
//            if (_size != rhs.size()){
//                throw std::runtime_error("operator/=: arrays must have the same size");
//            }
//            for (size_t i=0; i<_size; ++i)
//                _data[i] /= rhs[i];
//            return *this;
//        }
//
//        Array<dtype> &operator/=(const  dtype &rhs) {
//            for (size_t i=0; i<_size; ++i)
//                _data[i] /= rhs;
//            return *this;
//        }
//
//    /*SOME OTHER ARITHMETIC UTILITIES*/
//
//        //returns the sum of all elements (reduces the array)
//        const dtype sum() const {
//            if (_data == NULL)
//                throw std::runtime_error("array::sum(): array is empty, can't sum array elements");
//
//            dtype sum_array=0;
//            for (size_t i=0; i<_size; ++i)
//                sum_array += _data[i] ;
//            return sum_array;
//        }
//
//        //returns the product of all elements (reduces the array)
//        const dtype prod() const {
//            if (_data == NULL)
//                throw std::runtime_error("array::prod(): array is empty, can't take product of array elements");
//
//            dtype prod_array=1;
//            for (size_t i=0; i<_size; ++i)
//                prod_array *= _data[i] ;
//            return prod_array;
//        }
//    };


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

    /**
     * compute the dot product of two Arrays
     */
    inline double dot(Array<double> const &v1, Array<double> const &v2)
    {
      assert(v1.size() == v2.size());
      double dot = 0.;
      for (size_t i=0; i<v1.size(); ++i) {
        dot += v1[i] * v2[i];
      }
      return dot;
    }

    /**
     * compute the L2 norm of an Array
     */
    inline double norm(Array<double> const &v)
    {
      return sqrt(dot(v, v));
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

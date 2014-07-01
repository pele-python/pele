#ifndef PYGMIN_ARRAY_H
#define PYGMIN_ARRAY_H

#include <assert.h>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <math.h>
#include <memory>
#include <algorithm>

namespace pele {

template<typename dtype>
class _ArrayMemory{
    std::vector<dtype> _vector;
    dtype *_data;
    size_t _size;

public:
    _ArrayMemory()
        : _vector(), _data(_vector.data()), _size(_vector.size())
    {}

    _ArrayMemory(size_t size)
        : _vector(size), _data(_vector.data()), _size(_vector.size())
    {}

    _ArrayMemory(size_t size, dtype const & val)
        : _vector(size, val), _data(_vector.data()), _size(_vector.size())
    {}

    /**
     * wrap some data that is passed.  Do not take ownership of the data.
     */
    _ArrayMemory(dtype * data, size_t size)
        : _vector(), _data(data), _size(size)
    {}

    /**
     * resize the vector
     * 
     * This function does not guarantee that the data is preserved.
     * If it was previously wrapping externally allocated memory, that is
     * definitely lost.
     */
    void resize(size_t size) {
        _vector.resize(size);
        _data = _vector.data();
        _size = size;
    }

    /*
     * return the size of the array
     */
    inline size_t size() const {
        return _size;
    }

    /** 
     * return pointer to data
     */
    inline dtype *data() { return _data; }
    inline dtype const *data() const { return _data; }
};

template<typename dtype>
class Array
{
    std::shared_ptr<_ArrayMemory<dtype> > _memory;
public:
    /**
     * default constructor
     *
     * create an array of size 0
     */
    Array()
        : _memory(new _ArrayMemory<dtype>())
    {}

    /**
     * construct an array with a given size
     */
    Array(size_t size)
        : _memory(new _ArrayMemory<dtype>(size))
    {}

    /**
     * construct an array with a given size, each element is a copy of val
     */
    Array(size_t size, dtype const & val)
    : _memory(new _ArrayMemory<dtype>(size, val))
    {}

    /**
     * wrap some data that is passed.  Do not take ownership of the data.
     */
    Array(dtype *data, size_t size)
    : _memory(new _ArrayMemory<dtype>(data, size))
    {}

    /**
     * wrap a vector.  This memory should never be deleted.
     */
    Array(std::vector<dtype> &x)
    : _memory(new _ArrayMemory<dtype>(x.data(), x.size()))
    {}

    /**
     * wrap another array is implemented by shared_ptr parent class
     *
     * Note, the input here is a const array, but this constructs a modifiable array.  This
     * is a loophole around the const declaration, but i'm not sure there is any way around it.
     * The compiler complains if this is not here.  e.g. for constructions like Array(x.copy())
     *
     * This is commented because it just duplicates the default copy constructor
     *
     Array(Array<dtype> const & x)
         : _memory(x._memory)
     {}
     */

    /*
     * wrap another array
     */
    inline void wrap(Array<dtype> x)
    {
        _memory = x._memory;
    }

    /**
     * return pointer to data
     */
    inline dtype *data() { return _memory->data(); }
    inline dtype const *data() const { return _memory->data(); }

    /**
     * return the size of the array
     */
    inline size_t size() const { return _memory->size(); }

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
    inline iterator begin() { return data(); }
    inline iterator end() { return data() + size(); }
    inline const_iterator begin() const { return data(); }
    inline const_iterator end() const { return data() + size(); }

    /*
     * Assignment operator: WRAP the data
     *
     * This is commented because it just duplicates the default assignment operator
     *
    Array<dtype> &operator=(const Array<dtype> & rhs){
        _memory = rhs._memory;
    }
    */


    /**
     * equality operator
     *
     * test if they wrap the same data
     */
    inline bool operator==(Array<dtype> const rhs) const {
        return _memory.get() == rhs._memory.get();
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
            std::copy(rhs.begin(), rhs.end(), begin());
        }
        return *this;
    }

    /**
     * assign each element of the array to be d
     */
    Array<dtype> &assign(dtype const & d) {
        std::fill(begin(), end(), d);
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
    inline void resize(size_t size){
        _memory->resize(size);
    }

    /**
     * wrap a new empty array.
     */
    inline void free(){
        _memory = std::make_shared<_ArrayMemory<dtype> >();
    }

    /**
     * Returns whether the array is empty (i.e. whether its size is 0).
     */
    inline bool empty() const
    {
        return size() == 0;
    }

    inline long int reference_count() const
    {
        return _memory.use_count();
    }

    /*
     * Compound Assignment Operators += -= *=
     */
    Array<dtype> &operator+=(const Array<dtype> & rhs){
        if (size() != rhs.size()){
            throw std::runtime_error("operator+=: arrays must have the same size");
        }
        const_iterator iter = rhs.begin();
        for (dtype & val : *this){
            val += *iter;
            ++iter;
        }
        return *this;
    }

    Array<dtype> &operator+=(const dtype &rhs) {
        for (dtype & val : (*this)){
            val += rhs;
        }
        return *this;
    }

    Array<dtype> &operator-=(const Array<dtype> & rhs){
        if (size() != rhs.size()){
            throw std::runtime_error("operator-=: arrays must have the same size");
        }
        const_iterator iter = rhs.begin();
        for (dtype & val : *this){
            val -= *iter;
            ++iter;
        }
        return *this;
    }

    Array<dtype> &operator-=(const dtype &rhs) {
        for (dtype & val : (*this)){
            val -= rhs;
        }
        return *this;
   }

    Array<dtype> &operator*=(const Array<dtype> & rhs){
        if (size() != rhs.size()){
            throw std::runtime_error("operator*=: arrays must have the same size");
        }
        const_iterator iter = rhs.begin();
        for (dtype & val : *this){
            val *= *iter;
            ++iter;
        }
        return *this;
    }

    Array<dtype> &operator*=(const dtype &rhs) {
        for (dtype & val : (*this)){
            val *= rhs;
        }
        return *this;
    }


    Array<dtype> &operator/=(const Array<dtype> & rhs){
        if (size() != rhs.size()){
            throw std::runtime_error("operator/=: arrays must have the same size");
        }
        const_iterator iter = rhs.begin();
        for (dtype & val : *this){
            val /= *iter;
            ++iter;
        }
        return *this;
    }

    Array<dtype> &operator/=(const  dtype &rhs) {
        for (dtype & val : (*this)){
            val /= rhs;
        }
        return *this;
    }

/*SOME OTHER ARITHMETIC UTILITIES*/

    /**
     * returns the sum of all elements (reduces the array)
     */
    const dtype sum() const {
        if (empty()){
            throw std::runtime_error("array::sum(): array is empty, can't sum array elements");
        }
        dtype sum_array=0;
        for (dtype const & val : (*this)){
            sum_array += val;
        }
        return sum_array;
    }

    /**
     * returns the product of all elements (reduces the array)
     */
    const dtype prod() const {
        if (empty())
            throw std::runtime_error("array::prod(): array is empty, can't take product of array elements");
        dtype prod_array = 1;
        for (dtype const & val : (*this)){
            prod_array *= val;
        }
        return prod_array;
    }

    /**
     * return a view of the array.
     *
     * This is commented in case we ever want to add the functionality again.
     * This is an outline, it will need some modifications to actually work.
     *
    Array<dtype> view(size_t ibegin, size_t iend)
    {
        if (iend <= ibegin) {
            throw std::invalid_argument("iend must larger than ibegin");
        }
        if (iend > size()) {
            throw std::invalid_argument("iend cannot be larger than array size");
        }
        Array<dtype> newarray(*this);
        newarray._memory._data += ibegin;
        newarray._memory._size = iend - ibegin;
        return newarray;
    }
    */

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

#endif

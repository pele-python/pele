#ifndef PYGMIN_ARRAY_H
#define PYGMIN_ARRAY_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace pele {

/**
 * This manages the data of the Array class.  This can act as a simple wrapper for a
 * vector, or wrap an externally allocated block of memory.
 */
template<typename dtype>
class _ArrayMemory {
    std::vector<dtype> _vector;
    dtype *_data; /** _data will either point to the beginning of _vector, or to the beginning of the
                       block of externally allocated memory.  _data is a simple pointer and will not
                       represent allocated memory. */
    size_t _size;  /** the size of the block of memory, whether in vector or external. */

public:
    _ArrayMemory()
        : _vector(),
          _data(_vector.data()),
          _size(_vector.size())
    {}

    _ArrayMemory(size_t size)
        : _vector(size),
          _data(_vector.data()),
          _size(_vector.size())
    {}

    _ArrayMemory(size_t size, dtype const & val)
        : _vector(size, val),
          _data(_vector.data()),
          _size(_vector.size())
    {}

    /**
     * wrap some data that is passed.  Do not take ownership of the data.
     */
    _ArrayMemory(dtype * data, size_t size)
        : _vector(),
          _data(data),
          _size(size)
    {}

    /**
     * return the size of the array
     */
    inline size_t size() const { return _size; }

    /** 
     * return pointer to data
     */
    inline dtype *data() { return _data; }
    inline dtype const *data() const { return _data; }
};



/** An Array class which acts in many ways like a numpy array
 *
 * This copy constructor and assignment operator act to wrap existing
 * memory rather than copy the memory.  
 */
template<typename dtype>
class Array
{
protected:
    std::shared_ptr<_ArrayMemory<dtype> > _memory;
    dtype * _data; /**< _data will usually be a copy of memory->data().  If this
                        is a view of another array then _data will be
                        _memory->data() + ibegin */
    size_t _size;   /**< The size of the array. */
public:

    /** create an array of size 0
     */
    Array()
        : _memory(new _ArrayMemory<dtype>()),
          _data(_memory->data()),
          _size(_memory->size())
    {}

    /**
     * construct an array with a given size
     */
    Array(size_t size)
        : _memory(new _ArrayMemory<dtype>(size)),
          _data(_memory->data()),
          _size(_memory->size())
    {}

    /**
     * construct an array with a given size, each element is a copy of val
     */
    Array(size_t size, dtype const & val)
        : _memory(new _ArrayMemory<dtype>(size, val)),
          _data(_memory->data()),
          _size(_memory->size())

    {}

    /**
     * wrap some data that is passed.  Do not take ownership of the data.
     */
    Array(dtype *data, size_t size)
        : _memory(new _ArrayMemory<dtype>(data, size)),
          _data(_memory->data()),
          _size(_memory->size())
    {}

    /**
     * wrap the data in a vector.  This memory should never be deleted.
     */
    Array(std::vector<dtype> &x)
        : _memory(new _ArrayMemory<dtype>(x.data(), x.size())),
          _data(_memory->data()),
          _size(_memory->size())
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
           _data(x._data),
           _size(x._size)
     {}
     */

    /**
     * wrap another array
     */
    inline void wrap(Array<dtype> x)
    {
        _memory = x._memory;
        _data = x._data;
        _size = x._size;
    }

    /**
     * return pointer to data
     */
    inline dtype *data() { return _data; }
    inline dtype const *data() const { return _data; }

    /** return the size of the array
     * 
     */
    inline size_t size() const { return _size; }

    /**
     * access an element in the array
     */
    inline dtype &operator[](const size_t i) { return data()[i]; }
    inline dtype const &operator[](const size_t i) const { return data()[i]; }

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
        _data = rhs._data;
        _size = rhs._size;
    }
    */


    /**
     * return true if the two arrays wrap the same data
     */
    inline bool operator==(Array<dtype> const rhs) const
    {
        return data() == rhs.data() and size() == rhs.size();
    }
    inline bool operator!=(Array<dtype> const rhs) const
    {
        return !operator==(rhs);
    }


    /**
     * Assignment function: copy the data into the existing array
     *
     * arrays must be of same size
     */
    Array<dtype> &assign(const Array<dtype> & rhs)
    {
        //check for self assignment
        if ((*this) != rhs) {
            if (size() != rhs.size()) {
                throw std::runtime_error("arrays must have the same size during assignment");
            }
            std::copy(rhs.begin(), rhs.end(), begin());
        }
        return *this;
    }

    /**
     * assign each element of the array to be d
     */
    Array<dtype> &assign(dtype const & d)
    {
        std::fill(begin(), end(), d);
        return *this;
    }

//    /**
//     * assign each element of the array to be
//     */
//    Array<dtype> &assign(dtype const * const d) {
//        std::copy(rhs.begin(), rhs.end(), begin());
//        for (size_t i = 0; i < size(); ++i) {
//
//        }
//        return *this;
//    }

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
     * wrap a new empty array.
     */
    inline void free()
    {
        _memory = std::make_shared<_ArrayMemory<dtype> >();
        _data = _memory->data();
        _size = _memory->size();
    }

    /**
     * Returns whether the array is empty (whether its size is 0).
     */
    inline bool empty() const
    {
        return size() == 0;
    }

    inline long int reference_count() const
    {
        return _memory.use_count();
    }

    /**
     * Compound Assignment Operators += -= *=
     */
    Array<dtype> &operator+=(const Array<dtype> & rhs)
    {
        if (size() != rhs.size()) {
            throw std::runtime_error("operator+=: arrays must have the same size");
        }
        const_iterator iter = rhs.begin();
        for (dtype & val : *this) {
            val += *iter;
            ++iter;
        }
        return *this;
    }

    Array<dtype> &operator+=(const dtype &rhs)
    {
        for (dtype & val : (*this)) {
            val += rhs;
        }
        return *this;
    }

    Array<dtype> &operator-=(const Array<dtype> & rhs)
    {
        if (size() != rhs.size()) {
            throw std::runtime_error("operator-=: arrays must have the same size");
        }
        const_iterator iter = rhs.begin();
        for (dtype & val : *this) {
            val -= *iter;
            ++iter;
        }
        return *this;
    }

    Array<dtype> &operator-=(const dtype &rhs)
    {
        for (dtype & val : (*this)) {
            val -= rhs;
        }
        return *this;
    }

    Array<dtype> &operator*=(const Array<dtype> & rhs)
    {
        if (size() != rhs.size()) {
            throw std::runtime_error("operator*=: arrays must have the same size");
        }
        const_iterator iter = rhs.begin();
        for (dtype & val : *this) {
            val *= *iter;
            ++iter;
        }
        return *this;
    }

    Array<dtype> &operator*=(const dtype &rhs)
    {
        for (dtype & val : (*this)) {
            val *= rhs;
        }
        return *this;
    }
    
    Array<dtype> operator*(const dtype rhs)
    {
        Array<dtype> result = this->copy();
        return (result *= rhs).copy();
    }

    Array<dtype> &operator/=(const Array<dtype> & rhs)
    {
        if (size() != rhs.size()) {
            throw std::runtime_error("operator/=: arrays must have the same size");
        }
        const_iterator iter = rhs.begin();
        for (dtype & val : *this) {
            val /= *iter;
            ++iter;
        }
        return *this;
    }

    Array<dtype> &operator/=(const  dtype &rhs)
    {
        for (dtype & val : (*this)) {
            val /= rhs;
        }
        return *this;
    }

/*SOME OTHER ARITHMETIC UTILITIES*/

    /**
     * returns the sum of all elements (reduces the array)
     */
    dtype sum() const
    {
        if (empty()) {
            throw std::runtime_error("array::sum(): array is empty, can't sum array elements");
        }
        return std::accumulate(begin(), end(), dtype(0));
    }

    /**
     * returns the product of all elements (reduces the array)
     *
     * References:
     * http://www.cplusplus.com/reference/functional/multiplies/
     * http://en.cppreference.com/w/cpp/algorithm/accumulate
     * http://rosettacode.org/wiki/Sum_and_product_of_an_array
     */
    dtype prod() const
    {
        if (empty()) {
            throw std::runtime_error("array::prod(): array is empty, can't take product of array elements");
        }
        return std::accumulate(begin(), end(), dtype(1), std::multiplies<dtype>());
    }

    /**
     * return an array that wraps the data from index ibegin to index iend-1.
     */
    Array<dtype> view(size_t ibegin, size_t iend) const
    {
        if (iend <= ibegin) {
            throw std::invalid_argument("iend must larger than ibegin");
        }
        if (iend > size()) {
            throw std::invalid_argument("iend cannot be larger than array size");
        }
        Array<dtype> newarray(*this);
        newarray._data += ibegin;
        newarray._size = iend - ibegin;
        return newarray;
    }
    
    /**
     * Get maximum and minimum elements of array.
     */
    dtype get_max() const { return *std::max_element(begin(), end()); }
    dtype get_min() const { return *std::min_element(begin(), end()); }
};



// for array printing
template<class dtype>
inline std::ostream &operator<<(std::ostream &out, const Array<dtype> &a)
{
    out << "[ ";
    for(size_t i = 0; i < a.size(); ++i) {
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
  return std::inner_product(v1.begin(), v1.end(), v2.begin(), double(0));
}

/**
 * compute the L2 norm of an Array
 */
inline double norm(Array<double> const &v)
{
  return sqrt(dot(v, v));
}

template<class T, class U>
Array<T> operator*(const U rhs, const Array<T>& lhs)
{
    Array<T> result = lhs.copy();
    return (result *= rhs).copy();
}

} // namespace pele

#endif // #ifndef PYGMIN_ARRAY_H

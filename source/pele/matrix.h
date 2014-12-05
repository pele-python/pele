#ifndef _PELE_MATRIX_H_
#define _PELE_MATRIX_H_

namespace pele{

/**
 * This is a very minimal implementation of a matrix.  It's primary function is
 * to act as a wrapper for pele::Array, so a pele array can be act as a matrix
 * temporarily.  The idea is to redo somthing like the reshape() function in
 * numpy.
 */
template<class dtype>
class MatrixAdapter : public pele::Array<dtype> {
public:
    /**
     * the second dimension of the matrix, e.g. the number of colums
     *
     * note: if we make this const we can only use the assignment operator on
     * matrices with the same first dimension
     */
    size_t _dim2;

    /**
     * Construct a matrix and allocate memory for it.
     */
    MatrixAdapter(size_t dim1, size_t dim2, dtype val=0)
        : pele::Array<dtype>(dim1 * dim2, val),
          _dim2(dim2)
    {}

    /**
     * wrap a pele::Array to act like a matrix
     *
     * This is like numpy.reshape.  v.size() must be divisable by dim2
     */
    MatrixAdapter(pele::Array<double> v, size_t dim2)
        : pele::Array<dtype>(v),
          _dim2(dim2)
    {
        if (v.size() % dim2 != 0) {
            throw std::invalid_argument("v.size() is not divisible by dim2");
        }
    }

    /**
     * wrap an existing block of memory
     */
    MatrixAdapter(double * data, size_t dim1, size_t dim2)
        : pele::Array<dtype>(data, dim1*dim2),
          _dim2(dim2)
    {}

    /**
     * access the element of the array at row i and column j
     */
    inline dtype const & operator()(size_t i, size_t j) const
    {
        return this->operator[](i * _dim2 + j);
    }
    inline dtype & operator()(size_t i, size_t j)
    {
        return this->operator[](i * _dim2 + j);
    }

    /**
     * return a pair giving the shape of the array.
     */
    inline std::pair<size_t, size_t> shape() const
    {
        return std::pair<size_t, size_t>(this->size() / _dim2, _dim2);
    }
};

/**
 * multiply two matrices.  Note, this is a very inefficient way of doing matrix
 * multiplication.  If you have large matrices or care about speed you should
 * use something else.
 */
template<class dtype>
MatrixAdapter<dtype> hacky_mat_mul(MatrixAdapter<dtype> const & A, MatrixAdapter<dtype> const & B)
{
    assert(A.shape().second == B.shape().first);
    size_t const L = A.shape().second;
    size_t const N = A.shape().first;
    size_t const M = B.shape().second;

    MatrixAdapter<dtype> C(N, M, 0);
    for (size_t i = 0; i<N; ++i){
        for (size_t j = 0; j<M; ++j){
            double val = 0;
            for (size_t k = 0; k<L; ++k){
                val += A(i,k) * B(k,j);
            }
            C(i,j) = val;
        }
    }
    return C;
}

///**
// * multiply a matrix times an vector
// */
//template<class dtype>
//pele::Array<dtype> hacky_mat_mul(MatrixAdapter<dtype> const & A, pele::Array<dtype> const & v)
//{
//    assert(A.shape().second == v.size());
//    size_t const L = A.shape().second;
//    size_t const n = A.shape().first;
//
//    pele::Array<dtype> C(n, 0);
//    for (size_t i = 0; i<n; ++i){
//        dtype val = 0;
//        for (size_t k = 0; k<L; ++k){
//            val += A(i,k) * v[k];
//        }
//        C(i) = val;
//    }
//    return C;
//}

// for matrix printing
template<class dtype>
std::ostream &operator<<(std::ostream &out, const pele::MatrixAdapter<dtype> &a) {
    out << "[ ";
    size_t const N = a.shape().first;
    size_t const M = a.shape().second;
    for(size_t n=0; n<N;++n) {
        for(size_t m=0; m<M;++m) {
            if(m>0) out << ", ";
            out << a(n,m);
        }
        if (n < N-1) out << ",\n  ";
    }
    out << " ]";
    return out;
}


}
#endif

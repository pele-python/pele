#include "pele/array.h"
#include "pele/lowest_eig_potential.h"

#include <iostream>
#include <stdexcept>
#include <vector>
#include <gtest/gtest.h>
#include <cmath>
#include <memory>

using pele::Array;
using pele::BasePotential;

class OrthogonalizeTranslationalTest :  public ::testing::Test
{
public:

    pele::Array<double> ev0, ev1, ev2, _vector, _coords;
    size_t _natoms, _bdim, _ndim;
    std::shared_ptr<pele::Orthogonalize> _orthog;

    virtual void SetUp(){
        _natoms = 10;
        _bdim = 3;
        _ndim = _natoms*_bdim;
        _orthog = std::shared_ptr<pele::Orthogonalize> (new pele::OrthogonalizeTranslational(_natoms, _bdim));
        ev0.resize(_ndim);
        ev1.resize(_ndim);
        ev2.resize(_ndim);
        ev0.assign(0.0);
        ev1.assign(0.0);
        ev2.assign(0.0);
        double v = 1/sqrt(_natoms);
        for(size_t j=0;j<_ndim;j+=_bdim)
            ev0[j] = v;
        for(size_t j=1;j<_ndim;j+=_bdim)
            ev1[j] = v;
        for(size_t j=2;j<_ndim;j+=_bdim)
            ev2[j] = v;
        _vector.resize(_ndim);
        _coords.resize(_ndim);
        for(size_t j=0;j<_ndim;++j)
            _vector[j] = (double) j;
        _coords.assign(10);
};
};

TEST_F(OrthogonalizeTranslationalTest, orthogonalize_works0){
    _orthog->orthogonalize(_coords, _vector);
    double dot_prod = dot(_vector,ev0);
    EXPECT_NEAR(dot_prod, 0, 1e-10);
}

TEST_F(OrthogonalizeTranslationalTest, orthogonalize_works1){
    _orthog->orthogonalize(_coords, _vector);
    double dot_prod = dot(_vector,ev1);
    std::cout<<_vector<<std::endl;
    EXPECT_NEAR(dot_prod, 0, 1e-10);
}

TEST_F(OrthogonalizeTranslationalTest, orthogonalize_works2){
    _orthog->orthogonalize(_coords, _vector);
    double dot_prod = dot(_vector,ev2);
    EXPECT_NEAR(dot_prod, 0, 1e-10);
}

TEST_F(OrthogonalizeTranslationalTest, orthogonalize_works3){
    _orthog->orthogonalize(_coords, _vector);
    for(size_t i=0;i<_ndim;++i)
        ASSERT_EQ(_coords[i],10);
}



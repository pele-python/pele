#include "pele/array.h"
#include "pele/lowest_eig_potential.h"
#include "pele/lbfgs.h"
#include "pele/lj.h"

#include <iostream>
#include <stdexcept>
#include <vector>
#include <gtest/gtest.h>
#include <cmath>
#include <memory>

using pele::Array;
using pele::BasePotential;

class LowestEigPotentialTest :  public ::testing::Test
{
    public:
    std::shared_ptr<pele::BasePotential> _potential;
    std::shared_ptr<pele::BasePotential> _lowesteigpot;
    std::shared_ptr<pele::GradientOptimizer> _lbfgs;
    size_t _natoms, _bdim, _ndim, _M;
    pele::Array<double> _x, _g, _ranvec;
    pele::Array<double> ev0, ev1, ev2;
    double _c6,_c12, _tol;

    virtual void SetUp(){
        _natoms = 4;
        _bdim=3;
        _ndim=_bdim*_natoms;
        _c6 = 1.2;
        _c12 = 2.3;
        _x = Array<double>(_bdim*_natoms);
        _x[0] = 0.1;
        _x[1] = 0.2;
        _x[2] = 0.3;
        _x[3] = 0.44;
        _x[4] = 0.55;
        _x[5] = 1.66;
        _x[6] = 0;
        _x[7] = 0;
        _x[8] = -0.5;
        _x[9] = 0;
        _x[10] = +0.5;
        _x[11] = 0.;
        ev0 = Array<double>(_ndim);
        ev1 = Array<double>(_ndim);
        ev2 = Array<double>(_ndim);
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
        //setup potential
        _potential = std::shared_ptr<pele::BasePotential> (new pele::LJ(_c6, _c12));
        //setup LBFGS
        _tol = 1e-5;
        _M = 4;
        _lbfgs = std::make_shared<pele::LBFGS>(_potential, _x, _tol, _M);
        //setup lowesteigtest
        _g = Array<double>(_ndim);
        _ranvec = Array<double>(_ndim);
        for(size_t j=0;j<_ndim;++j)
            _ranvec[j] = (double) j;
        _ranvec /= norm(_ranvec);
    }

};

TEST_F(LowestEigPotentialTest, lbfgs_quench_works){
    _lbfgs->run(100);
    bool success = _lbfgs->success();
    ASSERT_TRUE(success);
}

TEST_F(LowestEigPotentialTest, lowesteigtest_works){
    _lbfgs->run(100);
    bool success = _lbfgs->success();
    ASSERT_TRUE(success);
    _x = _lbfgs->get_x();
    _g = _lbfgs->get_g();
    _lowesteigpot = std::shared_ptr<pele::BasePotential> (new pele::LowestEigPotential(_potential, _x, _bdim));
    std::shared_ptr<pele::LBFGS> lbfgs = std::shared_ptr<pele::LBFGS> (new pele::LBFGS(_lowesteigpot, _ranvec, _tol, _M));
    lbfgs->run(100);
    ASSERT_TRUE(success);
}

TEST_F(LowestEigPotentialTest, lowesteigtest_works2){
    _lbfgs->run(100);
    bool success = _lbfgs->success();
    ASSERT_TRUE(success);
    _x = _lbfgs->get_x();
    _g = _lbfgs->get_g();
    _lowesteigpot = std::shared_ptr<pele::BasePotential> (new pele::LowestEigPotential(_potential, _x, _bdim));
    std::shared_ptr<pele::LBFGS> lbfgs = std::shared_ptr<pele::LBFGS> (new pele::LBFGS(_lowesteigpot, _ranvec, _tol, _M));
    lbfgs->run(100);
    ASSERT_TRUE(success);
    double lowesteigenvalue = lbfgs->get_f();
    ASSERT_NEAR(lowesteigenvalue, 0, 1e-5);
}

TEST_F(LowestEigPotentialTest, LowestEig_ranvecorth){
    _lbfgs->run(100);
    bool success = _lbfgs->success();
    ASSERT_TRUE(success);
    _x = _lbfgs->get_x();
    _g = _lbfgs->get_g();
    _lowesteigpot = std::shared_ptr<pele::BasePotential> (new pele::LowestEigPotential(_potential, _x, _bdim));
    std::shared_ptr<pele::LBFGS> lbfgs = std::shared_ptr<pele::LBFGS> (new pele::LBFGS(_lowesteigpot, _ranvec, _tol, _M));
    lbfgs->run(100);
    ASSERT_TRUE(success);
    double lowesteigenvalue = lbfgs->get_f();
    ASSERT_NEAR(lowesteigenvalue, 0, 1e-5);
    _ranvec = lbfgs->get_x();
    double dot_prod = dot(_ranvec,ev0);
    EXPECT_NEAR(dot_prod, 0, 1e-10);
    dot_prod = dot(_ranvec,ev1);
    EXPECT_NEAR(dot_prod, 0, 1e-10);
    dot_prod = dot(_ranvec,ev2);
    EXPECT_NEAR(dot_prod, 0, 1e-10);
}

class OrthogonalizeTranslationalTest :  public ::testing::Test
{
public:

    pele::Array<double> ev0, ev1, ev2, _vector, _coords;
    size_t _natoms, _bdim, _ndim;
    std::shared_ptr<pele::Orthogonalize> _orthog;

    virtual void SetUp(){
        _natoms = 4;
        _bdim = 3;
        _ndim = _natoms*_bdim;
        _orthog = std::shared_ptr<pele::Orthogonalize> (new pele::OrthogonalizeTranslational(_natoms, _bdim));
        ev0 = Array<double>(_ndim);
        ev1 = Array<double>(_ndim);
        ev2 = Array<double>(_ndim);
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
        _vector = Array<double>(_ndim);
        _coords = Array<double>(_ndim);
        for(size_t j=0;j<_ndim;++j)
            _vector[j] = (double) j;
        _vector /= norm(_vector);
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



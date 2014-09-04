#include <gtest/gtest.h>
#include <limits>

#include "pele/aatopology.h"
#include "pele/lj.h"
#include "pele/vecn.h"

using pele::Array;
using pele::VecN;
using pele::MatrixNM;
using pele::CoordsAdaptor;
using pele::norm;

static double const EPS = std::numeric_limits<double>::min();
#define EXPECT_NEAR_RELATIVE(A, B, T)  ASSERT_NEAR(A/(fabs(A)+fabs(B) + EPS), B/(fabs(A)+fabs(B) + EPS), T)


TEST(AA2RotMat, Works)
{
    VecN<3> p;
    for (size_t i = 0; i < p.size(); ++i) p[i] = i+1;
    p /= norm<3>(p);
    auto mx = pele::aa_to_rot_mat(p);
    ASSERT_NEAR(mx(0,0), 0.57313786, 1e-5);
    ASSERT_NEAR(mx(1,0), 0.74034884, 1e-5);
    ASSERT_NEAR(mx(2,0), -0.35127851, 1e-4);
    ASSERT_NEAR(mx(0,1), -0.60900664, 1e-5);
    ASSERT_NEAR(mx(1,1), 0.6716445, 1e-5);
    ASSERT_NEAR(mx(2,1), 0.42190588, 1e-5);
    ASSERT_NEAR(mx(0,2), 0.54829181, 1e-5);
    ASSERT_NEAR(mx(1,2), -0.02787928, 1e-5);
    ASSERT_NEAR(mx(2,2), 0.83582225, 1e-5);
}

TEST(RotMatDerivs, Works)
{
    VecN<3> p;
    for (size_t i = 0; i < p.size(); ++i) p[i] = i+1;
    p /= norm<3>(p);
    MatrixNM<3,3> mx(1.);
    MatrixNM<3,3> drm1(0.);
    MatrixNM<3,3> drm2(0.);
    MatrixNM<3,3> drm3(0.);
    pele::rot_mat_derivatives(p, mx, drm1, drm2, drm3);

    // mx
    ASSERT_NEAR(mx(0,0), 0.57313786, 1e-5);
    ASSERT_NEAR(mx(1,0), 0.74034884, 1e-5);
    ASSERT_NEAR(mx(2,0), -0.35127851, 1e-4);
    ASSERT_NEAR(mx(0,1), -0.60900664, 1e-5);
    ASSERT_NEAR(mx(1,1), 0.6716445, 1e-5);
    ASSERT_NEAR(mx(2,1), 0.42190588, 1e-5);
    ASSERT_NEAR(mx(0,2), 0.54829181, 1e-5);
    ASSERT_NEAR(mx(1,2), -0.02787928, 1e-5);
    ASSERT_NEAR(mx(2,2), 0.83582225, 1e-5);

    // drm1
    ASSERT_NEAR(drm1(0,0), 0.01933859, 1e-5);
    ASSERT_NEAR(drm1(0,2), 0.32109128, 1e-5);
    ASSERT_NEAR(drm1(1,1), -0.23084292, 1e-5);
    ASSERT_NEAR(drm1(2,0), 0.40713948, 1e-5);
    ASSERT_NEAR(drm1(2,2), -0.23828083, 1e-5);

    // drm2
    ASSERT_NEAR(drm2(0,0), -0.45276033, 1e-5);
    ASSERT_NEAR(drm2(0,2), 0.74649729, 1e-5);
    ASSERT_NEAR(drm2(1,1), 0.02975168, 1e-5);
    ASSERT_NEAR(drm2(2,0), -0.76434829, 1e-5);
    ASSERT_NEAR(drm2(2,2), -0.47656167, 1e-5);

    // drm3
    ASSERT_NEAR(drm3(0,0), -0.67914049, 1e-5);
    ASSERT_NEAR(drm3(0,2), -0.01960117, 1e-5);
    ASSERT_NEAR(drm3(1,1), -0.69252875, 1e-5);
    ASSERT_NEAR(drm3(2,0), 0.23854341, 1e-5);
    ASSERT_NEAR(drm3(2,2), 0.02231376, 1e-5);
}

TEST(RotMatDerivsSmallTheta, Works)
{
    VecN<3> p;
    for (size_t i = 0; i < p.size(); ++i) p[i] = i+1;
    p /= norm<3>(p) * 1e7;
    MatrixNM<3,3> mx(1.);
    MatrixNM<3,3> drm1(0.);
    MatrixNM<3,3> drm2(0.);
    MatrixNM<3,3> drm3(0.);
    pele::rot_mat_derivatives_small_theta(p, mx, drm1, drm2, drm3, true);

    // mx
    EXPECT_NEAR_RELATIVE(mx(0,0), 1., 1e-5);
    EXPECT_NEAR_RELATIVE(mx(1,0), 8.01783726e-08, 1e-5);
    EXPECT_NEAR_RELATIVE(mx(2,0), -5.34522484e-08, 1e-4);
    EXPECT_NEAR_RELATIVE(mx(0,1), -8.01783726e-08, 1e-5);
    EXPECT_NEAR_RELATIVE(mx(1,1), 1, 1e-5);
    EXPECT_NEAR_RELATIVE(mx(2,1), 2.67261242e-08, 1e-5);
    EXPECT_NEAR_RELATIVE(mx(0,2), 5.34522484e-08, 1e-5);
    EXPECT_NEAR_RELATIVE(mx(1,2), -2.67261242e-08, 1e-5);
    EXPECT_NEAR_RELATIVE(mx(2,2), 1, 1e-5);

    // drm1
    EXPECT_NEAR_RELATIVE(drm1(0,0), 0, 1e-5);
    EXPECT_NEAR_RELATIVE(drm1(0,2), 4.00891863e-08, 1e-5);
    EXPECT_NEAR_RELATIVE(drm1(1,1), -2.67261242e-08, 1e-5);
    EXPECT_NEAR_RELATIVE(drm1(2,0), 4.00891863e-08, 1e-5);
    EXPECT_NEAR_RELATIVE(drm1(2,2), -2.67261242e-08, 1e-5);

    // drm2
    EXPECT_NEAR_RELATIVE(drm2(0,0), -5.34522484e-08, 1e-5);
    ASSERT_NEAR(drm2(0,2), 1, 1e-5);
    ASSERT_NEAR(drm2(1,1), 0., 1e-5);
    ASSERT_NEAR(drm2(2,0), -1, 1e-5);
    EXPECT_NEAR_RELATIVE(drm2(2,2), -5.34522484e-08, 1e-5);

    // drm3
    EXPECT_NEAR_RELATIVE(drm3(0,0), -8.01783726e-08, 1e-5);
    EXPECT_NEAR_RELATIVE(drm3(0,2), 1.33630621e-08, 1e-5);
    EXPECT_NEAR_RELATIVE(drm3(1,1), -8.01783726e-08, 1e-5);
    EXPECT_NEAR_RELATIVE(drm3(2,0), 1.33630621e-08, 1e-5);
    ASSERT_NEAR(drm3(2,2), 0., 1e-5);
}

TEST(RotMatDerivs, SmallTheta_Works)
{
    VecN<3> p;
    for (size_t i = 0; i < p.size(); ++i) p[i] = i+1;
    p /= norm<3>(p) * 1e7;
    MatrixNM<3,3> mx(1.);
    MatrixNM<3,3> drm1(0.);
    MatrixNM<3,3> drm2(0.);
    MatrixNM<3,3> drm3(0.);

    MatrixNM<3,3> mx_2(1.);
    MatrixNM<3,3> drm1_2(0.);
    MatrixNM<3,3> drm2_2(0.);
    MatrixNM<3,3> drm3_2(0.);

    pele::rot_mat_derivatives(p, mx, drm1, drm2, drm3);
    pele::rot_mat_derivatives_small_theta(p, mx_2, drm1_2, drm2_2, drm3_2, true);

    for (size_t i=0; i<9; ++i) {
        ASSERT_DOUBLE_EQ(mx.data()[i], mx_2.data()[i]);
        ASSERT_DOUBLE_EQ(drm1.data()[i], drm1_2.data()[i]);
        ASSERT_DOUBLE_EQ(drm2.data()[i], drm2_2.data()[i]);
        ASSERT_DOUBLE_EQ(drm3.data()[i], drm3_2.data()[i]);
    }
}

TEST(RotMat, SmallTheta_Works)
{
    VecN<3> p;
    for (size_t i = 0; i < p.size(); ++i) p[i] = i+1;
    p /= norm<3>(p) * 1e7;

    auto mx = pele::aa_to_rot_mat(p);

    MatrixNM<3,3> mx_2(1.);
    MatrixNM<3,3> drm1_2(0.);
    MatrixNM<3,3> drm2_2(0.);
    MatrixNM<3,3> drm3_2(0.);
    pele::rot_mat_derivatives_small_theta(p, mx_2, drm1_2, drm2_2, drm3_2, false);

    for (size_t i=0; i<9; ++i) {
        ASSERT_DOUBLE_EQ(mx.data()[i], mx_2.data()[i]);
    }
}



class AATopologyTest :  public ::testing::Test {
public:
    Array<double> x0;
    size_t nrigid;
    pele::RBTopology rbtopology;

    Array<double> make_otp_x()
    {
        Array<double> x(3*3);
        x[0] = 0;
        x[1] = -0.52890223;
        x[2] = 0;
        x[3] = 0.60876143;
        x[4] = 0.26445111;
        x[5] = 0;
        x[6] = -0.60876143;
        x[7] = 0.26445111;
        x[8] = 0;
        return x;
    }

    virtual void SetUp(){
        nrigid = 3;
        x0 = Array<double>(nrigid*6);
        for (size_t i = 0; i < x0.size(); ++i) {
            x0[i] = i;
        }
        for (size_t i = 0; i < nrigid; ++i) {
            auto p = x0.view(3*nrigid + 3*i, 3*nrigid + 3*i + 3);
            p /= norm(p);
        }
        for (size_t i = 0; i < nrigid; ++i) {
            rbtopology.add_site(make_otp_x());
        }

    }
};


TEST_F(AATopologyTest, CoordsAdaptorGetRBPositions_Works)
{
    auto ca = CoordsAdaptor(nrigid, 0, x0);
    Array<double> aapos = ca.get_rb_positions();
    ASSERT_EQ(aapos.size(), nrigid*3);
    for (size_t i=0; i<aapos.size(); ++i){
        ASSERT_DOUBLE_EQ(aapos[i], x0[i]);
    }
    // change aapos and make sure x0 changes also
    aapos[0] = 100;
    ASSERT_DOUBLE_EQ(aapos[0], x0[0]);
}

TEST_F(AATopologyTest, CoordsAdaptorGetRBRotations_Works)
{
    auto ca = CoordsAdaptor(nrigid, 0, x0);
    Array<double> aarot = ca.get_rb_rotations();
    ASSERT_EQ(aarot.size(), nrigid*3);
    for (size_t i=0; i<aarot.size(); ++i){
        ASSERT_DOUBLE_EQ(aarot[i], x0[i + 3*nrigid]);
    }
    // change aapos and make sure x0 changes also
    aarot[0] = 100;
    ASSERT_DOUBLE_EQ(aarot[0], x0[0 + 3*nrigid]);
}

TEST_F(AATopologyTest, CoordsAdaptorGetAtomPositions_Works)
{
    auto ca = CoordsAdaptor(nrigid, 0, x0);
    Array<double> pos = ca.get_atom_positions();
    ASSERT_EQ(pos.size(), 0u);
}

TEST_F(AATopologyTest, SiteToAtomistic_Works)
{
    auto rf = pele::RigidFragment(make_otp_x());
    Array<double> com(3);
    VecN<3> p;
    for (size_t i=0; i<3; ++i){
        com[i] = i+4;
        p[i] = i+1;
    }
    p /= norm<3>(p);
    auto pos = rf.to_atomistic(com, p);

    ASSERT_NEAR(pos[0], 4.32210497, 1e-5);
    ASSERT_NEAR(pos[1], 4.64476573, 1e-5);
    ASSERT_NEAR(pos[2], 5.77685304, 1e-5);
    ASSERT_NEAR(pos[3], 4.18785174, 1e-5);
    ASSERT_NEAR(pos[4], 5.62831296, 1e-5);
    ASSERT_NEAR(pos[5], 5.89772867, 1e-5);
    ASSERT_NEAR(pos[6], 3.4900433, 1e-5);
    ASSERT_NEAR(pos[7], 4.72692132, 1e-5);
    ASSERT_NEAR(pos[8], 6.32541829, 1e-5);
}

TEST_F(AATopologyTest, ToAtomisticOneMolecule_Works)
{
    pele::RBTopology rbtop;
    rbtop.add_site(make_otp_x());
    pele::RigidFragment rf(make_otp_x());
    Array<double> com(3);
    VecN<3> p;
    Array<double> x(6);
    for (size_t i=0; i<3; ++i){
        com[i] = i+4;
        p[i] = i+1;
    }
    p /= norm<3>(p);
    for (size_t i=0; i<3; ++i){
        x[i] = com[i];
        x[i+3] = p[i];
    }
    auto pos1 = rbtop.to_atomistic(x);
    auto pos2 = rf.to_atomistic(com, p);

    for (size_t i = 0; i < pos1.size(); ++i){
        ASSERT_DOUBLE_EQ(pos1[i], pos2[i]);
    }
}

TEST_F(AATopologyTest, ToAtomistic_Works)
{
//    std::cout << "x0 " << x0 << "\n";
    auto x = rbtopology.to_atomistic(x0);
//    std::cout << x << "\n";
    ASSERT_EQ(x.size(), 27u);
    ASSERT_NEAR(x[0], 0.20925348, 1e-5);
    ASSERT_NEAR(x[2], 1.68095005, 1e-5);
    ASSERT_NEAR(x[4], 1.59078223, 1e-5);
    ASSERT_NEAR(x[14], 4.95902545181, 1e-5);
    ASSERT_NEAR(x[23], 7.9605436832, 1e-5);
    ASSERT_NEAR(x[26], 8.36592352, 1e-5);
}

TEST_F(AATopologyTest, SiteTransformGrad_Works)
{
    auto rf = pele::RigidFragment(make_otp_x());
    VecN<3> p;
    for (size_t i=0; i<3; ++i){
        p[i] = i+1;
    }
    p /= norm<3>(p);
    Array<double> g = x0.view(0,9).copy();
    VecN<3> g_com;
    VecN<3> g_rot;
//    std::cout << g << "\n";
    rf.transform_grad(p, g, g_com, g_rot);

    ASSERT_NEAR(g_com[0], 9., 1e-5);
    ASSERT_NEAR(g_com[1], 12., 1e-5);
    ASSERT_NEAR(g_com[2], 15., 1e-5);

    ASSERT_NEAR(g_rot[0], 1.00790482, 1e-5);
    ASSERT_NEAR(g_rot[1], 3.63361269, 1e-5);
    ASSERT_NEAR(g_rot[2], -3.20618434, 1e-5);
}

TEST_F(AATopologyTest, TransformGradient_Works)
{
    auto x = rbtopology.to_atomistic(x0);
    auto lj = pele::LJ(4., 4.);
    Array<double> g_atom(rbtopology.natoms_total() * 3);
    lj.get_energy_gradient(x, g_atom);
    Array<double> grb(x0.size());
    rbtopology.transform_gradient(x0, g_atom, grb);

    ASSERT_EQ(grb.size(), 18u);
    ASSERT_NEAR(grb[0], -1.45358337e-03, 1e-8);
    ASSERT_NEAR(grb[2], -1.54473759e-03, 1e-8);
    ASSERT_NEAR(grb[4], -6.72613718e-06, 1e-8);
    ASSERT_NEAR(grb[8], 1.55119284e-03, 1e-8);
    ASSERT_NEAR(grb[15], -2.32088836e-04, 1e-8);
    ASSERT_NEAR(grb[17], 6.21604179e-05, 1e-8);

}

TEST_F(AATopologyTest, RBPotential_Works)
{
    pele::RBPotentialWrapper rbpot(std::make_shared<pele::LJ>(4,4));
    for (size_t i=0; i<nrigid; ++i){
        rbpot.add_site(make_otp_x());
    }
    double e = rbpot.get_energy(x0);
    ASSERT_NEAR(e, -2.55718209697, 1e-5);

    Array<double> grad(x0.size());
    double eg = rbpot.get_energy_gradient(x0, grad);
    ASSERT_DOUBLE_EQ(e, eg);

    Array<double> ngrad(x0.size());
    rbpot.numerical_gradient(x0, ngrad);

    for (size_t i=0; i<grad.size(); ++i) {
        ASSERT_NEAR(grad[i], ngrad[i], 1e-5);
    }
}


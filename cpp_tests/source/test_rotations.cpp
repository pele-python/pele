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


TEST(Rotations, AA2RotMat_Works)
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

TEST(Rotations, RotMat2q_Works){
    MatrixNM<3,3> mx;
    for (size_t i = 0; i < mx.size(); ++i) mx.data()[i] = i;
    VecN<4> q = pele::rot_mat_to_quaternion(mx);
    ASSERT_NEAR(q[0], 1.80277564, 1e-5);
    ASSERT_NEAR(q[1], 0.2773501, 1e-5);
    ASSERT_NEAR(q[2], -0.5547002, 1e-5);
    ASSERT_NEAR(q[3], 0.2773501, 1e-5);
}

TEST(Rotations, Quaternion2AA_Works){
    VecN<3> v;
    for (size_t i = 0; i < v.size(); ++i) v[i] = i+1;
    v /= norm<3>(v);

    VecN<4> q(4);
    for (size_t i = 0; i < 3; ++i) q[i+1] = v[i];
//    std::cout << q << std::endl;

    VecN<3> aa = pele::quaternion_to_aa(q);
    ASSERT_NEAR(aa[0], 0.1309466, 1e-5);
    ASSERT_NEAR(aa[1], 0.26189321, 1e-5);
    ASSERT_NEAR(aa[2], 0.39283981, 1e-5);
}

TEST(Rotations, RotMat2aa_Works){
    MatrixNM<3,3> mx;
    for (size_t i = 0; i < mx.size(); ++i) mx.data()[i] = i;
    VecN<3> p = pele::rot_mat_to_aa(mx);
    ASSERT_NEAR(p[0], 0.29425463, 1e-5);
    ASSERT_NEAR(p[1], -0.58850926, 1e-5);
    ASSERT_NEAR(p[2], 0.29425463, 1e-5);
}

TEST(Rotations, QuaternionMultiply_Works){
    VecN<4> q1, q2;
    for (size_t i = 0; i < q1.size(); ++i) q1[i] = i+1;
    for (size_t i = 0; i < q2.size(); ++i) q2[i] = i+2;

    VecN<4> q3 = pele::quaternion_multiply(q1, q2);
    ASSERT_NEAR(q3[0], -36., 1e-5);
    ASSERT_NEAR(q3[1], 6., 1e-5);
    ASSERT_NEAR(q3[2], 12., 1e-5);
    ASSERT_NEAR(q3[3], 12., 1e-5);
}

TEST(Rotations, AA2Q_Works)
{
    VecN<3> p;
    for (size_t i = 0; i < p.size(); ++i) p[i] = i+1;
    p /= norm<3>(p);
    VecN<4> q = pele::aa_to_quaternion(p);
    ASSERT_NEAR(q[0], 0.87758256, 1e-5);
    ASSERT_NEAR(q[1], 0.12813186, 1e-5);
    ASSERT_NEAR(q[2], 0.25626373, 1e-5);
    ASSERT_NEAR(q[3], 0.38439559, 1e-5);
}

TEST(Rotations, AA2QAndBack_Works)
{
    VecN<3> p;
    for (size_t i = 0; i < p.size(); ++i) p[i] = i+1;
    p /= norm<3>(p);
    VecN<3> pnew = pele::quaternion_to_aa(pele::aa_to_quaternion(p));

    ASSERT_NEAR(p[0], pnew[0], 1e-5);
    ASSERT_NEAR(p[1], pnew[1], 1e-5);
    ASSERT_NEAR(p[2], pnew[2], 1e-5);
}

TEST(Rotations, RotateAA_Works)
{
    VecN<3> p1, p2;
    for (size_t i = 0; i < p1.size(); ++i) p1[i] = i+1;
    p2 = p1;
    p2 += 1;
    VecN<3> p3 = pele::rotate_aa(p1, p2);

    ASSERT_NEAR(p3[0], 0.74050324, 1e-5);
    ASSERT_NEAR(p3[1], 1.64950785, 1e-5);
    ASSERT_NEAR(p3[2], 2.20282887, 1e-5);
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

TEST(RotMatDerivs, SmallTheta_Works)
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

TEST(RotMatDerivs, SmallTheta2_Works)
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

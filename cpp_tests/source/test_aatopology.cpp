#include <gtest/gtest.h>
#include <limits>

#include <memory>
#include "pele/aatopology.h"
#include "pele/lj.h"
#include "pele/vecn.h"
#include "pele/matrix.h"
#include "pele/distance.h"

using pele::Array;
using pele::VecN;
using pele::MatrixNM;
using pele::CoordsAdaptor;
using pele::norm;

static double const EPS = std::numeric_limits<double>::min();
#define EXPECT_NEAR_RELATIVE(A, B, T)  ASSERT_NEAR(A/(fabs(A)+fabs(B) + EPS), B/(fabs(A)+fabs(B) + EPS), T)

class AATopologyTest :  public ::testing::Test {
public:
    Array<double> x0;
    size_t nrigid;
    std::shared_ptr<pele::RBTopology> rbtopology;
    VecN<3> p0;


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

    pele::RigidFragment make_otp()
    {
        auto x = make_otp_x();
        Array<double> cog(3,0);
        double M = 3.;
        double W = 3.;
        Array<double> S(9,0);
        S[0] = 0.74118095;
        S[4] = 0.41960635;
        bool invertible = true;
        Array<double> eye(9,0);
        eye[0] = 1;
        eye[4] = 1;
        eye[8] = 1;
        Array<double> inversion = eye.copy();
        inversion[4] = -1;
        inversion[8] = -1;
        auto distance_function = std::make_shared<pele::CartesianDistanceWrapper<3> >();
        pele::RigidFragment rf(x, cog, M, W, S, inversion, invertible, distance_function);

        auto rot = eye.copy();
        rot[0] = -1;
        rot[8] = -1;
        rf.add_symmetry_rotation(rot);

        rf.add_symmetry_rotation(eye.copy());

        return rf;
    }

    void add_otp(pele::RBTopology * top)
    {
        top->add_site(make_otp());
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
        for (size_t i = 0; i < 3; ++i) p0[i] = i+1;
        p0 /= norm(p0);
        rbtopology = std::make_shared<pele::RBTopology>();
        for (size_t i = 0; i < nrigid; ++i) {
            add_otp(rbtopology.get());
        }
        rbtopology->assign_atom_indices();

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
    pele::RigidFragment rf = make_otp();
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
    add_otp(&rbtop);
    rbtop.assign_atom_indices();
    pele::RigidFragment rf = make_otp();
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
    auto x = rbtopology->to_atomistic(x0);
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
    auto rf = make_otp();
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
    auto x = rbtopology->to_atomistic(x0);
    auto lj = pele::LJ(4., 4.);
    Array<double> g_atom(rbtopology->natoms_total() * 3);
    lj.get_energy_gradient(x, g_atom);
    Array<double> grb(x0.size());
    rbtopology->transform_gradient(x0, g_atom, grb);

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
    pele::RBPotentialWrapper rbpot(std::make_shared<pele::LJ>(4,4), rbtopology);
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

TEST_F(AATopologyTest, TransformRotate_Works)
{
    auto x = x0.copy();
    pele::TransformAACluster transform(rbtopology.get());

    VecN<3> p;
    for (size_t i = 0; i < p.size(); ++i) p[i] = i+1;
    p /= norm<3>(p);

    transform.rotate(x, pele::aa_to_rot_mat(p));
//    std::cout << p << std::endl;
//    std::cout << pele::aa_to_rot_mat(p) << std::endl;
//    std::cout << x0 << std::endl;
//    std::cout << x << std::endl;

    ASSERT_NEAR(x[0], 0.48757698, 1e-5);
    ASSERT_NEAR(x[4], 4.76822812, 1e-5);
    ASSERT_NEAR(x[8], 7.53224809, 1e-5);
    ASSERT_NEAR(x[12], 0.72426504, 1e-5);
    ASSERT_NEAR(x[16], 1.25159032, 1e-5);
}

TEST_F(AATopologyTest, AlignPath_Works)
{
    auto x1 = x0.copy();
    auto x2 = x1.copy();
    x2 += 5;
    auto x20 = x2.copy();

    std::list<Array<double> > path;
    path.push_back(x1);
    path.push_back(x2);
    rbtopology->align_path(path);
//    std::cout << x2 << std::endl;

    for (size_t i = 0; i < x1.size(); ++i) {
        ASSERT_EQ(x1[i], x0[i]);
    }

    ASSERT_NEAR(x2[9], 1.92786071, 1e-5);
    ASSERT_NEAR(x2[13], 1.94869267, 1e-5);
    ASSERT_NEAR(x2[17], 1.96164668, 1e-5);
}

TEST_F(AATopologyTest, ZeroEV_Works)
{
    std::vector<pele::Array<double> > zev;
    rbtopology->get_zero_modes(x0, zev);
    ASSERT_EQ(zev.size(), 6u);
    for (auto const & v : zev) {
        ASSERT_EQ(v.size(), x0.size());
    }

    // check translational
    for (size_t i = 0; i<3; ++i) {
        auto v = zev[i];
        for (size_t j = 0; j<v.size(); ++j) {
            if ( j<3*nrigid && j % 3 == i ){
                ASSERT_NEAR(v[j], 1./sqrt(nrigid), 1e-5);
            } else {
                ASSERT_EQ(v[j], 0);
            }
        }
    }

    pele::Array<double> v;

    v = zev[3];
    ASSERT_NEAR(v[0], 0., 1e-5);
    ASSERT_NEAR(v[4], -0.39260962, 1e-5);
    ASSERT_NEAR(v[13], -0.0223232, 1e-5);
    ASSERT_NEAR(v[17], 0.02484174, 1e-5);

    v = zev[4];
    ASSERT_NEAR(v[0], 1.68325526e-01, 1e-5);
    ASSERT_NEAR(v[4], 0., 1e-5);
    ASSERT_NEAR(v[13], 7.93978746e-02, 1e-5);
    ASSERT_NEAR(v[17], -2.02278328e-02, 1e-5);

    v = zev[5];
    ASSERT_NEAR(v[0], -9.35924237e-02, 1e-5);
    ASSERT_NEAR(v[4],  2.80775399e-01, 1e-5);
    ASSERT_NEAR(v[13], 2.77268394e-02, 1e-5);
    ASSERT_NEAR(v[17], 8.86371673e-02, 1e-5);
}

TEST_F(AATopologyTest, SiteDistance_Works)
{
    pele::RigidFragment rf = make_otp();
    VecN<3> com1(0);
    VecN<3> com2(1);
    auto p1 = p0;
    auto p2 = p0;
    p2 += 1;
    double dist2 = rf.distance_squared(com1, p1, com2, p2);
    ASSERT_NEAR(dist2, 10.9548367929, 1e-5);

}

TEST_F(AATopologyTest, DistanceSquared_Works)
{
    auto x1 = x0.copy();
    auto x2 = x0.copy();
    x2 += 1.1;
    double d2 = rbtopology->distance_squared(x1, x2);
    ASSERT_NEAR(d2, 38.9401810973, 1e-5);
}

TEST_F(AATopologyTest, DistanceSquaredGrad_Works)
{
    auto x1 = x0.copy();
    auto x2 = x0.copy();
    x2 += 1.1;
    auto grad = x1.copy();
    rbtopology->distance_squared_grad(x1, x2, grad);
//    std::cout << grad << std::endl;

    for (size_t i =0; i < nrigid*3; ++i) {
        ASSERT_DOUBLE_EQ(grad[i], -6.6);
    }

    ASSERT_NEAR(grad[nrigid*3 + 0], -1.21579025, 1e-5);
    ASSERT_NEAR(grad[nrigid*3 + 4], -0.06984532, 1e-5);
    ASSERT_NEAR(grad[nrigid*3 + 8], -1.28362943, 1e-5);

}

TEST_F(AATopologyTest, MeasureAlign_Works)
{
    auto x1 = x0.copy();
    auto x2 = x0.copy();
    x2 += 5.1;
    x2[x2.size()-1] = x1[x2.size()-1] + .1;
    auto x20 = x2.copy();
    pele::MeasureAngleAxisCluster measure(rbtopology.get());
    measure.align(x1, x2);

    // assert the com coordinates didn't change
    for (size_t i =0; i < nrigid*3; ++i) {
        ASSERT_DOUBLE_EQ(x2[i], x20[i]);
    }

    ASSERT_NEAR(x2[nrigid*3 + 0], 1.34332528, 1e-5);
    ASSERT_NEAR(x2[nrigid*3 + 4], 0.47517882, 1e-5);
    ASSERT_NEAR(x2[nrigid*3 + 8], 0.15531234, 1e-5);
}



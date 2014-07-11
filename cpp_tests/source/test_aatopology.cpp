#include <gtest/gtest.h>

#include "pele/aatopology.h"

using pele::Array;
using pele::CoordsAdaptor;


TEST(AA2RotMat, Works)
{
    Array<double> p(3);
    for (size_t i = 0; i < p.size(); ++i) p[i] = i+1;
    p /= norm(p);
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


class AATopologyTest :  public ::testing::Test {
public:
    Array<double> x0;
    size_t nrigid;

    virtual void SetUp(){
        nrigid = 3;
        x0 = Array<double>(nrigid*6);
        for (size_t i = 0; i < x0.size(); ++i) {
            x0[i] = i;
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
    ASSERT_EQ(pos.size(), 0);
}

#include <gtest/gtest.h>

#include "pele/aatopology.h"

using pele::Array;
using pele::CoordsAdaptor;

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

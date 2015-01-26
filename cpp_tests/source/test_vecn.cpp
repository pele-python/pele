#include "pele/vecn.h"
#include <iostream>
#include <stdexcept>
#include <vector>
#include <gtest/gtest.h>
#include "pele/aatopology.h"
#include "pele/matrix.h"

using pele::VecN;
using pele::MatrixNM;

// to silence "unused variable" warnings from the compiler
#define UNUSED(var) ((void) var)

TEST(VecNTest, DefaultConstructor_IsOK){
    pele::VecN<3> v;
    EXPECT_EQ(3u, v.size());
}

TEST(VecNTest, IteratorConstructor_IsOK){
    std::vector<double> v(3);
    v[0] = 0; v[1] = 1; v[2] = 2.1;
    pele::VecN<3> v2(v.begin(), v.end());
    for (size_t i=0; i<3; ++i) {
        ASSERT_DOUBLE_EQ(v2[i], v[i]);
    }
}

TEST(VecNTest, Constructor_SetsValue){
    pele::VecN<3> v(10.);
    for (int i = 0; i < 3; ++i){
        EXPECT_EQ(v[i], 10.);
    }
}

TEST(VecNTest, CopyConstructor_CoppiesCorrectly){
    pele::VecN<3> v(1.);
    pele::VecN<3> v2(v);
    EXPECT_NE(v.data(), v2.data());
    for (int i = 0; i < 3; ++i){
        EXPECT_EQ(v[i], v2[i]);
    }
}

TEST(VecNTest, AssignmentOperator_CoppiesCorrectly){
    // should wrap v
    pele::VecN<3> v(6);
    pele::VecN<3> v2;
    double * data = v2.data();
    v2 = v;
    EXPECT_NE(v.data(), v2.data());
    EXPECT_EQ(v2.data(), data);
    for (int i = 0; i < 3; ++i){
        EXPECT_EQ(v[i], v2[i]);
    }
}

TEST(VecNTest, RangeBasedFor_Works){
    VecN<3> v(0);
    for (double & val : v){
        val += 1;
    }
    for (size_t i=0; i<v.size(); ++i){
        EXPECT_EQ(1,v[i]);
    }
}

TEST(VecNTest, SumOperator_VecN){
    pele::VecN<3> v(-1);
    pele::VecN<3> v2(1);
    for (size_t i=0; i<v2.size(); ++i) v2[i] = 1;
    v += v2;
    EXPECT_NE(v.data(), v2.data());
    for (int i = 0; i < 3; ++i){
        EXPECT_EQ(v[i], 0);
        EXPECT_EQ(v2[i], 1);
    }
}

TEST(VecNTest, SumOperator_VecNSelf){
    pele::VecN<3> v(1);
    v += v;
    for (int i = 0; i < 3; ++i){
        EXPECT_EQ(v[i], 2);
    }
}

TEST(VecNTest, SumOperator_Const){
    double c = 1;
    pele::VecN<3> v(-1);
    v += c;
    for (int i = 0; i < 3; ++i){
        EXPECT_EQ(v[i], 0);
    }
    EXPECT_EQ(c,1);
}

///////////

TEST(VecNTest, DifOperator_VecN){
    pele::VecN<3> v(1);
    pele::VecN<3> v2(1);
    v -= v2;
    EXPECT_NE(v.data(), v2.data());
    for (int i = 0; i < 3; ++i){
        EXPECT_EQ(v[i], 0);
        EXPECT_EQ(v2[i], 1);
    }
}

TEST(VecNTest, DifOperator_VecNSelf){
    pele::VecN<3> v(1);
    v -= v;
    for (int i = 0; i < 3; ++i){
        EXPECT_EQ(v[i], 0);
    }
}

TEST(VecNTest, DifOperator_Const){
    double c = 1;
    pele::VecN<3> v(1);
    v -= c;
    for (int i = 0; i < 3; ++i){
        EXPECT_EQ(v[i], 0);
    }
    EXPECT_EQ(c,1);
}

///////////

TEST(VecNTest, ProdOperator_VecN){
    pele::VecN<3> v(2);
    pele::VecN<3> v2(10);
    v *= v2;
    EXPECT_NE(v.data(), v2.data());
    for (int i = 0; i < 3; ++i){
        EXPECT_EQ(v[i], 20);
        EXPECT_EQ(v2[i], 10);
    }
}

TEST(VecNTest, ProdOperator_VecNSelf){
    pele::VecN<3> v(2);
    v *= v;
    for (int i = 0; i < 3; ++i){
        EXPECT_EQ(v[i], 4);
    }
}

TEST(VecNTest, ProdOperator_Const){
    double c = 10;
    pele::VecN<3> v(2);
    v *= c;
    for (int i = 0; i < 3; ++i){
        EXPECT_EQ(v[i], 20);
    }
    EXPECT_EQ(c, 10);
}

///////////

TEST(VecNTest, DivOperator_VecN){
    pele::VecN<3> v(2);
    pele::VecN<3> v2(2);
    v /= v2;
    for (int i = 0; i < 3; ++i){
        EXPECT_EQ(v[i], 1);
        EXPECT_EQ(v2[i], 2);
    }
}

TEST(VecNTest, DivOperator_VecNSelf){
    pele::VecN<3> v(2);
    v /= v;
    for (int i = 0; i < 3; ++i){
        EXPECT_EQ(v[i], 1);
    }
}

TEST(VecNTest, DivOperator_Const){
    double c = 2;
    pele::VecN<3> v(2);
    v /= c;
    for (int i = 0; i < 3; ++i){
        EXPECT_EQ(v[i], 1);
    }
    EXPECT_EQ(c,2);
}

//////////////

TEST(VecNTest, DifferenceOperator_Works){
    pele::VecN<3> v1(1);
    pele::VecN<3> v2(-1);
    pele::VecN<3> v3 = v1 - v2;
    for (int i = 0; i < 3; ++i){
        EXPECT_EQ(v3[i], 2);
    }
}


TEST(VecNTest, SumFunction){
    pele::VecN<3> v(2);
    EXPECT_EQ(v.sum(), 6);
}

TEST(VecNTest, ProdFunction){
    pele::VecN<3> v(2);
    EXPECT_EQ(v.prod(), 8);
}

TEST(VecNTest, Iterator_Works){
    pele::VecN<3> v(-1);
    size_t count = 0;
    for (VecN<3>::iterator iter = v.begin(); iter != v.end(); ++iter){
        *iter = count;
        count++;
    }
    for (count=0; count<v.size(); ++count){
        EXPECT_EQ(v[count], count);
    }
}

//
//TEST(VecNTest, ConstVecN_NotModifiable){
//    pele::VecN<3> const v(6,0);
//
//    double x = v[0];
//    v.empty();
//    double const * d = v.data();
//    v.data();
//    size_t count = 0;
//    for (VecN<3>::const_iterator iter = v.begin(); iter != v.end(); ++iter){
//        x = *iter;
//        count++;
//    }
//    EXPECT_EQ(count, v.size());
//    UNUSED(d);
//    UNUSED(x);
//}
//

TEST(VecNDotTest, Dot_Works){
    pele::VecN<3> v1(5);
    pele::VecN<3> v2(2);

    double d = pele::dot<3>(v1, v2);
    EXPECT_EQ(5*2*3, d);
}

TEST(VecNNormTest, Norm_Works){
    pele::VecN<3> v1(5);

    double d = pele::norm<3>(v1);
    EXPECT_NEAR(sqrt(5*5*3), d, 1e-10);
}

TEST(MatrixNMTest, MatrixMultiplication_Works){
    MatrixNM<2,3> v1;
    MatrixNM<3,2> v2;
    for (size_t i = 0; i < 3*2; ++i){
        v1.data()[i] = i;
        v2.data()[i] = i+1;
    }
//    std::cout << v1 << std::endl;
//    std::cout << v2 << std::endl;
    MatrixNM<2,2> v3 = pele::dot(v1, v2);
//    std::cout << v3 << std::endl;
    EXPECT_EQ(v3(0,0), 13);
    EXPECT_EQ(v3(1,0), 40);
    EXPECT_EQ(v3(0,1), 16);
    EXPECT_EQ(v3(1,1), 52);
}


TEST(MatrixNMTest, MatrixVectorMultiplication_Works){
    MatrixNM<2,3> v1;
    VecN<3> v2;
    for (size_t i = 0; i < 3*2; ++i){
        v1.data()[i] = i;
    }
    for (size_t i = 0; i < 3; ++i){
        v2.data()[i] = i+1;
    }
//    std::cout << v1 << std::endl;
//    std::cout << v2 << std::endl;
    VecN<2> v3 = pele::dot(v1, v2);
//    std::cout << v3 << std::endl;
    EXPECT_EQ(v3[0], 8);
    EXPECT_EQ(v3[1], 26);
}


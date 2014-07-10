#include "pele/array.h"
#include <iostream>
#include <stdexcept>
#include <vector>
#include <gtest/gtest.h>

using pele::Array;

// to silence "unused variable" warnings from the compiler
#define UNUSED(var) ((void) var)

TEST(ArrayTest, DefaultConstructor_IsOK){
    pele::Array<double> v;
    EXPECT_EQ(0u, v.size());
    EXPECT_FALSE(v.data());
    EXPECT_TRUE(v.empty());
}

TEST(ArrayTest, Constructor_HandlesPositiveInput){
    pele::Array<double> v(6);
    EXPECT_EQ (6u, v.size());
    EXPECT_TRUE(v.data());
}

TEST(ArrayTest, Constructor_HandlesZero){
    pele::Array<double> v(0);
    EXPECT_EQ (0u, v.size());
    EXPECT_TRUE(v.empty());
}

TEST(ArrayTest, Constructor_SetsValue){
    pele::Array<double> v(6u, 10.);
    EXPECT_EQ (6u, v.size());
    for (int i = 0; i < 6; ++i){
        EXPECT_NEAR(v[i], 10., 1e-10);
    }
}

TEST(ArrayTest, CopyConstructor_WrapsCorrectly){
    // should wrap v
    pele::Array<double> v(6);
    pele::Array<double> v2(v);
    EXPECT_EQ(v.data(), v2.data());
    EXPECT_EQ(v.size(), v2.size());
}

TEST(ArrayTest, FreeWrappedArray_Works){
    // I discovered an error with this test
    pele::Array<double> v(6);
    pele::Array<double> v2(v);
    pele::Array<double> v3(v);
    EXPECT_EQ(3, v.reference_count());
    EXPECT_EQ(3, v2.reference_count());
    v.free();
    EXPECT_EQ(1, v.reference_count());
    EXPECT_EQ(2, v2.reference_count());
    EXPECT_EQ(2, v3.reference_count());
    //    ASSERT_THROW(pele::Array<double> v3(v), std::runtime_error);
}


//TEST(ArrayTest, CopyConstructorEmptyArray_Fails){
//    // should wrap v
//    pele::Array<double> v;
//    ASSERT_THROW(pele::Array<double> v2(v), std::runtime_error);
//}

TEST(ArrayTest, AssignmentOperator_WrapsCorrectly){
    // should wrap v
    pele::Array<double> v(6);
    pele::Array<double> v2;
    v2 = v;
    EXPECT_EQ(v.data(), v2.data());
    EXPECT_EQ(v.size(), v2.size());
}

//TEST(ArrayTest, AssignmentOperatorEmptyArray_Fails){
//    // This throwing an error is consistent with the logic, but I wish it didn't fail
//    // Maybe someday we can change the logic so this won't have to fail
//    pele::Array<double> v;
//    ASSERT_THROW(v = Array<double>(), std::runtime_error);
//    pele::Array<double> v2;
//    ASSERT_THROW(v = v2, std::runtime_error);
//}

TEST(ArrayTest, Wrap_Works){
    pele::Array<double> v(6);
    EXPECT_EQ(v.reference_count(), 1);
    pele::Array<double> v2;
    v2.wrap(v);
    EXPECT_EQ(v.data(), v2.data());
    EXPECT_EQ(v.reference_count(), 2);
    EXPECT_EQ(v2.reference_count(), v.reference_count());
}

//TEST(ArrayTest, WrapEmptyArray_Fails){
//    // can't wrap an array with no data
//    pele::Array<double> v;
//    pele::Array<double> v2;
//    ASSERT_THROW(v2.wrap(v), std::runtime_error);
//}

TEST(ArrayTest, WrapSelf_DoesNothing){
    pele::Array<double> v(6);
    EXPECT_EQ(v.reference_count(), 1u);
    v.wrap(v);
    EXPECT_EQ(v.reference_count(), 1u);
}

TEST(ArrayTest, Free_Works){
    pele::Array<double> v(6);
    v.free();
    EXPECT_EQ (0u, v.size());
//    EXPECT_FALSE(v.data());
    EXPECT_TRUE(v.empty());
    v.free();
}

TEST(ArrayTest, Free_Unwraps){
    pele::Array<double> v1(6);
    pele::Array<double> v2(v1);
    EXPECT_EQ(2, v1.reference_count());
    EXPECT_TRUE(v1 == v2);
    v2.free();
    EXPECT_EQ(1u, v1.reference_count());
    EXPECT_EQ(1u, v2.reference_count());
    EXPECT_TRUE(v1 != v2);
    EXPECT_EQ(0u, v2.size());
    EXPECT_EQ(6u, v1.size());
}

TEST(ArrayTest, Copy_WorksNotWrap){
    pele::Array<double> v(6);
    for (size_t i=0; i<v.size(); ++i) v[i] = i;
    pele::Array<double> v2(6);
    v2 = v.copy();
    EXPECT_NE(v.data(), v2.data());
    for (int i = 0; i < 6; ++i){
        EXPECT_NEAR(v[i], v2[i], 1e-10);
    }
}

TEST(ArrayTest, Assign_CopiesNoWrap){
    pele::Array<double> v(6);
    for (size_t i=0; i<v.size(); ++i) v[i] = i;
    pele::Array<double> v2(6);
    v2.assign(v);
    EXPECT_NE(v.data(), v2.data());
    for (int i = 0; i < 6; ++i){
        EXPECT_NEAR(v[i], v2[i], 1e-10);
    }
}

TEST(ArrayTest, AssignWrongSize_Fails){
    pele::Array<double> v(6);
    pele::Array<double> v2;
    EXPECT_THROW(v2.assign(v), std::runtime_error);
    v2 = Array<double>(1);
    EXPECT_THROW(v2.assign(v), std::runtime_error);
    v2 = Array<double>(8);
    EXPECT_THROW(v2.assign(v), std::runtime_error);
    v2.free();
    v.free();
    v2.assign(v); //this is pointless, but it shouldn't fail
}

TEST(ArrayTest, ConstructFromVector_Wraps){
    std::vector<double> vec(6);
    pele::Array<double> v(vec);
    EXPECT_FALSE(v.empty());
    EXPECT_EQ(v.reference_count(), 1);
    EXPECT_EQ(v.data(), vec.data());
    EXPECT_EQ(v.size(), vec.size());
    pele::Array<double> v2(v);
    EXPECT_EQ(v.data(), v2.data());
    v2.free();
    v2 = v;
    EXPECT_EQ(v.data(), v2.data());
    v2.free();
    v2.wrap(v);
    EXPECT_EQ(v.data(), v2.data());

    v.free();
    v2.free();
    vec[0] = 0; // the data in vec should not have been deallocated
}

TEST(ArrayTest, RangeBasedFor_Works){
    Array<double> v(3,0);
    for (double & val : v){
        val += 1;
    }
    for (size_t i=0; i<v.size(); ++i){
        EXPECT_EQ(1,v[i]);
    }
}

TEST(ArrayTest, EqualityOperator_Works){
    Array<double> v1(3);
    Array<double> v2(3);
    EXPECT_FALSE(v1 == v2);
    v1.wrap(v2);
    EXPECT_TRUE(v1 == v2);
}

TEST(ArrayTest, InEqualityOperator_Works){
    Array<double> v1(3);
    Array<double> v2(3);
    EXPECT_TRUE(v1 != v2);
    v1.wrap(v2);
    EXPECT_FALSE(v1 != v2);
}

TEST(ArrayTest, SumOperator_Array){
    pele::Array<double> v(6,-1);
    pele::Array<double> v2(6,1);
    for (size_t i=0; i<v2.size(); ++i) v2[i] = 1;
    v += v2;
    EXPECT_NE(v.data(), v2.data());
    for (int i = 0; i < 6; ++i){
        EXPECT_EQ(v[i], 0);
        EXPECT_EQ(v2[i], 1);
    }
}

TEST(ArrayTest, SumOperator_ArraySelf){
    pele::Array<double> v(6,1);
    double* old_v_data = v.data();

    v += v;

    EXPECT_EQ(v.data(), old_v_data); //check that memory has not moved

    for (int i = 0; i < 6; ++i){
        EXPECT_EQ(v[i], 2);
    }
}

TEST(ArrayTest, SumOperator_Const){
    double c = 1;
    pele::Array<double> v(6, -1);
    v += c;
    for (int i = 0; i < 6; ++i){
        EXPECT_EQ(v[i], 0);
    }
    EXPECT_EQ(c,1);
}

///////////

TEST(ArrayTest, DifOperator_Array){
    pele::Array<double> v(6,1);
    pele::Array<double> v2(6,1);
    v -= v2;
    EXPECT_NE(v.data(), v2.data());
    for (int i = 0; i < 6; ++i){
        EXPECT_EQ(v[i], 0);
        EXPECT_EQ(v2[i], 1);
    }
}

TEST(ArrayTest, DifOperator_ArraySelf){
    pele::Array<double> v(6,1);
    double* old_v_data = v.data();

    v -= v;

    EXPECT_EQ(v.data(), old_v_data); //check that memory has not moved

    for (int i = 0; i < 6; ++i){
        EXPECT_EQ(v[i], 0);
    }
}

TEST(ArrayTest, DifOperator_Const){
    double c = 1;
    pele::Array<double> v(6,1);
    v -= c;
    for (int i = 0; i < 6; ++i){
        EXPECT_EQ(v[i], 0);
    }
    EXPECT_EQ(c,1);
}

///////////

TEST(ArrayTest, ProdOperator_Array){
    pele::Array<double> v(6, 2);
    pele::Array<double> v2(6 ,10);
    v *= v2;
    EXPECT_NE(v.data(), v2.data());
    for (int i = 0; i < 6; ++i){
        EXPECT_EQ(v[i], 20);
        EXPECT_EQ(v2[i], 10);
    }
}

TEST(ArrayTest, ProdOperator_ArraySelf){
    pele::Array<double> v(6, 2);

    double* old_v_data = v.data();

    v *= v;

    EXPECT_EQ(v.data(), old_v_data); //check that memory has not moved

    for (int i = 0; i < 6; ++i){
        EXPECT_EQ(v[i], 4);
    }
}

TEST(ArrayTest, ProdOperator_Const){
    double c = 10;
    pele::Array<double> v(6, 2);
    v *= c;
    for (int i = 0; i < 6; ++i){
        EXPECT_EQ(v[i], 20);
    }
    EXPECT_EQ(c, 10);
}

///////////

TEST(ArrayTest, DivOperator_Array){
    pele::Array<double> v(6,2);
    pele::Array<double> v2(6,2);
    v /= v2;
    EXPECT_NE(v.data(), v2.data());
    for (int i = 0; i < 6; ++i){
        EXPECT_EQ(v[i], 1);
        EXPECT_EQ(v2[i], 2);
    }
}

TEST(ArrayTest, DivOperator_ArraySelf){
    pele::Array<double> v(6,2);

    double* old_v_data = v.data();

    v /= v;

    EXPECT_EQ(v.data(), old_v_data); //check that memory has not moved

    for (int i = 0; i < 6; ++i){
        EXPECT_EQ(v[i], 1);
    }
}

TEST(ArrayTest, DivOperator_Const){
    double c = 2;
    pele::Array<double> v(6,2);
    v /= c;
    for (int i = 0; i < 6; ++i){
        EXPECT_EQ(v[i], 1);
    }
    EXPECT_EQ(c,2);
}

TEST(ArrayTest, SumFunction){
    pele::Array<double> v(6,2);
    EXPECT_EQ(v.sum(), 12);
}

TEST(ArrayTest, ProdFunction){
    pele::Array<double> v(6,2);
    EXPECT_EQ(v.prod(), 64);
}

TEST(ArrayTest, Iterator_Works){
    pele::Array<double> v(6,-1);
    size_t count = 0;
    for (Array<double>::iterator iter = v.begin(); iter != v.end(); ++iter){
        *iter = count;
        count++;
    }
    for (count=0; count<v.size(); ++count){
        EXPECT_EQ(v[count], count);
    }
}


TEST(ArrayTest, ConstArray_NotModifiable){
    pele::Array<double> const v(6,0);

    double x = v[0];
    v.empty();
    double const * d = v.data();
    v.data();
    size_t count = 0;
    for (Array<double>::const_iterator iter = v.begin(); iter != v.end(); ++iter){
        x = *iter;
        count++;
    }
    EXPECT_EQ(count, v.size());
    UNUSED(d);
    UNUSED(x);
}

TEST(ArrayTest, View_Works){
    pele::Array<size_t> v(6);
    for (size_t i=0; i<v.size(); ++i){
        v[i] = i;
    }
    pele::Array<size_t> v2 = v.view(1,3);
    EXPECT_EQ(v[0], 0u);
    EXPECT_EQ(v2[0], 1u);
    EXPECT_EQ(v2[1], 2u);
    EXPECT_EQ(v2.size(), 2u);
    EXPECT_FALSE(v == v2);

    pele::Array<size_t> v3 = v2.view(1,2);
    EXPECT_EQ(v3[0], 2u);
    EXPECT_EQ(v3.size(), 1u);

    EXPECT_EQ(v.size(), 6u);
    EXPECT_EQ(v[0], 0u);
}

TEST(ArrayTest, FullView_IsSame){
    pele::Array<size_t> v(6);
    pele::Array<size_t> v2 = v.view(0,6);
    EXPECT_TRUE(v == v2);
}

TEST(ArrayTest, BadInput_Fails){
    pele::Array<size_t> v(6);
    pele::Array<size_t> v2;
    ASSERT_THROW(v2 = v.view(-1,2), std::invalid_argument);
    ASSERT_THROW(v2 = v.view(3,2), std::invalid_argument);
    ASSERT_THROW(v2 = v.view(0,7), std::invalid_argument);
}

TEST(ArrayDotTest, Dot_Works){
    pele::Array<double> v1(6, 3);
    pele::Array<double> v2(6, 2);

    double d = pele::dot(v1, v2);
    EXPECT_EQ(3*2*6, d);
}

TEST(ArrayNormTest, Norm_Works){
    pele::Array<double> v1(6, 3);

    double d = pele::norm(v1);
    EXPECT_NEAR(sqrt(3*3*6), d, 1e-10);
}



#include "pele/array.h"
#include <iostream>
#include <stdexcept>
#include <vector>
#include <gtest/gtest.h>

using pele::Array;

TEST(ArrayTest, DefaultConstructor_IsOK){
    pele::Array<double> v;
    EXPECT_EQ(0, v.size());
    EXPECT_FALSE(v.data());
    EXPECT_TRUE(v.empty());
}

TEST(ArrayTest, Constructor_HandlesPositiveInput){
    pele::Array<double> v(6);
    EXPECT_EQ (6, v.size());
    EXPECT_TRUE(v.data());
}

TEST(ArrayTest, Constructor_HandlesZero){
    pele::Array<double> v(0);
    EXPECT_EQ (0, v.size());
    EXPECT_TRUE(v.empty());
}

TEST(ArrayTest, Constructor_SetsValue){
    pele::Array<double> v(6, 10.);
    EXPECT_EQ (6, v.size());
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
    v.free();
    EXPECT_EQ(0, v.reference_count());
    EXPECT_EQ(1, v2.reference_count());
//    ASSERT_THROW(pele::Array<double> v3(v), std::runtime_error);
}


TEST(ArrayTest, CopyConstructorEmptyArray_Fails){
    // should wrap v
    pele::Array<double> v;
    ASSERT_THROW(pele::Array<double> v2(v), std::runtime_error);
}

TEST(ArrayTest, AssignmentOperator_WrapsCorrectly){
    // should wrap v
    pele::Array<double> v(6);
    pele::Array<double> v2;
    v2 = v;
    EXPECT_EQ(v.data(), v2.data());
    EXPECT_EQ(v.size(), v2.size());
}

TEST(ArrayTest, AssignmentOperatorEmptyArray_Fails){
    // This throwing an error is consistent with the logic, but I wish it didn't fail
    // Maybe someday we can change the logic so this won't have to fail
    pele::Array<double> v;
    ASSERT_THROW(v = Array<double>(), std::runtime_error);
    pele::Array<double> v2;
    ASSERT_THROW(v = v2, std::runtime_error);
}

TEST(ArrayTest, Wrap_Works){
    pele::Array<double> v(6);
    EXPECT_EQ(v.reference_count(), 1);
    pele::Array<double> v2;
    v2.wrap(v);
    EXPECT_EQ(v.data(), v2.data());
    EXPECT_EQ(v.reference_count(), 2);
    EXPECT_EQ(v2.reference_count(), v.reference_count());
}

TEST(ArrayTest, WrapEmptyArray_Fails){
    // can't wrap an array with no data
    pele::Array<double> v;
    pele::Array<double> v2;
    ASSERT_THROW(v2.wrap(v), std::runtime_error);
}

TEST(ArrayTest, WrapSelf_DoesNothing){
    pele::Array<double> v(6);
    EXPECT_EQ(v.reference_count(), 1);
    v.wrap(v);
    EXPECT_EQ(v.reference_count(), 1);
}

TEST(ArrayTest, Free_Works){
    pele::Array<double> v(6);
    v.free();
    EXPECT_EQ (0, v.size());
    EXPECT_FALSE(v.data());
    EXPECT_TRUE(v.empty());
    v.free();
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
    EXPECT_EQ(v.reference_count(), 0);
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

TEST(ArrayTest, Resize_Works){
    Array<double> v(3);
    double * old_data = v.data();
    v.resize(6);
    EXPECT_NE(v.data(), old_data);
    v.resize(0);
    EXPECT_TRUE(v.empty());
    v.resize(6);
    // you can't resize a wrapped array
    Array<double> v2(v);
    EXPECT_THROW(v.resize(7), std::runtime_error);
}

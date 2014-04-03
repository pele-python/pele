#include "pele/array.h"
#include <iostream>
#include <stdexcept>
#include <gtest/gtest.h>

TEST(array_test, constructor1){
    pele::Array<double> v;
    ASSERT_EQ(0, v.size());
    ASSERT_FALSE(v.data());
}

TEST(array_test, constructor2){
    pele::Array<double> v(6);
    ASSERT_EQ (6, v.size());
    ASSERT_TRUE(v.data());
}


TEST(array_test, constructor3){
    pele::Array<double> v(6, 10.);
    ASSERT_EQ (6, v.size());
    for (int i = 0; i < 6; ++i){
        ASSERT_NEAR(v[i], 10., 1e-10);
    }
}

TEST(array_test, copy_constructor){
    // should wrap v
    pele::Array<double> v(6);
    pele::Array<double> v2(v);
    ASSERT_EQ(v.data(), v2.data());
    ASSERT_EQ(v.size(), v2.size());
}

TEST(array_test, copy_constructor_fail1){
    // I discovered an error with this test
    pele::Array<double> v(6);
    pele::Array<double> v2(v);
    v.free();
    ASSERT_EQ(0, v.reference_count());
    ASSERT_EQ(1, v2.reference_count());
//    ASSERT_THROW(pele::Array<double> v3(v), std::runtime_error);
}


TEST(array_test, copy_constructor_fail){
    // should wrap v
    pele::Array<double> v;
    ASSERT_THROW(pele::Array<double> v2(v), std::runtime_error);
    pele::Array<double> v3(6);
    pele::Array<double> v6(v3);
//    v3.free();
//    ASSERT_THROW(pele::Array<double> v4(v3), std::runtime_error);
}

TEST(array_test, assignment_operator){
    // should wrap v
    pele::Array<double> v(6);
    pele::Array<double> v2;
    v2 = v;
    ASSERT_EQ(v.data(), v2.data());
    ASSERT_EQ(v.size(), v2.size());
}

TEST(array_test, wrap_ok){
    pele::Array<double> v(6);
    ASSERT_EQ(v.reference_count(), 1);
    pele::Array<double> v2;
    v2.wrap(v);
    ASSERT_EQ(v.data(), v2.data());
    ASSERT_EQ(v.reference_count(), 2);
    ASSERT_EQ(v2.reference_count(), v.reference_count());
}

TEST(array_test, wrap_fail){
    // can't wrap an array with no data
    pele::Array<double> v;
    pele::Array<double> v2;
    ASSERT_THROW(v2.wrap(v), std::runtime_error);
}

TEST(array_test, wrap_self){
    pele::Array<double> v(6);
    ASSERT_EQ(v.reference_count(), 1);
    v.wrap(v);
    ASSERT_EQ(v.reference_count(), 1);
}

TEST(array_test, free){
    pele::Array<double> v(6);
    v.free();
    ASSERT_EQ (0, v.size());
    ASSERT_FALSE(v.data());
    v.free();
}



TEST(array_test, wrapping_copying){
    pele::Array<double> v(6);
    pele::Array<double> v3(6);
    ASSERT_NE(v.data(), v3.data());
    v3 = v;
    ASSERT_EQ(v.data(), v3.data());
    v3 = v.copy();
    ASSERT_NE(v.data(), v3.data());
    v3.assign(v);
    ASSERT_NE(v.data(), v3.data());
    for (int i = 0; i < 6; ++i){
        ASSERT_NEAR(v[i], v3[i], 1e-10);
    }
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

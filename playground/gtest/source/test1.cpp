#include "pele/array.h"
#include <iostream>
#include <gtest/gtest.h>

TEST(array_test, PositiveNos){
    pele::Array<double> v(6);
    std::cout << v.size() << std::endl;
    ASSERT_EQ (6, v.size());
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

#include "pele/vecn.h"
#include <iostream>
#include <stdexcept>
#include <vector>
#include <gtest/gtest.h>
#include "pele/matrix.h"


TEST(HackyMatrixTest, MatrixMultiplication_Works)
{
    pele::HackyMatrix<double> v1(2,3);
    pele::HackyMatrix<double> v2(3,2);
    for (size_t i = 0; i < 3*2; ++i){
        v1.data()[i] = i;
        v2.data()[i] = i+1;
    }
//    std::cout << v1 << std::endl;
//    std::cout << v2 << std::endl;
    pele::HackyMatrix<double> v3 = pele::hacky_mat_mul(v1, v2);
//    std::cout << v3 << std::endl;
    EXPECT_EQ(v3.shape().first, 2u);
    EXPECT_EQ(v3.shape().second, 2u);
    EXPECT_EQ(v3(0,0), 13);
    EXPECT_EQ(v3(1,0), 40);
    EXPECT_EQ(v3(0,1), 16);
    EXPECT_EQ(v3(1,1), 52);

}

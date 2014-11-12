#include "pele/vecn.h"
#include <iostream>
#include <stdexcept>
#include <vector>
#include <gtest/gtest.h>
#include "pele/matrix.h"

using pele::MatrixAdapter;
using pele::Array;

TEST(MatrixAdapterTest, ArrayWrap_Works)
{
    size_t i = 2;
    size_t j = 3;
    size_t ncol = 6;
    Array<double> v(5*ncol, 0);
    MatrixAdapter<double> m(v, ncol);
    ASSERT_FLOAT_EQ(v[i*ncol + j], 0);
    ASSERT_FLOAT_EQ(m(i,j), 0);
    m(i,j) = 1;
    ASSERT_FLOAT_EQ(v[i*ncol + j], 1);
}


TEST(MatrixAdapterTest, MatrixMultiplication_Works)
{
    pele::MatrixAdapter<double> v1(2,3);
    pele::MatrixAdapter<double> v2(3,2);
    for (size_t i = 0; i < 3*2; ++i){
        v1.data()[i] = i;
        v2.data()[i] = i+1;
    }
//    std::cout << v1 << std::endl;
//    std::cout << v2 << std::endl;
    pele::MatrixAdapter<double> v3 = pele::hacky_mat_mul(v1, v2);
//    std::cout << v3 << std::endl;
    EXPECT_EQ(v3.shape().first, 2u);
    EXPECT_EQ(v3.shape().second, 2u);
    EXPECT_EQ(v3(0,0), 13);
    EXPECT_EQ(v3(1,0), 40);
    EXPECT_EQ(v3(0,1), 16);
    EXPECT_EQ(v3(1,1), 52);

}

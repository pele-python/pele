#include <random>
#include <numeric>

#include <gtest/gtest.h>

#include "pele/combination.h"
#include "pele/array.h"

using pele::make_combination_generator;

static double const EPS = std::numeric_limits<double>::min();
#define EXPECT_NEAR_RELATIVE(A, B, T)  ASSERT_NEAR(A/(fabs(A)+fabs(B) + EPS), B/(fabs(A)+fabs(B) + EPS), T)

class CombinationTest :  public ::testing::Test {
public:

    pele::Array<size_t> data;

    virtual void SetUp()
    {
        size_t n = 5;
        data = pele::Array<size_t>(n);
        for (size_t i = 0; i < n; ++i) {
            data[i] = i;
        }
    }

    virtual void TearDown()
    {}
};

TEST_F(CombinationTest, Combination_Works)
{
    size_t r = 4;
    auto combination_generator = make_combination_generator(data.data(), data.data()+data.size(), r);
    size_t it[r];
    size_t j = 0;

    bool go_on = true;
    while(go_on){
        go_on = combination_generator(it);
        for(size_t j=0; j<r; ++j){
            std::cout<<it[j]<<" ";
        }
        std::cout<<"\n";
        ++j;
    }
    ASSERT_EQ(j, (size_t) 5);
}

TEST_F(CombinationTest, Combination_Works2)
{
    size_t r = 4;
    auto combination_generator = make_combination_generator(data.data(), data.data()+data.size(), r);
    std::ostringstream oss;
    std::ostream_iterator<size_t> it(oss, "");

    bool go_on = true;
    while(go_on){
        go_on = combination_generator(it);
    }
    ASSERT_EQ(oss.str().compare(std::string("12340234013401240123")), 0);
}



#include <stdexcept>
#include <gtest/gtest.h>
#include "test_array.hpp"
#include "test_base_potential.hpp"
#include "test_lj.hpp"
#include "test_lbfgs.hpp"
#include "test_morse.hpp"

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

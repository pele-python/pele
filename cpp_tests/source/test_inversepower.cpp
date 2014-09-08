#include <iostream>
#include <stdexcept>
#include <gtest/gtest.h>

#include "pele/array.h"
#include "pele/inversepower.h"

using pele::Array;
using pele::InversePower;

static double const EPS = std::numeric_limits<double>::min();
#define EXPECT_NEAR_RELATIVE(A, B, T)  EXPECT_NEAR(A/(fabs(A)+fabs(B) + EPS), B/(fabs(A)+fabs(B) + EPS), T)

/*
 * InversePower tests
 */

class InversePowerTest :  public ::testing::Test
{
public:
    double pow, eps, etrue;
    Array<double> x, g, gnum, radii;
    virtual void SetUp(){
    	pow = 2.5;
    	eps = 1;
        x = Array<double>(9);
        x[0] = 0.1;
        x[1] = 0.2;
        x[2] = 0.3;
        x[3] = 0.44;
        x[4] = 0.55;
        x[5] = 1.66;
        x[6] = 0.88;
        x[7] = 1.1;
        x[8] = 3.32;
        radii = Array<double>(3);
        double f = 1.;
        radii[0] = .91 * f;
        radii[1] = 1.1 * f;
        radii[2] = 1.13 * f;
        etrue = 0.023173380132354017;
        g = Array<double>(x.size());
        gnum = Array<double>(x.size());
    }
};

TEST_F(InversePowerTest, Energy_Works){
    InversePower<3> pot(pow, eps, radii);
    double e = pot.get_energy(x);
    ASSERT_NEAR(e, etrue, 1e-10);
}

TEST_F(InversePowerTest, EnergyGradient_AgreesWithNumerical){
	InversePower<3> pot(pow, eps, radii);
	double e = pot.get_energy_gradient(x, g);
	std::cout<<"energy"<<e<<std::endl;
    double ecomp = pot.get_energy(x);
    ASSERT_NEAR(e, ecomp, 1e-10);
    pot.numerical_gradient(x, gnum, 1e-6);
    for (size_t k = 0; k < 6; ++k) {
        ASSERT_NEAR(g[k], gnum[k], 1e-6);
    }
}

TEST_F(InversePowerTest, EnergyGradientHessian_AgreesWithNumerical){
	InversePower<3> pot(pow, eps, radii);
    Array<double> h(x.size()*x.size());
    Array<double> hnum(h.size());
    double e = pot.get_energy_gradient_hessian(x, g, h);
    double ecomp = pot.get_energy(x);
    pot.numerical_gradient(x, gnum);
    pot.numerical_hessian(x, hnum);
    EXPECT_NEAR(e, ecomp, 1e-10);
    for (size_t i = 0; i < g.size(); ++i) {
        ASSERT_NEAR(g[i], gnum[i], 1e-6);
    }
    for (size_t i = 0; i < h.size(); ++i) {
        ASSERT_NEAR(h[i], hnum[i], 1e-3);
    }
}

TEST_F(InversePowerTest, MetaPowFunctionsBasic_Work){
    const double op = 42.42;
    const int POW = 5;
    double true_result_direct = op * op * op * op * op;
    double true_result_std = std::pow(op, POW);
    EXPECT_DOUBLE_EQ(true_result_direct, true_result_std);
    EXPECT_DOUBLE_EQ(true_result_direct, pele::pos_int_pow<POW>(op));
    true_result_direct = double(1) / true_result_direct;
    true_result_std = std::pow(op, - POW);
    EXPECT_DOUBLE_EQ(true_result_direct, true_result_std);
    EXPECT_DOUBLE_EQ(true_result_direct, pele::neg_int_pow<- POW>(op));
    true_result_direct = std::sqrt(op * op * op * op * op);
    true_result_std = std::pow(op, 0.5 * POW);
    EXPECT_DOUBLE_EQ(true_result_direct, true_result_std);
    EXPECT_DOUBLE_EQ(true_result_direct, pele::pos_half_int_pow<POW>(op));
    true_result_direct = double(1) / std::sqrt(op * op * op * op * op);
    true_result_std = std::pow(op, - 0.5 * POW);
    EXPECT_DOUBLE_EQ(true_result_direct, true_result_std);
    EXPECT_DOUBLE_EQ(true_result_direct, pele::neg_half_int_pow<- POW>(op));
}

TEST_F(InversePowerTest, InverseIntPower_AgreesWithInversePower){
    const int pow = 4;
    pele::InversePower<3> pot(pow, eps, radii);
    pele::InverseIntPower<3, 4> pot_int(eps, radii);
    const double e = pot.get_energy(x);
    const double e_int = pot_int.get_energy(x);
    ASSERT_NEAR(e, e_int, 1e-10);
    pot.numerical_gradient(x, gnum, 1e-6);
    pele::Array<double> gnum_int(gnum.size());
    pot_int.numerical_gradient(x, gnum_int, 1e-6);
    for (size_t k = 0; k < 6; ++k) {
        ASSERT_NEAR(gnum[k], gnum_int[k], 1e-6);
    }
    pele::Array<double> h(x.size()*x.size());
    pele::Array<double> hnum(h.size());
    pot.get_energy_gradient_hessian(x, g, h);
    pele::Array<double>g_int(g.size());
    pele::Array<double>h_int(h.size());
    pot_int.get_energy_gradient_hessian(x, g_int, h_int);
    pot.numerical_gradient(x, gnum);
    pot.numerical_hessian(x, hnum);
    pele::Array<double> hnum_int(hnum.size());
    pot_int.numerical_gradient(x, gnum_int);
    pot_int.numerical_hessian(x, hnum_int);
    for (size_t i = 0; i < gnum.size(); ++i) {
        ASSERT_NEAR(gnum[i], gnum_int[i], 1e-10);
    }
    for (size_t i = 0; i < hnum.size(); ++i) {
        ASSERT_NEAR(hnum[i], hnum_int[i], 1e-10);
    }
}

TEST_F(InversePowerTest, InverseHalfIntPower_AgreesWithInversePower){
    const double pow = 2.5;
    pele::InversePower<3> pot(pow, eps, radii);
    pele::InverseHalfIntPower<3, 5> pot_int(eps, radii);
    const double e = pot.get_energy(x);
    const double e_int = pot_int.get_energy(x);
    ASSERT_NEAR(e, e_int, 1e-10);
    pot.numerical_gradient(x, gnum, 1e-6);
    pele::Array<double> gnum_int(gnum.size());
    pot_int.numerical_gradient(x, gnum_int, 1e-6);
    for (size_t k = 0; k < 6; ++k) {
        ASSERT_NEAR(gnum[k], gnum_int[k], 1e-6);
    }
    pele::Array<double> h(x.size() * x.size());
    pele::Array<double> hnum(h.size());
    pot.get_energy_gradient_hessian(x, g, h);
    pele::Array<double>g_int(g.size());
    pele::Array<double>h_int(h.size());
    pot_int.get_energy_gradient_hessian(x, g_int, h_int);
    pot.numerical_gradient(x, gnum);
    pot.numerical_hessian(x, hnum);
    pele::Array<double> hnum_int(hnum.size());
    pot_int.numerical_gradient(x, gnum_int);
    pot_int.numerical_hessian(x, hnum_int);
    for (size_t i = 0; i < gnum.size(); ++i) {
        ASSERT_NEAR(gnum[i], gnum_int[i], 1e-10);
    }
    for (size_t i = 0; i < hnum.size(); ++i) {
        ASSERT_NEAR(hnum[i], hnum_int[i], 1e-10);
    }
}

//BEGIN: TEST_F(InversePowerTest, MetaPowFunctionsLoop_Work)

template<int N>
struct perform_tests{
    /**
     *This test should work fine.
     *The comparison is to Mathematica.
     */
    static void f(const double op)
    {
        std::vector<double> true_result_pos_int_pow = {
                1.,3.141592653589793238462643,9.869604401089358618834491,
                   31.00627668029982017547632,97.40909103400243723644033,306.0196847852814532627413,
                   961.3891935753044370302194,3020.293227776792067514206,9488.531016070574007128576,
                   29809.0993334462116665094,93648.04747608302097371669,294204.0179738905971056956,
                   924269.1815233741862225792,2.903677270613283404988596e6,
                   9.122171181754353170204375e6,2.865814596938799845337882e7,
                   9.003222084293327956713077e7,2.828445635865330531542306e8,
                   8.885824030712633806702436e8,2.791563949597845565289898e9,
                   8.769956796082699474752256e9,2.755163184287328907022773e10,
                   8.655600419198134152251136e10,2.719237068936159300609803e11,
                   8.5427351991388802218689e11,2.683779414317764549009928e12,
                   8.431341691876207066643299e12,2.648784111910363022711004e13,
                   8.321400706922961226063773e13,2.614245132844608709628724e14,
                   8.212893304027495815865036e14,2.580156526864958510404037e15,
                   8.105800789910709653155358e15,2.546512421304582847058354e16,
                   8.000104715045633955294925e16,2.513307020073642986160889e17,
                   7.895786870479011810888161e17,2.480534602660760780433711e18,
                   7.792829284694322855623028e18,2.448189523147508805464168e19,
                   7.691214220515712725726519e19,2.416266209235751111302893e20,
                   7.590924172052293908738764e20
        };
        std::vector<double> true_result_neg_int_pow = {
                1.,0.3183098861837906715377675,0.1013211836423377714438795,
                   0.03225153443319948918442205,0.01026598225468433518915278,
                   0.00326776364305338547262825,0.001040161473295852296089838,
                   0.0003310936801775667643259528,0.0001053903916534936663317287,
                   0.00003354680357208869128739854,0.00001067827922686153366204078,
                   3.399001845341031027722038e-6,1.081935890528998049270004e-6,
                   3.443908901724435726384944e-7,1.096230250535248669242748e-7,
                   3.489409262791033367941994e-8,1.110713465287678744850282e-8,
                   3.535510767185247521491798e-9,1.125388029904302577577821e-9,
                   3.582221357114389661349794e-10,1.140256472468225530884818e-10,
                   3.62954907971691510158457e-11,1.155321354463173420997521e-11,
                   3.677502088448756104146607e-12,1.17058527121477605077025e-12,
                   3.726088644487970993249919e-13,1.186050852337680909412383e-13,
                   3.775317118157951264085158e-14,1.201720762188574065278781e-14,
                   3.825195990369431871049614e-15,1.217597700325186296592143e-15,
                   3.875733854081553120026146e-16,1.233684401971363536081775e-16,
                   3.926939415782225871387328e-17,1.249983638488281750262295e-17,
                   3.978821496988055085908549e-18,1.266498217851887432617611e-18,
                   4.031389035764080115217285e-19,1.283230985136645962482754e-19,
                   4.084651088263593552365809e-20,1.300184823005681168113334e-20,
                   4.138616830288303917179299e-21,1.31736265220739053348689e-21
        };
        std::vector<double> true_result_pos_half_int_pow = {
                1.,1.772453850905516027298167,3.141592653589793238462643,
                   5.568327996831707845284818,9.869604401089358618834491,17.49341832762486284626282,
                   31.00627668029982017547632,54.95719450423931590518358,97.40909103400243723644033,
                   172.6531185164235939249475,306.0196847852814532627413,542.4057687505642644109464,
                   961.3891935753044370302194,1704.017978371496937589994,3020.293227776792067514206,
                   5353.330362436825965607015,9488.531016070574007128576,16817.98333887071768121264,
                   29809.0993334462116665094,52835.25290559178884197817,93648.04747608302097371669,
                   165986.8423787659413592678,294204.0179738905971056956,521463.0446096980418896652,
                   924269.1815233741862225792,1.638224470064393998768002e6,
                   2.903677270613283404988596e6,5.146633960085332338954873e6,
                   9.122171181754353170204375e6,1.616862743971982523962448e7,
                   2.865814596938799845337882e7,5.079524118325415049041655e7,
                   9.003222084293327956713077e7,1.595779565386329556084115e8,
                   2.828445635865330531542306e8,5.013289359366406037772319e8,
                   8.885824030712633806702436e8,1.574971302170538210989138e9,
                   2.791563949597845565289898e9,4.947918272513713221492129e9,
                   8.769956796082699474752256e9,1.554434369549178203998298e10,
                   2.755163184287328907022773e10
        };
        std::vector<double> true_result_neg_half_int_pow = {
                1.,0.5641895835477562869480795,0.3183098861837906715377675,
                   0.179587122125166561689082,0.1013211836423377714438795,0.05716435640373628375718308,
                   0.03225153443319948918442205,0.01819597978064294190827957,
                   0.01026598225468433518915278,0.005791960252979011188701099,
                   0.00326776364305338547262825,0.001843638208906788476283242,
                   0.001040161473295852296089838,0.000586848268441207528999498,
                   0.0003310936801775667643259528,0.0001867996055346754036057575,
                   0.0001053903916534936663317287,0.00005946016117691952171296072,
                   0.00003354680357208869128739854,0.00001892675713669510174007699,
                   0.00001067827922686153366204078,6.024573910009665655827734e-6,
                   3.399001845341031027722038e-6,1.917681435601011418454231e-6,
                   1.081935890528998049270004e-6,6.104169595029262449374847e-7,
                   3.443908901724435726384944e-7,1.943017529040320125975439e-7,
                   1.096230250535248669242748e-7,6.184816885219344851747217e-8,
                   3.489409262791033367941994e-8,1.968688358801756393426931e-8,
                   1.110713465287678744850282e-8,6.266529674215407297042577e-9,
                   3.535510767185247521491798e-9,1.994698347366853133200588e-9,
                   1.125388029904302577577821e-9,6.34932203921338369785651e-10,
                   3.582221357114389661349794e-10,2.02105197564624585593554e-10,
                   1.140256472468225530884818e-10,6.433208243394817945475693e-11,
                   3.62954907971691510158457e-11
        };
        if (N < 16) {
            EXPECT_DOUBLE_EQ(true_result_pos_int_pow.at(N), pele::pos_int_pow<N>(op));
            EXPECT_DOUBLE_EQ(true_result_neg_int_pow.at(N), pele::neg_int_pow<- N>(op));
            EXPECT_DOUBLE_EQ(true_result_pos_half_int_pow.at(N), pele::pos_half_int_pow<N>(op));
            EXPECT_DOUBLE_EQ(true_result_neg_half_int_pow.at(N), pele::neg_half_int_pow<- N>(op));
        }
        else {
            // for large exponents, there seem to be some issues in the last digit or two
            const double prec = 1e-15;
            EXPECT_NEAR_RELATIVE(true_result_pos_int_pow.at(N), pele::pos_int_pow<N>(op), prec);
            EXPECT_NEAR_RELATIVE(true_result_neg_int_pow.at(N), pele::neg_int_pow<- N>(op), prec);
            EXPECT_NEAR_RELATIVE(true_result_pos_half_int_pow.at(N), pele::pos_half_int_pow<N>(op), prec);
            EXPECT_NEAR_RELATIVE(true_result_neg_half_int_pow.at(N), pele::neg_half_int_pow<- N>(op), prec);
        }
        return perform_tests<N - 1>::f(op);
    }
    /**
     *It seems that this test fails on some platforms.
     *The comparison is to std::pow.
     */
    static void f_std_comparison(const double op)
    {
        double true_result = std::pow(op, N);
        EXPECT_DOUBLE_EQ(true_result, pele::pos_int_pow<N>(op));
        true_result = std::pow(op, - N);
        EXPECT_DOUBLE_EQ(true_result, pele::neg_int_pow<- N>(op));
        true_result = std::pow(op, 0.5 * N);
        EXPECT_DOUBLE_EQ(true_result, pele::pos_half_int_pow<N>(op));
        true_result = std::pow(op, - 0.5 * N);
        EXPECT_DOUBLE_EQ(true_result, pele::neg_half_int_pow<- N>(op));
        return perform_tests<N - 1>::f(op);
    }
};

template<>
struct perform_tests<0>{
  static void f(const double op)
  {
      const double true_result_direct = 1;
      const double true_result_std = std::pow(op, 0);
      EXPECT_DOUBLE_EQ(true_result_direct, true_result_std);
      EXPECT_DOUBLE_EQ(true_result_direct, pele::pos_int_pow<0>(op));
  }
};

template<int N>
inline void perform_tests_zero_to(const double op)
{
    return perform_tests<N>::f(op);
}

TEST_F(InversePowerTest, MetaPowFunctionsLoop_Work){
    const double op = M_PI;
    const int POW_MAXIMUM = 42;
    perform_tests_zero_to<POW_MAXIMUM>(op);
}

//END: TEST_F(InversePowerTest, MetaPowFunctionsLoop_Work)

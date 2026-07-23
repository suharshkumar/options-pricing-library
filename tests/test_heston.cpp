#include "pricer/BlackScholes.h"
#include "pricer/Heston.h"

#include <gtest/gtest.h>

#include <cmath>

using namespace pricer;

// v0 = theta and (almost) no vol-of-vol -> volatility is effectively constant,
// so Heston must collapse to Black-Scholes with sigma = sqrt(v0). This is an
// INDEPENDENT correctness check on the characteristic-function integration.
TEST(Heston, ReducesToBlackScholesWhenVolOfVolTiny) {
    HestonParams p{100.0, 100.0, 0.03, 0.0, 1.0, 0.04, 2.0, 0.04, 1e-3, 0.0};
    const double heston = hestonPrice(p, true);
    Option bs{OptionType::Call, 100.0, 100.0, 0.03, 0.0, std::sqrt(0.04), 1.0};
    EXPECT_NEAR(heston, blackScholesPrice(bs), 2e-2);
}

TEST(Heston, PutCallParityHolds) {
    HestonParams p{100.0, 105.0, 0.02, 0.01, 1.0, 0.05, 1.5, 0.04, 0.5, -0.6};
    const double call = hestonPrice(p, true);
    const double put  = hestonPrice(p, false);
    const double rhs  = p.S * std::exp(-p.q * p.T) - p.K * std::exp(-p.r * p.T);
    EXPECT_NEAR(call - put, rhs, 1e-6);
}

TEST(Heston, AnalyticalAgreesWithMonteCarlo) {
    HestonParams p{100.0, 100.0, 0.02, 0.0, 1.0, 0.04, 1.5, 0.04, 0.5, -0.7};
    const double ana = hestonPrice(p, true);
    const double mc  = hestonMonteCarlo(p, true, 300000, 200, 123);
    EXPECT_NEAR(mc, ana, 0.10);
}

TEST(Heston, CallPriceRisesWithInitialVariance) {
    HestonParams lo{100.0, 100.0, 0.02, 0.0, 1.0, 0.02, 1.5, 0.04, 0.5, -0.7};
    HestonParams hi = lo;
    hi.v0 = 0.09;
    EXPECT_GT(hestonPrice(hi, true), hestonPrice(lo, true));
}

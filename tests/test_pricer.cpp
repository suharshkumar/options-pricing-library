#include "pricer/BlackScholes.h"
#include "pricer/Binomial.h"
#include "pricer/MonteCarlo.h"

#include <gtest/gtest.h>

#include <cmath>

using namespace pricer;

static Option atmCall() { return Option{OptionType::Call, 100, 100, 0.05, 0.0, 0.20, 1.0}; }
static Option atmPut()  { return Option{OptionType::Put,  100, 100, 0.05, 0.0, 0.20, 1.0}; }

// --- Black-Scholes ----------------------------------------------------------

TEST(BlackScholes, MatchesKnownTextbookValues) {
    EXPECT_NEAR(blackScholesPrice(atmCall()), 10.4506, 1e-3);
    EXPECT_NEAR(blackScholesPrice(atmPut()),   5.5735, 1e-3);
}

TEST(BlackScholes, PutCallParityHolds) {
    const Option c = atmCall(), p = atmPut();
    const double lhs = blackScholesPrice(c) - blackScholesPrice(p);
    const double rhs = c.S * std::exp(-c.q * c.T) - c.K * std::exp(-c.r * c.T);
    EXPECT_NEAR(lhs, rhs, 1e-10);
}

TEST(BlackScholes, GreeksHaveCorrectSignsAndRanges) {
    const Greeks gc = blackScholesGreeks(atmCall());
    const Greeks gp = blackScholesGreeks(atmPut());
    EXPECT_GT(gc.delta, 0.0); EXPECT_LT(gc.delta, 1.0);
    EXPECT_LT(gp.delta, 0.0); EXPECT_GT(gp.delta, -1.0);
    EXPECT_GT(gc.gamma, 0.0);          // gamma always positive for long options
    EXPECT_GT(gc.vega,  0.0);          // vega  always positive for long options
}

TEST(BlackScholes, DeltaMatchesFiniteDifference) {
    const Option o = atmCall();
    const double h = 1e-4;
    Option up = o, dn = o; up.S += h; dn.S -= h;
    const double fd = (blackScholesPrice(up) - blackScholesPrice(dn)) / (2.0 * h);
    EXPECT_NEAR(blackScholesGreeks(o).delta, fd, 1e-6);
}

// --- Binomial ---------------------------------------------------------------

TEST(Binomial, ConvergesToBlackScholes) {
    const Option o = atmCall();
    EXPECT_NEAR(binomialPrice(o, 2000), blackScholesPrice(o), 0.01);
}

TEST(Binomial, AmericanPutIsWorthAtLeastEuropean) {
    const Option o = atmPut();
    const double european = binomialPrice(o, 500, false);
    const double american = binomialPrice(o, 500, true);
    EXPECT_GE(american, european);     // early-exercise right can only add value
}

// --- Monte Carlo ------------------------------------------------------------

TEST(MonteCarlo, AgreesWithBlackScholes) {
    const Option o = atmCall();
    const MCResult r = monteCarlo(o, 500'000, 12345, true);
    EXPECT_NEAR(r.price, blackScholesPrice(o), 0.05);
}

TEST(MonteCarlo, AntitheticReducesStandardError) {
    const Option o = atmCall();
    const MCResult plain = monteCarlo(o, 400'000, 999, false);  // 400k evaluations
    const MCResult anti  = monteCarlo(o, 200'000, 999, true);   // 200k pairs = 400k evals
    EXPECT_LT(anti.stdError, plain.stdError);
}

TEST(MonteCarlo, ParallelAgreesWithBlackScholes) {
    const Option o = atmCall();
    const MCResult par = monteCarloParallel(o, 500'000, 2024, 4);
    EXPECT_NEAR(par.price, blackScholesPrice(o), 0.10);
}

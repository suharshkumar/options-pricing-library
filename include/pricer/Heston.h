#pragma once

#include <cstddef>

namespace pricer {

// Heston (1993) stochastic-volatility model. Unlike Black-Scholes, volatility is
// itself random and mean-reverting, which lets the model reproduce the implied-
// volatility SMILE/SKEW that real option markets show and Black-Scholes cannot.
//
//   dS = (r - q) S dt + sqrt(v) S dW1
//   dv = kappa (theta - v) dt + sigma sqrt(v) dW2,   corr(dW1, dW2) = rho
struct HestonParams {
    double S;      // spot
    double K;      // strike
    double r;      // risk-free rate
    double q;      // dividend yield
    double T;      // maturity (years)
    double v0;     // initial variance
    double kappa;  // mean-reversion speed of variance
    double theta;  // long-run variance
    double sigma;  // volatility of variance ("vol of vol")
    double rho;    // correlation between spot and variance shocks
};

// Semi-analytical price via characteristic-function (Fourier) integration —
// the numerically stable "little trap" formulation, integrated with Simpson's rule.
double hestonPrice(const HestonParams& p, bool isCall = true);

// Monte Carlo price using a full-truncation Euler scheme for the variance path.
// Cross-validates the analytical price and extends to path-dependent payoffs.
double hestonMonteCarlo(const HestonParams& p, bool isCall,
                        std::size_t paths, std::size_t steps, unsigned seed = 42);

} // namespace pricer

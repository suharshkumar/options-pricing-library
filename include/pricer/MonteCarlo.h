#pragma once

#include "pricer/Option.h"

#include <cstddef>

namespace pricer {

struct MCResult {
    double      price;      // discounted expected payoff
    double      stdError;   // standard error of the estimate (uncertainty)
    std::size_t samples;    // number of estimator samples used
};

// Single-threaded Monte Carlo under geometric Brownian motion.
// If `antithetic` is true, each draw Z is paired with -Z (variance reduction).
MCResult monteCarlo(const Option& o, std::size_t numSamples,
                    unsigned seed = 42, bool antithetic = true);

// Multithreaded plain Monte Carlo. Each thread gets an independent RNG stream,
// then partial sums are combined. `numThreads == 0` uses hardware concurrency.
MCResult monteCarloParallel(const Option& o, std::size_t numSamples,
                            unsigned seed = 42, unsigned numThreads = 0);

} // namespace pricer

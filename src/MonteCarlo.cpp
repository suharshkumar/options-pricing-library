#include "pricer/MonteCarlo.h"

#include <algorithm>
#include <cmath>
#include <random>
#include <thread>
#include <vector>

namespace pricer {

namespace {
double payoff(const Option& o, double terminalSpot) {
    return (o.type == OptionType::Call) ? std::max(terminalSpot - o.K, 0.0)
                                        : std::max(o.K - terminalSpot, 0.0);
}

// Turn running sum / sum-of-squares into (mean, standard error).
MCResult finalize(double sum, double sumSq, std::size_t n) {
    const double dn   = static_cast<double>(n);
    const double mean = sum / dn;
    const double var  = (sumSq - sum * sum / dn) / (dn - 1.0);  // sample variance
    return MCResult{mean, std::sqrt(var / dn), n};
}
} // namespace

MCResult monteCarlo(const Option& o, std::size_t numSamples, unsigned seed, bool antithetic) {
    // Under GBM, S_T = S * exp(drift + vol * Z), Z ~ N(0,1).
    const double drift = (o.r - o.q - 0.5 * o.sigma * o.sigma) * o.T;
    const double vol   = o.sigma * std::sqrt(o.T);
    const double disc  = std::exp(-o.r * o.T);

    std::mt19937_64                  rng(seed);
    std::normal_distribution<double> Z(0.0, 1.0);

    double sum = 0.0, sumSq = 0.0;
    for (std::size_t i = 0; i < numSamples; ++i) {
        double sample;
        if (antithetic) {
            const double z  = Z(rng);
            const double up = disc * payoff(o, o.S * std::exp(drift + vol * z));
            const double dn = disc * payoff(o, o.S * std::exp(drift - vol * z));
            sample = 0.5 * (up + dn);   // averaging +Z and -Z cancels much of the noise
        } else {
            sample = disc * payoff(o, o.S * std::exp(drift + vol * Z(rng)));
        }
        sum   += sample;
        sumSq += sample * sample;
    }
    return finalize(sum, sumSq, numSamples);
}

MCResult monteCarloParallel(const Option& o, std::size_t numSamples, unsigned seed, unsigned numThreads) {
    if (numThreads == 0)
        numThreads = std::max(1u, std::thread::hardware_concurrency());

    const double drift = (o.r - o.q - 0.5 * o.sigma * o.sigma) * o.T;
    const double vol   = o.sigma * std::sqrt(o.T);
    const double disc  = std::exp(-o.r * o.T);

    struct Partial { double sum = 0.0, sumSq = 0.0; std::size_t n = 0; };
    std::vector<Partial>     partials(numThreads);
    std::vector<std::thread> threads;
    threads.reserve(numThreads);

    for (unsigned t = 0; t < numThreads; ++t) {
        const std::size_t begin = numSamples * t / numThreads;
        const std::size_t end   = numSamples * (t + 1) / numThreads;
        threads.emplace_back([&, t, begin, end]() {
            std::mt19937_64                  rng(seed + t);   // independent stream per thread
            std::normal_distribution<double> Z(0.0, 1.0);
            double sum = 0.0, sumSq = 0.0;
            for (std::size_t i = begin; i < end; ++i) {
                const double s = disc * payoff(o, o.S * std::exp(drift + vol * Z(rng)));
                sum += s;
                sumSq += s * s;
            }
            partials[t] = Partial{sum, sumSq, end - begin};
        });
    }
    for (auto& th : threads) th.join();

    double sum = 0.0, sumSq = 0.0;
    std::size_t n = 0;
    for (const auto& p : partials) { sum += p.sum; sumSq += p.sumSq; n += p.n; }
    return finalize(sum, sumSq, n);
}

} // namespace pricer

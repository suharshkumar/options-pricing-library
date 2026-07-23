#include "pricer/BlackScholes.h"
#include "pricer/Binomial.h"
#include "pricer/MonteCarlo.h"

#include <chrono>
#include <cstdio>

using namespace pricer;
using Clock = std::chrono::steady_clock;

int main() {
    const Option call{OptionType::Call, 100.0, 100.0, 0.05, 0.0, 0.20, 1.0};
    const double bs = blackScholesPrice(call);

    std::printf("European CALL  S=100  K=100  r=5%%  q=0  sigma=20%%  T=1y\n\n");
    std::printf("Black-Scholes (exact closed form):   %.4f\n", bs);

    std::printf("\nBinomial tree -> converges to Black-Scholes as steps grow:\n");
    for (int steps : {10, 50, 200, 1000}) {
        const double bp = binomialPrice(call, steps);
        std::printf("   steps = %5d   price = %.4f   (error %+.4f)\n", steps, bp, bp - bs);
    }

    std::printf("\nMonte Carlo -> error shrinks ~ 1/sqrt(N):\n");
    for (std::size_t n : {1'000ul, 10'000ul, 100'000ul, 1'000'000ul}) {
        const MCResult r = monteCarlo(call, n, 42, true);
        std::printf("   N = %8zu   price = %.4f   (stderr %.4f)\n", n, r.price, r.stdError);
    }

    // Variance reduction: same number of payoff evaluations, lower error.
    const MCResult plain = monteCarlo(call, 2'000'000, 7, false);  // 2M evals
    const MCResult anti  = monteCarlo(call, 1'000'000, 7, true);   // 1M pairs = 2M evals
    std::printf("\nVariance reduction (same 2M evaluations):\n");
    std::printf("   plain        stderr = %.5f\n", plain.stdError);
    std::printf("   antithetic   stderr = %.5f   (%.1fx tighter)\n",
                anti.stdError, plain.stdError / anti.stdError);

    // Parallel speedup.
    const auto t0 = Clock::now();
    const MCResult serial = monteCarlo(call, 10'000'000, 1, false);
    const auto t1 = Clock::now();
    const MCResult par    = monteCarloParallel(call, 10'000'000, 1, 0);
    const auto t2 = Clock::now();
    const double msSerial = std::chrono::duration<double, std::milli>(t1 - t0).count();
    const double msPar    = std::chrono::duration<double, std::milli>(t2 - t1).count();
    std::printf("\nParallel Monte Carlo (10M paths):\n");
    std::printf("   1 thread    price = %.4f   %.0f ms\n", serial.price, msSerial);
    std::printf("   N threads   price = %.4f   %.0f ms   (%.1fx faster)\n",
                par.price, msPar, msSerial / msPar);

    const Greeks g = blackScholesGreeks(call);
    std::printf("\nGreeks:  delta=%.4f  gamma=%.4f  vega=%.4f  theta=%.4f  rho=%.4f\n",
                g.delta, g.gamma, g.vega, g.theta, g.rho);
    return 0;
}

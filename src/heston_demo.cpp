// Heston demo: cross-validate the analytical price against Monte Carlo, then show
// the payoff Heston delivers that Black-Scholes cannot — a volatility SMILE. We
// price options across strikes with Heston, then back out the Black-Scholes
// implied vol of each; a flat line would mean "no smile", a curve means Heston
// reproduces the skew real markets show.

#include "pricer/BlackScholes.h"
#include "pricer/Heston.h"

#include <cmath>
#include <cstdio>
#include <initializer_list>

using namespace pricer;

// Invert Black-Scholes for the implied volatility that reproduces `price`.
static double impliedVol(double price, double S, double K, double r, double q, double T, bool isCall) {
    double lo = 1e-4, hi = 3.0;
    for (int i = 0; i < 100; ++i) {
        const double mid = 0.5 * (lo + hi);
        Option o{isCall ? OptionType::Call : OptionType::Put, S, K, r, q, mid, T};
        if (blackScholesPrice(o) > price) hi = mid;
        else                              lo = mid;
    }
    return 0.5 * (lo + hi);
}

int main() {
    HestonParams p{100.0, 100.0, 0.02, 0.0, 1.0,
                   /*v0*/ 0.04, /*kappa*/ 1.5, /*theta*/ 0.04, /*sigma*/ 0.5, /*rho*/ -0.7};

    const double ana = hestonPrice(p, true);
    const double mc  = hestonMonteCarlo(p, true, 400000, 200, 7);
    std::printf("ATM call  (S=K=100, T=1, v0=theta=0.04, vol-of-vol=0.5, rho=-0.7)\n");
    std::printf("   Heston analytical : %.4f\n", ana);
    std::printf("   Heston Monte Carlo: %.4f   (agree -> both correct)\n\n", mc);

    std::printf("Implied-volatility smile Heston generates (Black-Scholes can't):\n");
    std::printf("   strike |  Heston price |  BS implied vol\n");
    for (double K : {80.0, 90.0, 100.0, 110.0, 120.0}) {
        HestonParams pk = p;
        pk.K = K;
        const bool isCall = K >= p.S;               // price the OTM option each side
        const double price = hestonPrice(pk, isCall);
        const double iv = impliedVol(price, pk.S, K, pk.r, pk.q, pk.T, isCall);
        std::printf("   %6.0f | %13.4f | %8.2f%%\n", K, price, iv * 100.0);
    }
    std::printf("\nThe implied vol changes with strike (a skew) — the market smile that a\n"
                "single Black-Scholes sigma cannot produce. That's why Heston exists.\n");
    return 0;
}

#include "pricer/BlackScholes.h"

#include <cmath>

namespace pricer {

namespace {
constexpr double kInvSqrt2   = 0.7071067811865476;  // 1/sqrt(2)
constexpr double kInvSqrt2Pi = 0.3989422804014327;  // 1/sqrt(2*pi)

// d1 and d2 from the Black-Scholes formula.
void computeD(const Option& o, double& d1, double& d2) {
    const double volSqrtT = o.sigma * std::sqrt(o.T);
    d1 = (std::log(o.S / o.K) + (o.r - o.q + 0.5 * o.sigma * o.sigma) * o.T) / volSqrtT;
    d2 = d1 - volSqrtT;
}
} // namespace

double normalCdf(double x) { return 0.5 * std::erfc(-x * kInvSqrt2); }
double normalPdf(double x) { return kInvSqrt2Pi * std::exp(-0.5 * x * x); }

double blackScholesPrice(const Option& o) {
    double d1, d2;
    computeD(o, d1, d2);
    const double dfR = std::exp(-o.r * o.T);   // risk-free discount factor
    const double dfQ = std::exp(-o.q * o.T);   // dividend discount factor

    if (o.type == OptionType::Call)
        return o.S * dfQ * normalCdf(d1) - o.K * dfR * normalCdf(d2);
    else
        return o.K * dfR * normalCdf(-d2) - o.S * dfQ * normalCdf(-d1);
}

Greeks blackScholesGreeks(const Option& o) {
    double d1, d2;
    computeD(o, d1, d2);
    const double sqrtT = std::sqrt(o.T);
    const double dfR   = std::exp(-o.r * o.T);
    const double dfQ   = std::exp(-o.q * o.T);
    const double pdf   = normalPdf(d1);

    Greeks g{};
    g.gamma = dfQ * pdf / (o.S * o.sigma * sqrtT);   // same for call and put
    g.vega  = o.S * dfQ * pdf * sqrtT;               // same for call and put

    if (o.type == OptionType::Call) {
        g.delta = dfQ * normalCdf(d1);
        g.theta = -(o.S * dfQ * pdf * o.sigma) / (2.0 * sqrtT)
                  - o.r * o.K * dfR * normalCdf(d2)
                  + o.q * o.S * dfQ * normalCdf(d1);
        g.rho   = o.K * o.T * dfR * normalCdf(d2);
    } else {
        g.delta = -dfQ * normalCdf(-d1);
        g.theta = -(o.S * dfQ * pdf * o.sigma) / (2.0 * sqrtT)
                  + o.r * o.K * dfR * normalCdf(-d2)
                  - o.q * o.S * dfQ * normalCdf(-d1);
        g.rho   = -o.K * o.T * dfR * normalCdf(-d2);
    }
    return g;
}

} // namespace pricer

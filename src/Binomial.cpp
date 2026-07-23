#include "pricer/Binomial.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace pricer {

namespace {
double intrinsic(const Option& o, double spot) {
    return (o.type == OptionType::Call) ? std::max(spot - o.K, 0.0)
                                        : std::max(o.K - spot, 0.0);
}
} // namespace

double binomialPrice(const Option& o, int steps, bool american) {
    const double dt   = o.T / steps;
    const double u    = std::exp(o.sigma * std::sqrt(dt));   // up factor
    const double d    = 1.0 / u;                             // down factor
    const double disc = std::exp(-o.r * dt);                 // per-step discount
    const double p    = (std::exp((o.r - o.q) * dt) - d) / (u - d);  // risk-neutral prob

    // Option value at each terminal node (i = number of down moves).
    std::vector<double> value(steps + 1);
    for (int i = 0; i <= steps; ++i) {
        const double spot = o.S * std::pow(u, steps - i) * std::pow(d, i);
        value[i] = intrinsic(o, spot);
    }

    // Roll the tree backwards, discounting the risk-neutral expectation.
    for (int step = steps - 1; step >= 0; --step) {
        for (int i = 0; i <= step; ++i) {
            value[i] = disc * (p * value[i] + (1.0 - p) * value[i + 1]);
            if (american) {
                const double spot = o.S * std::pow(u, step - i) * std::pow(d, i);
                value[i] = std::max(value[i], intrinsic(o, spot));  // early exercise
            }
        }
    }
    return value[0];
}

} // namespace pricer

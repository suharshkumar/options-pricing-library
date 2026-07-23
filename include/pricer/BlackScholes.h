#pragma once

#include "pricer/Option.h"

namespace pricer {

// The five first-order sensitivities a trader watches.
struct Greeks {
    double delta;   // d(price)/d(spot)
    double gamma;   // d(delta)/d(spot)
    double vega;    // d(price)/d(vol), per 1.00 change in sigma
    double theta;   // d(price)/d(time), per year (usually negative)
    double rho;     // d(price)/d(rate), per 1.00 change in r
};

// Closed-form Black-Scholes price and Greeks for a European option.
double blackScholesPrice(const Option& o);
Greeks blackScholesGreeks(const Option& o);

// Standard normal CDF and PDF (exposed — the other pricers/tests reuse them).
double normalCdf(double x);
double normalPdf(double x);

} // namespace pricer

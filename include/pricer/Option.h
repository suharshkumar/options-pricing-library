#pragma once

namespace pricer {

enum class OptionType { Call, Put };

// A vanilla option contract plus the market data needed to price it.
// Rates and yields are continuously compounded and annualised; T is in years.
struct Option {
    OptionType type;
    double S;      // spot price of the underlying
    double K;      // strike price
    double r;      // risk-free rate
    double q;      // continuous dividend yield
    double sigma;  // volatility
    double T;      // time to maturity (years)
};

} // namespace pricer

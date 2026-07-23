#include "pricer/Heston.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <random>

namespace pricer {

namespace {
using cd = std::complex<double>;
constexpr double kPi = 3.14159265358979323846;

// Heston characteristic function of log-spot for j = 1, 2 ("little trap" form).
cd charFunc(double phi, int j, const HestonParams& p) {
    const double kappa = p.kappa, theta = p.theta, sig = p.sigma, rho = p.rho;
    const double x = std::log(p.S);
    const double a = kappa * theta;
    const double u = (j == 1) ? 0.5 : -0.5;
    const double b = (j == 1) ? (kappa - rho * sig) : kappa;

    const cd I(0.0, 1.0);
    const cd rspi = rho * sig * phi * I;
    const cd d = std::sqrt((rspi - b) * (rspi - b) - sig * sig * (2.0 * u * phi * I - phi * phi));
    const cd g = (b - rspi - d) / (b - rspi + d);          // reciprocal-g "trap"
    const cd edt = std::exp(-d * p.T);

    const cd D = ((b - rspi - d) / (sig * sig)) * ((1.0 - edt) / (1.0 - g * edt));
    const cd C = (p.r - p.q) * phi * I * p.T
               + (a / (sig * sig)) * ((b - rspi - d) * p.T - 2.0 * std::log((1.0 - g * edt) / (1.0 - g)));
    return std::exp(C + D * p.v0 + I * phi * x);
}

// P_j = 1/2 + (1/pi) * integral_0^inf Re[ exp(-i phi ln K) f_j(phi) / (i phi) ] dphi
double probability(int j, const HestonParams& p) {
    const cd I(0.0, 1.0);
    const double lnK = std::log(p.K);
    auto integrand = [&](double phi) {
        return (std::exp(-I * phi * lnK) * charFunc(phi, j, p) / (I * phi)).real();
    };
    // Simpson's rule over a truncated range; the integrand decays quickly.
    const double lo = 1e-6, hi = 200.0;
    const int    N = 4000;
    const double h = (hi - lo) / N;
    double s = integrand(lo) + integrand(hi);
    for (int k = 1; k < N; ++k)
        s += (k % 2 ? 4.0 : 2.0) * integrand(lo + k * h);
    return 0.5 + (s * h / 3.0) / kPi;
}
} // namespace

double hestonPrice(const HestonParams& p, bool isCall) {
    const double P1 = probability(1, p);
    const double P2 = probability(2, p);
    const double call = p.S * std::exp(-p.q * p.T) * P1 - p.K * std::exp(-p.r * p.T) * P2;
    if (isCall) return call;
    return call - p.S * std::exp(-p.q * p.T) + p.K * std::exp(-p.r * p.T);   // put-call parity
}

double hestonMonteCarlo(const HestonParams& p, bool isCall,
                        std::size_t paths, std::size_t steps, unsigned seed) {
    const double dt = p.T / steps;
    const double sqrtDt = std::sqrt(dt);
    const double disc = std::exp(-p.r * p.T);
    const double corr = std::sqrt(std::max(0.0, 1.0 - p.rho * p.rho));

    std::mt19937_64 rng(seed);
    std::normal_distribution<double> Z(0.0, 1.0);

    double sum = 0.0;
    for (std::size_t i = 0; i < paths; ++i) {
        double x = std::log(p.S);   // log-spot
        double v = p.v0;
        for (std::size_t t = 0; t < steps; ++t) {
            const double z1 = Z(rng);
            const double z2 = p.rho * z1 + corr * Z(rng);   // correlate the shocks
            const double vPos = std::max(v, 0.0);           // full-truncation scheme
            x += (p.r - p.q - 0.5 * vPos) * dt + std::sqrt(vPos) * sqrtDt * z2;
            v += p.kappa * (p.theta - vPos) * dt + p.sigma * std::sqrt(vPos) * sqrtDt * z1;
        }
        const double ST = std::exp(x);
        sum += isCall ? std::max(ST - p.K, 0.0) : std::max(p.K - ST, 0.0);
    }
    return disc * sum / static_cast<double>(paths);
}

} // namespace pricer

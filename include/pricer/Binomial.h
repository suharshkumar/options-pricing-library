#pragma once

#include "pricer/Option.h"

namespace pricer {

// Cox-Ross-Rubinstein binomial tree.
//   steps    — number of time steps (more steps -> converges to Black-Scholes)
//   american — if true, allow early exercise at every node
double binomialPrice(const Option& o, int steps, bool american = false);

} // namespace pricer

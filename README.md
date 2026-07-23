# Options Pricing Library (C++17)

Prices European and American options **three independent ways** — closed-form
**Black-Scholes**, a **Cox-Ross-Rubinstein binomial tree**, and **Monte Carlo**
simulation — plus the full set of **Greeks**. The three methods cross-validate
each other, and the Monte Carlo engine demonstrates variance reduction and
multithreading. It also includes the **Heston stochastic-volatility model**, which
reproduces the implied-volatility smile that Black-Scholes cannot.

## Methods

| Method | Handles | Notes |
|---|---|---|
| Black-Scholes | European | exact closed form + analytic Greeks |
| Binomial (CRR) | European **and American** | early exercise via backward induction |
| Monte Carlo | European | GBM paths, antithetic variates, multithreaded |
| **Heston** | European (stochastic vol) | characteristic-function integration + MC; reproduces the vol smile |

## Results (`./build/demo`, S=K=100, r=5%, σ=20%, T=1y)

All three agree on the true price of **10.4506**:

```
Black-Scholes (exact) : 10.4506

Binomial convergence  : 10.2534 (10 steps) -> 10.4486 (1000 steps)
Monte Carlo (1M paths): 10.4518  ± 0.0074

Antithetic variates   : 1.4x tighter standard error for the same compute
Parallel Monte Carlo  : 6.4x faster (10M paths across all cores)

Greeks : delta 0.637  gamma 0.019  vega 37.5  theta -6.41  rho 53.2
```

Monte Carlo standard error falls from **0.239** (1k paths) to **0.0074** (1M paths)
— a ~32× drop for 1000× the paths, exactly the **1/√N** convergence rate.

## Stochastic volatility (Heston)

Black-Scholes assumes constant volatility, so it prices every strike off one sigma
and **cannot** produce the implied-vol smile real markets show. The
[Heston model](include/pricer/Heston.h) makes volatility itself random and
mean-reverting. Priced two ways — semi-analytical characteristic-function
integration and Monte Carlo — which agree:

```
ATM call (v0=theta=0.04, vol-of-vol=0.5, rho=-0.7):
   Heston analytical : 8.1950
   Heston Monte Carlo: 8.1818

Implied-vol smile it generates:   strike 80 -> 23.7%   ...   strike 120 -> 14.2%
```

The implied vol varying with strike (a skew) is the whole point — and a correctness
test confirms Heston collapses to Black-Scholes when vol-of-vol -> 0. Run `./build/heston_demo`.

## Why these choices

- **Prices via `std::erfc`** for the normal CDF — accurate and branch-free.
- **Greeks validated against finite differences** in the tests, not just coded from formulas.
- **Antithetic variates**: pairing each random draw `Z` with `-Z` cancels much of
  the sampling noise for monotonic payoffs — free variance reduction.
- **Per-thread RNG streams** (`mt19937_64` seeded per thread) so parallel Monte
  Carlo is correct, not just fast.

## Build & run

```bash
cmake -S . -B build
cmake --build build

./build/demo            # convergence + variance-reduction + parallel demo
./build/heston_demo     # Heston vs Monte Carlo + the implied-vol smile
ctest --test-dir build  # 13 unit tests (parity, convergence, Greeks, variance, Heston)
```

Requires CMake ≥ 3.14 and a C++17 compiler. GoogleTest is fetched automatically.

## Project layout

```
include/pricer/   Option.h, BlackScholes.h, Binomial.h, MonteCarlo.h, Heston.h
src/              BlackScholes.cpp, Binomial.cpp, MonteCarlo.cpp, Heston.cpp,
                  main.cpp, heston_demo.cpp
tests/            test_pricer.cpp, test_heston.cpp   (GoogleTest)
```

## Validation the tests enforce

- Black-Scholes matches known textbook values and **put-call parity** to 1e-10
- Binomial **converges** to Black-Scholes; **American put ≥ European put**
- Monte Carlo agrees with Black-Scholes within error; **antithetic lowers std error**
- Delta matches a central finite-difference estimate

## Roadmap

- [ ] Implied-volatility solver (Newton-Raphson on vega)
- [ ] Exotic payoffs (Asian, barrier) via the Monte Carlo engine
- [ ] Control-variate variance reduction
- [ ] Longstaff-Schwartz for American Monte Carlo

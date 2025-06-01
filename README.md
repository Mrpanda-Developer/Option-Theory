# Option Trading Theory in C++

## My Journey into Finance and Programming

After completing a statistics class, I became fascinated with quantitative finance and option trading. The mathematical models, probability theory, and risk management aspects particularly intrigued me. I realized that to truly understand and work in this field, I needed to combine my statistical knowledge with programming skills.

Since I'm new to both finance and programming, I decided to:
1. Learn the fundamentals of option pricing theory
2. Teach myself C++ (as it's widely used in high-performance finance applications)
3. Implement the models from scratch to ensure deep understanding

This project represents my learning journey - implementing option pricing models and trading concepts using pure C++ without any external libraries.

## Project Overview

This C++ implementation covers fundamental option trading concepts:

- **Pricing Models**:
  - Black-Scholes (analytical)
  - Binomial Tree (discrete)
  - Monte Carlo (simulation)

- **Risk Measures**:
  - Greeks (Delta, Gamma, Theta, Vega, Rho)
  - Implied volatility calculation

- **Trading Strategies**:
  - Payoff diagrams for calls, puts, spreads, straddles, butterflies

## Why C++?

I chose C++ because:
- It's the industry standard for high-performance financial systems
- Manual memory management forces deeper understanding
- No "magic" - everything is implemented from first principles
- Excellent for numerical computations

## How to Use

1. Compile with any C++ compiler:
   ```bash
   g++ options.cpp -o options -std=c++11

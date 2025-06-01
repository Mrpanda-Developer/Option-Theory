#include <iostream>
#include <cmath>
#include <random>
#include <iomanip>

using namespace std;

// Function to price European option using Monte Carlo
double monteCarloOptionPrice(
    char optionType,    // 'C' for Call, 'P' for Put
    double S,          // Current stock price
    double K,          // Strike price
    double T,          // Time to maturity (years)
    double r,          // Risk-free rate
    double sigma,      // Volatility
    float numSimulations // Number of simulations
) {
    // Initialize random number generator
    random_device rd;
    mt19937 generator(rd());
    normal_distribution<double> distribution(0.0, 1.0);
    
    double sumPayoffs = 0.0;
    
    // Run simulations
    for (int i = 0; i < numSimulations; i++) {
        // Generate random number from standard normal distribution
        double z = distribution(generator);
        
        // Calculate stock price at expiration using Geometric Brownian Motion
        // This is the core of the financial model
        double ST = S * exp((r - 0.5 * sigma * sigma) * T + sigma * sqrt(T) * z);
        
        // Calculate payoff for this simulation
        double payoff;
        if (optionType == 'C') {
            payoff = max(ST - K, 0.0); // Call option payoff
        } else {
            payoff = max(K - ST, 0.0); // Put option payoff
        }
        
        sumPayoffs += payoff;
    }
    
    // Average the payoffs and discount back to present value
    double optionPrice = exp(-r * T) * (sumPayoffs / numSimulations);
    
    return optionPrice;
}

int main() {
    // Example parameters
    char optionType = 'C'; // Call option
    double S = 100.0;     // Stock price $100
    double K = 105.0;     // Strike price $105
    double T = 2.0;       // 1 year to expiration
    double r = 0.05;      // 5% risk-free rate
    double sigma = 0.2;   // 20% volatility
    float numSimulations = 99999;
    
    double price = monteCarloOptionPrice(optionType, S, K, T, r, sigma, numSimulations);
    
    cout << fixed << setprecision(2);
    cout << "European " << (optionType == 'C' ? "Call" : "Put") << " Option Price: $" << price << endl;
    
    return 0;
}

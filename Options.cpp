#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <limits>
#include <algorithm>

using namespace std;

// Constants
const double PI = 3.14159265358979323846;
const double E = 2.71828182845904523536;

// Normal cumulative distribution function (CDF)
double normCDF(double x) {
    // Abramowitz & Stegun approximation
    double a1 = 0.254829592;
    double a2 = -0.284496736;
    double a3 = 1.421413741;
    double a4 = -1.453152027;
    double a5 = 1.061405429;
    double p = 0.3275911;
    
    int sign = (x < 0) ? -1 : 1;
    x = fabs(x) / sqrt(2.0);
    
    double t = 1.0 / (1.0 + p * x);
    double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x * x);
    
    return 0.5 * (1.0 + sign * y);
}

// Black-Scholes option pricing model
double blackScholes(char optionType, double S, double K, double T, double r, double sigma) {
    if (T <= 0 || sigma <= 0 || S <= 0) {
        return (optionType == 'C') ? max(S - K, 0.0) : max(K - S, 0.0);
    }
    
    double d1 = (log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    double d2 = d1 - sigma * sqrt(T);
    
    if (optionType == 'C') {
        return S * normCDF(d1) - K * exp(-r * T) * normCDF(d2);
    } else if (optionType == 'P') {
        return K * exp(-r * T) * normCDF(-d2) - S * normCDF(-d1);
    }
    
    return 0.0;
}

// Binomial option pricing model
double binomialOptionPrice(char optionType, double S, double K, double T, double r, double sigma, int steps) {
    if (steps < 1) return blackScholes(optionType, S, K, T, r, sigma);
    
    double dt = T / steps;
    double u = exp(sigma * sqrt(dt));
    double d = 1.0 / u;
    double p = (exp(r * dt) - d) / (u - d);
    
    vector<double> prices(steps + 1);
    
    // Calculate terminal prices
    for (int i = 0; i <= steps; i++) {
        prices[i] = S * pow(u, steps - i) * pow(d, i);
    }
    
    // Calculate terminal option values
    vector<double> values(steps + 1);
    for (int i = 0; i <= steps; i++) {
        if (optionType == 'C') {
            values[i] = max(prices[i] - K, 0.0);
        } else {
            values[i] = max(K - prices[i], 0.0);
        }
    }
    
    // Backward induction
    for (int step = steps - 1; step >= 0; step--) {
        for (int i = 0; i <= step; i++) {
            values[i] = (p * values[i] + (1 - p) * values[i + 1]) * exp(-r * dt);
        }
    }
    
    return values[0];
}

// Implied volatility calculation using Newton-Raphson method
double impliedVolatility(char optionType, double S, double K, double T, double r, double marketPrice, double epsilon = 1e-6, int maxIter = 100) {
    if (marketPrice <= 0) return 0.0;
    
    double sigma = 0.5; // Initial guess
    double diff;
    int iter = 0;
    
    do {
        double price = blackScholes(optionType, S, K, T, r, sigma);
        double vega = S * sqrt(T) * (1.0 / sqrt(2 * PI)) * exp(-0.5 * pow((log(S/K) + (r + 0.5*sigma*sigma)*T) / (sigma*sqrt(T)), 2));
        
        if (vega < epsilon) break; // Avoid division by zero
        
        diff = (price - marketPrice) / vega;
        sigma -= diff;
        iter++;
        
    } while (fabs(diff) > epsilon && iter < maxIter);
    
    return sigma;
}

// Greeks calculation
struct Greeks {
    double delta;
    double gamma;
    double theta;
    double vega;
    double rho;
};

Greeks calculateGreeks(char optionType, double S, double K, double T, double r, double sigma) {
    Greeks greeks;
    
    if (T <= 0 || sigma <= 0) {
        greeks.delta = (optionType == 'C') ? ((S > K) ? 1.0 : 0.0) : ((S < K) ? -1.0 : 0.0);
        greeks.gamma = 0.0;
        greeks.theta = 0.0;
        greeks.vega = 0.0;
        greeks.rho = 0.0;
        return greeks;
    }
    
    double d1 = (log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    double d2 = d1 - sigma * sqrt(T);
    
    double normPDF = (1.0 / sqrt(2 * PI)) * exp(-0.5 * d1 * d1);
    
    if (optionType == 'C') {
        greeks.delta = normCDF(d1);
        greeks.theta = - (S * sigma * normPDF) / (2 * sqrt(T)) - r * K * exp(-r * T) * normCDF(d2);
        greeks.rho = K * T * exp(-r * T) * normCDF(d2);
    } else { // Put
        greeks.delta = normCDF(d1) - 1;
        greeks.theta = - (S * sigma * normPDF) / (2 * sqrt(T)) + r * K * exp(-r * T) * normCDF(-d2);
        greeks.rho = -K * T * exp(-r * T) * normCDF(-d2);
    }
    
    greeks.gamma = normPDF / (S * sigma * sqrt(T));
    greeks.vega = S * sqrt(T) * normPDF;
    
    return greeks;
}

// Option strategy payoff diagrams
void plotPayoffDiagram(char strategy, double S_min, double S_max, double K1, double K2 = 0, double K3 = 0, double premium = 0) {
    const int width = 60;
    const int height = 20;
    
    vector<vector<char>> grid(height, vector<char>(width, ' '));
    
    // Draw axes
    for (int i = 0; i < width; i++) {
        grid[height/2][i] = '-';
    }
    for (int i = 0; i < height; i++) {
        grid[i][width/2] = '|';
    }
    grid[height/2][width/2] = '+';
    
    // Calculate scale
    double S_range = S_max - S_min;
    double x_scale = width / S_range;
    double y_max = max(fabs(K1 - S_max), fabs(K1 - S_min));
    if (strategy == 'S' || strategy == 'V') { // Spread or butterfly
        y_max = max(fabs(K1 - S_max) - fabs(K2 - S_max), 
                   fabs(K1 - S_min) - fabs(K2 - S_min));
    }
    y_max = max(y_max, 0.2 * S_max); // Ensure some visibility
    double y_scale = (height/2 - 1) / y_max;
    
    // Plot payoff
    for (int i = 0; i < width; i++) {
        double S = S_min + (i / x_scale);
        double payoff = 0;
        
        switch(strategy) {
            case 'C': // Call
                payoff = max(S - K1, 0.0) - premium;
                break;
            case 'P': // Put
                payoff = max(K1 - S, 0.0) - premium;
                break;
            case 'S': // Vertical spread
                if (K1 < K2) { // Bull spread
                    payoff = max(S - K1, 0.0) - max(S - K2, 0.0) - premium;
                } else { // Bear spread
                    payoff = max(K1 - S, 0.0) - max(K2 - S, 0.0) - premium;
                }
                break;
            case 'V': // Straddle
                payoff = max(S - K1, 0.0) + max(K1 - S, 0.0) - premium;
                break;
            case 'B': // Butterfly
                payoff = max(S - K1, 0.0) - 2 * max(S - K2, 0.0) + max(S - K3, 0.0) - premium;
                break;
        }
        
        int y_pos = height/2 - static_cast<int>(payoff * y_scale);
        if (y_pos >= 0 && y_pos < height) {
            grid[y_pos][i] = '*';
        }
    }
    
    // Print grid
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            cout << grid[i][j];
        }
        cout << endl;
    }
    
    // Print labels
    cout << "\nStock Price: " << S_min << " to " << S_max << endl;
    cout << "Strike Price(s): " << K1;
    if (K2 != 0) cout << ", " << K2;
    if (K3 != 0) cout << ", " << K3;
    cout << endl;
}

// Monte Carlo option pricing
double monteCarloOptionPrice(char optionType, double S, double K, double T, double r, double sigma, int simulations = 10000) {
    if (simulations < 1) return 0.0;
    
    double sumPayoffs = 0.0;
    double dt = T;
    
    for (int i = 0; i < simulations; i++) {
        double z = sqrt(-2.0 * log((rand() + 1.0) / (RAND_MAX + 1.0))) * cos(2.0 * PI * rand() / (RAND_MAX + 1.0));
        double ST = S * exp((r - 0.5 * sigma * sigma) * dt + sigma * sqrt(dt) * z);
        
        if (optionType == 'C') {
            sumPayoffs += max(ST - K, 0.0);
        } else {
            sumPayoffs += max(K - ST, 0.0);
        }
    }
    
    double discountedPrice = (sumPayoffs / simulations) * exp(-r * T);
    return discountedPrice;
}

// Main menu
void displayMenu() {
    cout << "\n=== Option Trading Theory Calculator ===" << endl;
    cout << "1. Price European Option (Black-Scholes)" << endl;
    cout << "2. Price European Option (Binomial Model)" << endl;
    cout << "3. Price European Option (Monte Carlo)" << endl;
    cout << "4. Calculate Implied Volatility" << endl;
    cout << "5. Calculate Option Greeks" << endl;
    cout << "6. Plot Option Strategy Payoff Diagram" << endl;
    cout << "7. Compare Pricing Models" << endl;
    cout << "8. Exit" << endl;
    cout << "Enter your choice: ";
}

int main() {
    srand(static_cast<unsigned>(time(0)));
    
    int choice;
    char optionType;
    double S, K, T, r, sigma, marketPrice;
    int steps, simulations;
    
    do {
        displayMenu();
        cin >> choice;
        
        if (cin.fail()) {
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
            cout << "Invalid input. Please try again." << endl;
            continue;
        }
        
        switch(choice) {
            case 1: {
                cout << "\nEnter option type (C/P): ";
                cin >> optionType;
                cout << "Enter stock price: ";
                cin >> S;
                cout << "Enter strike price: ";
                cin >> K;
                cout << "Enter time to expiration (years): ";
                cin >> T;
                cout << "Enter risk-free rate: ";
                cin >> r;
                cout << "Enter volatility: ";
                cin >> sigma;
                
                optionType = toupper(optionType);
                double price = blackScholes(optionType, S, K, T, r, sigma);
                cout << "\nBlack-Scholes Option Price: " << fixed << setprecision(2) << price << endl;
                break;
            }
            
            case 2: {
                cout << "\nEnter option type (C/P): ";
                cin >> optionType;
                cout << "Enter stock price: ";
                cin >> S;
                cout << "Enter strike price: ";
                cin >> K;
                cout << "Enter time to expiration (years): ";
                cin >> T;
                cout << "Enter risk-free rate: ";
                cin >> r;
                cout << "Enter volatility: ";
                cin >> sigma;
                cout << "Enter number of steps: ";
                cin >> steps;
                
                optionType = toupper(optionType);
                double price = binomialOptionPrice(optionType, S, K, T, r, sigma, steps);
                cout << "\nBinomial Model Option Price: " << fixed << setprecision(2) << price << endl;
                break;
            }
            
            case 3: {
                cout << "\nEnter option type (C/P): ";
                cin >> optionType;
                cout << "Enter stock price: ";
                cin >> S;
                cout << "Enter strike price: ";
                cin >> K;
                cout << "Enter time to expiration (years): ";
                cin >> T;
                cout << "Enter risk-free rate: ";
                cin >> r;
                cout << "Enter volatility: ";
                cin >> sigma;
                cout << "Enter number of simulations: ";
                cin >> simulations;
                
                optionType = toupper(optionType);
                double price = monteCarloOptionPrice(optionType, S, K, T, r, sigma, simulations);
                cout << "\nMonte Carlo Option Price: " << fixed << setprecision(2) << price << endl;
                break;
            }
            
            case 4: {
                cout << "\nEnter option type (C/P): ";
                cin >> optionType;
                cout << "Enter stock price: ";
                cin >> S;
                cout << "Enter strike price: ";
                cin >> K;
                cout << "Enter time to expiration (years): ";
                cin >> T;
                cout << "Enter risk-free rate: ";
                cin >> r;
                cout << "Enter market option price: ";
                cin >> marketPrice;
                
                optionType = toupper(optionType);
                double iv = impliedVolatility(optionType, S, K, T, r, marketPrice);
                cout << "\nImplied Volatility: " << fixed << setprecision(4) << iv << " (" << iv*100 << "%)" << endl;
                break;
            }
            
            case 5: {
                cout << "\nEnter option type (C/P): ";
                cin >> optionType;
                cout << "Enter stock price: ";
                cin >> S;
                cout << "Enter strike price: ";
                cin >> K;
                cout << "Enter time to expiration (years): ";
                cin >> T;
                cout << "Enter risk-free rate: ";
                cin >> r;
                cout << "Enter volatility: ";
                cin >> sigma;
                
                optionType = toupper(optionType);
                Greeks greeks = calculateGreeks(optionType, S, K, T, r, sigma);
                
                cout << "\nOption Greeks:" << endl;
                cout << "Delta: " << fixed << setprecision(4) << greeks.delta << endl;
                cout << "Gamma: " << greeks.gamma << endl;
                cout << "Theta: " << greeks.theta << " per year" << endl;
                cout << "Vega:  " << greeks.vega << " per 1% change in vol" << endl;
                cout << "Rho:   " << greeks.rho << " per 1% change in rate" << endl;
                break;
            }
            
            case 6: {
                char strategy;
                double K1, K2 = 0, K3 = 0, premium = 0;
                double S_min, S_max;
                
                cout << "\nEnter strategy (C=Call, P=Put, S=Spread, V=Straddle, B=Butterfly): ";
                cin >> strategy;
                cout << "Enter minimum stock price for diagram: ";
                cin >> S_min;
                cout << "Enter maximum stock price for diagram: ";
                cin >> S_max;
                cout << "Enter strike price 1: ";
                cin >> K1;
                
                strategy = toupper(strategy);
                
                if (strategy == 'S' || strategy == 'B') {
                    cout << "Enter strike price 2: ";
                    cin >> K2;
                }
                if (strategy == 'B') {
                    cout << "Enter strike price 3: ";
                    cin >> K3;
                }
                if (strategy != 'C' && strategy != 'P') {
                    cout << "Enter net premium paid: ";
                    cin >> premium;
                }
                
                cout << "\nPayoff Diagram:" << endl;
                plotPayoffDiagram(strategy, S_min, S_max, K1, K2, K3, premium);
                break;
            }
            
            case 7: {
                cout << "\nEnter option type (C/P): ";
                cin >> optionType;
                cout << "Enter stock price: ";
                cin >> S;
                cout << "Enter strike price: ";
                cin >> K;
                cout << "Enter time to expiration (years): ";
                cin >> T;
                cout << "Enter risk-free rate: ";
                cin >> r;
                cout << "Enter volatility: ";
                cin >> sigma;
                cout << "Enter number of steps for binomial model: ";
                cin >> steps;
                cout << "Enter number of simulations for Monte Carlo: ";
                cin >> simulations;
                
                optionType = toupper(optionType);
                double bsPrice = blackScholes(optionType, S, K, T, r, sigma);
                double binPrice = binomialOptionPrice(optionType, S, K, T, r, sigma, steps);
                double mcPrice = monteCarloOptionPrice(optionType, S, K, T, r, sigma, simulations);
                
                cout << "\nModel Comparison:" << endl;
                cout << "Black-Scholes: " << fixed << setprecision(4) << bsPrice << endl;
                cout << "Binomial (" << steps << " steps): " << binPrice << endl;
                cout << "Monte Carlo (" << simulations << " sims): " << mcPrice << endl;
                cout << "Difference (BS - Bin): " << (bsPrice - binPrice) << endl;
                cout << "Difference (BS - MC): " << (bsPrice - mcPrice) << endl;
                break;
            }
            
            case 8:
                cout << "Exiting program..." << endl;
                break;
                
            default:
                cout << "Invalid choice. Please try again." << endl;
        }
        
    } while (choice != 8);
    
    return 0;
}

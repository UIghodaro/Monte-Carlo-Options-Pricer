#include <iostream>
#include <cmath>
#include <random>                                          

// compilation: g++ -O3 -o simulation.out simulation.cpp 
// running:     ./simulation.out

// Quickly implement CDF for Black-Scholes computation
// The error function is included in cmath, use its relationship with the cdf to not need to use other library calls
double CDF(double x) {
    double out = (1 + std::erf(x/std::sqrt(2)))/2;
    return out;
}

// Black-Scholes closed-form equation
// Computes the value of a call option using:
//      S - Stock price, 'the value of the stock at a time t' 
//      K - Exercise/Strike price, 'how much we must pay for the option at expiry time T' [Fixed]
//      r - Risk-free interest rate                                                       [Assume Fixed]
//      T - Time till expiry date                                                         [Fixed under European]
//  Sigma - Volatility of a stock, that is, the standard deviation/variance of the stock from some mean price (Normal distribution) [Assume Fixed] 

// Full equation is C = S*CDF(d1) - K*exp(-rT)*CDF(d2) where d1 and d2 are measures affected mainly by stock price and its standard deviation sigma (I say mainly as K is fixed under European, based on T, so ln(S/K) only varies by dS)
double BlackScholesMertonCall(double S, double K, double r, double T, double sigma) {
    // Needs to calculate d1, d2, their CDFs, then you can just plug in numbers? Or am I allowing time variation?

    double call;                                                    // Eventual output

    double visor = sigma*std::sqrt(T);                              // Separate out here due to being involved in both

    // Complete sigma^2 as sigma*sigma to avoid excess library calls? After studying compiler design this just makes more sense
    // Can't do log this way sadly
    double d1 = (std::log(S/K) + (r + T*(sigma*sigma / 2)))/visor;
    double d2 = d1 - visor;

    // Use the CDF function with these to compute Black-Scholes
    call = S*CDF(d1) - K*std::exp(-r*T)*CDF(d2);
    
    return call;
}

// Letting N being number of Monte-Carlo epochs
// This needs to use the fundamental Brownian motion equation which the Black-Scholes-Merton equation finds the closed form expectation of
// Monte-Carlo methods just abuses the Law of Large Numbers

// We use the GBM equation: S(T) = S(0)*exp((r - sigma^2/2)T + sigma*sqrt{T}*Z(i)) over i from 1 -> N iterations
// Optimise I guess, S(T) = S(0) * exp((r - sigma^2/2)T + sigma*sqrt{T}*Z(i)) = S(0)*exp((r - sigma^2/2)*T)* exp(sigma*sqrt{T}*Z(i))

// Z here is a value picked from a normal distribution of mean 0 and variance 1
double monteCarloCall(double S, double K, double r, double T, double sigma, int N) {
    // Generate N values for Z, then compute S(T) at that value of Z

    double accum = 0.0;                                                     // Cumulative S(T) value, divide by N after
    double call;                                                            // Eventual output
    
    // Generates Z values
    std::random_device rd;
    std::mt19937 gen(rd());                                                 // 'Mersenne Twister' apparently good
    std::normal_distribution<double> dist(0, 1);                            // Holy acrobatics
    
    // Variables for compute with Z, reduce power of the for loop
    double outerMult = S*std::exp((r - (sigma*sigma/2))*T);                 // This optimises into a number, far easier compute overall
    double innerMult = sigma*std::sqrt(T);                                  // Same line of reasoning
    
    for (int i = 0; i < N; i++) {
        double Z = dist(gen);
        accum += std::fmax(outerMult*std::exp(innerMult*Z) - K, 0);         // accumulate MAX(S(T) - K, 0)
    }

    call = std::exp(-r * T)*(accum/N);

    return call;
}

int main() {
    // Configurable inputs
    double S = 1; double K = 1; double r = 1; double T = 1; double sigma = 1; 
    //int N = 100000; 

    std::cout << "The Black-Scholes-Merton output of this configuration is: " << BlackScholesMertonCall(S,K,r,T,sigma) << std::endl;
    //std::cout << "The Monte-Carlo output of this configuration over " << N << " epochs is: " << monteCarloCall(S,K,r,T,sigma, N) << std::endl;
    std::cout << "The Monte-Carlo output of this configuration over 1000 epochs is: " << monteCarloCall(S, K, r, T, sigma, 1000) << "\n";
    std::cout << "The Monte-Carlo output of this configuration over 10000 epochs is: " << monteCarloCall(S, K, r, T, sigma, 10000) << "\n";
    std::cout << "The Monte-Carlo output of this configuration over 100000 epochs is: " << monteCarloCall(S, K, r, T, sigma, 100000) << "\n";
    std::cout << "The Monte-Carlo output of this configuration over 1000000 epochs is: " << monteCarloCall(S, K, r, T, sigma, 1000000) << "\n";
    return 0;
}
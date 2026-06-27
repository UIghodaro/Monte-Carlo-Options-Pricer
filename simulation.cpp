#include <iostream>
#include <cmath>
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

    double BSM;                                                   // Eventual output

    double visor = sigma*std::sqrt(T);                              // Separate out here due to being involved in both

    // Complete sigma^2 as sigma*sigma to avoid excess library calls? After studying compiler design this just makes more sense
    // Can't do log this way sadly
    double d1 = (std::log(S/K) + (r + T*(sigma*sigma / 2)))/visor;
    double d2 = d1 - visor;

    // Use the CDF function with these to compute Black-Scholes
    BSM = S*CDF(d1) - K*std::exp(-r*T)*CDF(d2);
    
    return BSM;
}

// Letting N being number of Monte-Carlo epochs
// This needs to use the fundamental Brownian motion equation which the Black-Scholes-Merton equation finds the closed form expectation of
// Monte-Carlo methods just abuses the Law of Large Numbers

// We use the GBM equation: S(T) = S(0)*exp()
double monteCarloCall(double S, double K, double r, double T, double sigma, int N) {
    return 0;
}

int main(double S, double K, double r, double T, double sigma) {
    std::cout << "The Black-Scholes-Merton output of this configuration is: " << BlackScholesMertonCall(1,1,1,1,1) << std::endl;
    return 0;
}
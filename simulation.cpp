// Baseline
#include <iostream>

// Mathematics 
#include <cmath>        
#include <random>   

// Threadpool
#include <thread>
#include <vector>      

// Testing and Export
#include <chrono>
#include <fstream>
#include <sstream>
#include <string>

// compilation: g++ -O3 -o simulation.out simulation.cpp 
// running:     ./simulation.out

// =========================================
// 1. Black-Scholes closed-form equation
// =========================================

// The overarching options-pricing aim. The goal in any market is estimating the best cost of a stock to make profit - or at least to not get cooked

// Quickly implement CDF for Black-Scholes computation
// The error function is included in cmath, use its relationship with the cdf to not need to use other library calls
double CDF(double x) {
    double out = (1 + std::erf(x/std::sqrt(2)))/2;
    return out;
}


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
    double d1 = (std::log(S/K) + T*(r + (sigma*sigma / 2)))/visor;
    double d2 = d1 - visor;

    // Use the CDF function with these to compute Black-Scholes
    call = S*CDF(d1) - K*std::exp(-r*T)*CDF(d2);
    
    return call;
}

// =========================================
// 2. Monte-Carlo method approximation
// =========================================

// We will not always have access to an easy closed-form expected value equation. Monte-Carlo methods gets us close to its solution 

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
    std::mt19937 gen(rd());                                                 // 'Mersenne Twister' randomiser, apparently good
    std::normal_distribution<double> dist(0, 1);                            // Yapper
    
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

// =========================================
// Multithreading for speed and accuracy
// =========================================

// Increase the accuracy of Monte-Carlo while keeping the compute time the same (or even increase the compute time for drastically better accuracy)

// Now, we attempt the above but with use of multithreading.
// The idea: Increasing epochs N increases accuracy at the cost of compute time
//           Smaller number of epochs = Faster Compute
//           What if we do a smaller number of epochs many times at once and combine that to effectively get many epochs? 
//           A simple concept, the improvement is likely amazing

// Architecture:
// each thread completes N/thread_count epochs, creating their own partial cumulative sum
// join worker threads, then process finishes in the same way as before - divide cumulative sum by N, then discount it to get your answer

// This could be redone even easier by simply reusing monteCarloCall and taking an average of each thread's Monte Carlo. 
// In fact the difference in speed should be negligible I might test that
double threadPoolMonteCarlo(double S, double K, double r, double T, double sigma, int N){
    
    double accum = 0.0; 
    double call;

    unsigned int numThreads = std::thread::hardware_concurrency();  // Draw threads from available hardware
    std::vector<double> partials(numThreads, 0.0);                  // Initialise a vector with length numThreads, set each index to 0.0
    std::vector<std::thread> workers;                               // Start a thread vector

    int chunkSize = N / numThreads;

    // These will remain the same for every thread computation
    double outerMult = S*std::exp((r - (sigma*sigma/2))*T);                 
    double innerMult = sigma*std::sqrt(T);                                  

    // Worker logic
    for(unsigned int t = 0; t < numThreads; t++){
        
        // The following is implemented for the case where chunk size is not a clean division
        int start = t * chunkSize;                                  // Start of a worker thread's set of epochs to do                           
        int end = (t == numThreads - 1) ? N : start + chunkSize;    // End of a worker thread's set of epochs to do, if it is the final thread, take it all the way to N
        int threadEpochs = end - start;                             // actual number

        // Pass t and threadEpochs into each thread since each thread uses them, right?
        workers.push_back(std::thread([&, t, threadEpochs]() {
            std::random_device rd;
            std::mt19937 gen(rd());                                                 
            std::normal_distribution<double> dist(0, 1);                           
    
            double p_accum = 0.0;

            for (int i = 0; i < threadEpochs; i++) {
                double Z = dist(gen);
                p_accum += std::fmax(outerMult*std::exp(innerMult*Z) - K, 0);
            }

            partials[t] = p_accum;                  // Remember this is self contained
                                                    // each thread will create an instance of partials as partials[i] = 0.0 for i != t and partials[t] = accum 
        }));

    }

    // Join the processes, complete the partials vector to have each cumulative sum
    for (auto& worker : workers) {
        worker.join();
    }

    // Join the partial sums
    for (auto& partial : partials) {
        accum += partial;
    }

    // Return the call option price
    call = std::exp(-r * T)*(accum/N);
    return call;

}

// ==========================================
// Side thing - Visualise the convergence
// ==========================================

// Just export it to a CSV and then visualise using a matplotlib, no use using cpp

struct Entry {
    std::string epochs;
    std::string call;
};

// Ascending powers of 10 for epochs
void visualiseConverge(double S, double K, double r, double T, double sigma, int power) {

    int N = 1;
    std::vector<Entry> entries;

    for(int i = 0; i < power+1; i++) {
        for (int mult : {1, 2, 5}) {
            int N = mult * std::pow(10, i);
            entries.push_back({std::to_string(N), std::to_string(threadPoolMonteCarlo(S,K,r,T,sigma,N))});
        }
    }

    std::ofstream csvFile("out.csv");
    for (const auto& entry : entries) {
        csvFile << entry.epochs << "," << entry.call << "\n";
    }

    csvFile.close();
}

// ==========================================
// Main, printing output and presenting error
// ==========================================

int main() {
    // For displaying time performance
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    // Configurable inputs
    double S = 1; double K = 1; double r = 1; double T = 1; double sigma = 1; 
    int N = 100000000; // default of 100 million

    std::cout << "There are " << std::thread::hardware_concurrency() << " threads available" << std::endl;

    // Computing values, their compute times and error
    auto BSMs  = high_resolution_clock::now();
    double BSM = BlackScholesMertonCall(S,K,r,T,sigma);
    auto BSMf  = high_resolution_clock::now();
    auto ms1 = duration_cast<milliseconds>(BSMf - BSMs);

    auto MCMs  = high_resolution_clock::now();
    double MCM = monteCarloCall(S,K,r,T,sigma, N);
    auto MCMf  = high_resolution_clock::now();
    auto ms2   = duration_cast<milliseconds>(MCMf - MCMs);
    double acc1 = std::abs(BSM - MCM);

    auto TPMCMs  = high_resolution_clock::now();
    double TPMCM = threadPoolMonteCarlo(S,K,r,T,sigma, N);
    auto TPMCMf  = high_resolution_clock::now();
    auto ms3     = duration_cast<milliseconds>(TPMCMf - TPMCMs);
    double acc2  = std::abs(BSM - TPMCM);

    std::cout << "The Black-Scholes-Merton output of this configuration is: " << BSM << " in " << ms1.count() << "ms\n";
    std::cout << "The Monte-Carlo output of this configuration over " << N << " epochs is: " << MCM << " in " << ms2.count() << "ms [Error = " << acc1 << "]\n";
    std::cout << "The TP Monte-Carlo output of this configuration over " << N << " epochs is: " << TPMCM << " in " << ms3.count() << "ms [Error = " << acc2 << "]\n ... \n";
    
    // Get visualiser csv
    auto t1 = high_resolution_clock::now();
    visualiseConverge(S, K, r, T, sigma, 8);
    auto t2 = high_resolution_clock::now();
    auto ms = duration_cast<milliseconds>(t2-t1);

    std::cout << "Datapoints CSV generated in " << ms.count() << "ms\n";

    /* Getting number of milliseconds as a double. */
    // duration<double, std::milli> ms_double = t2 - t1;

    return 0;
}


// Previously used outputs, visualising the convergence increase over number of epochs
    //std::cout << "The Monte-Carlo output of this configuration over 1000 epochs is: " << monteCarloCall(S, K, r, T, sigma, 1000) << "\n";
    //std::cout << "The Monte-Carlo output of this configuration over 10000 epochs is: " << monteCarloCall(S, K, r, T, sigma, 10000) << "\n";
    //std::cout << "The Monte-Carlo output of this configuration over 100000 epochs is: " << monteCarloCall(S, K, r, T, sigma, 100000) << "\n";
    //std::cout << "The Monte-Carlo output of this configuration over 1000000 epochs is: " << monteCarloCall(S, K, r, T, sigma, 1000000) << "\n";
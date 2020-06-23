//LSN_Exercise_03
//Exercise 03.1 - Plain vanilla option pricing

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <numeric>
#include <iomanip>
#include "Parallel_Generator/random.h"

using namespace std;


//Function for statistical uncertainty estimation
double error(const std::vector<double>& AV, const std::vector<double>& AV2, const int& n) {

    return (n == 0) ? 0. : sqrt((AV2[n] - pow(AV[n], 2))/n);
}


//Print results to file
template <typename T, typename U>
void print_to_file(const std::string& name, const std::vector<T>& v1, const std::vector<U>& v2, int prec = 5) {
    ofstream fileout;
    int size = std::min(v1.size(), v2.size());

    if (v1.size() != v2.size())
        cout << "WARNING: function print_to_file received vectors with different size." << endl;

    fileout.open(name);
    if (fileout.is_open()) {
        for (int i = 0; i < size; i++) {
            fileout << std::setprecision(prec) << std::fixed;
            fileout << i << " " << v1[i] << " " << v2[i] << endl;
        }
    }
    else {
        cerr << "ERROR: Unable to write file " << name << endl;
    }

    fileout.close();
}

//Function for initializing random generator - if it fails then program is terminated
void init_random_gen(Random& gen) {
    int seed[4];
    int p1, p2;
    bool success = true;

    ifstream Primes("Parallel_Generator/Primes");
    if (Primes.is_open())
       Primes >> p1 >> p2;
    else {
        cerr << "PROBLEM: Unable to open Primes" << endl;
        success = false;
    }
    Primes.close();

    ifstream input("Parallel_Generator/seed.in");
    string property;
    if (input.is_open()) {
       while ( !input.eof() ) {
          input >> property;
          if( property == "RANDOMSEED" ) {
             input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
             gen.SetRandom(seed, p1, p2);
          }
       }
       input.close();
    }
    else {
        cerr << "PROBLEM: Unable to open seed.in" << endl;
        success = false;
    }

    if (!success) {
        cerr << "ERROR: Unable to initialize random generator." << endl;
        cerr << "---STOP---" << endl;
        exit (EXIT_FAILURE);
    }
}

/*************************************MAIN************************************/

int main() {
    const int M = 1000000;         //steps
    const int N = 100;             //blocks
    const int L = M/N;             //steps per block
    const int intervals = 100;     //number of intervals for dicretized GBR

    const double S0 = 100.;        //asset price at t = 0
    const double T = 1.;           //delivery time
    const double K = 100.;         //strike price
    const double r = 0.1;          //risk-free interest rate
    const double sigma = 0.25;     //volatility

    std::vector<double> call_put(N), call_put2(N);      //Vector with call or put option
    std::vector<double> sum_prog(N), su2_prog(N), err_prog(N);

    Random rnd;
    init_random_gen(rnd);

//Run simulation
/*--------------------European call-option prices - DIRECT--------------------*/
    for (int i = 0; i < N; i++) {
        double S = 0., W = 0., sum_call = 0.;

        for (int j = 0; j < L; j++) {
            W = rnd.Gauss(0., T);
            S = S0*exp((r - sigma*sigma/2)*T + sigma*W);
            sum_call += std::max(0., S - K)*exp(-r*T);
        }
        call_put[i] = sum_call/L;
        call_put2[i] = pow(call_put[i], 2);
    }

    for (int i = 0; i < N; i++) {
        sum_prog[i] = std::accumulate(call_put.begin(), call_put.begin()+i+1, 0., std::plus<double>())/(i+1);
        su2_prog[i] = std::accumulate(call_put2.begin(), call_put2.begin()+i+1, 0., std::plus<double>())/(i+1);
        err_prog[i] = error(sum_prog, su2_prog, i);    //Statistical uncertainty
    }

    cout << "END: call direct." << endl;

//Print results to file
    print_to_file("call_direct.dat", sum_prog, err_prog);

/*--------------------European put-option prices - DIRECT---------------------*/
    for (int i = 0; i < N; i++) {
        double S = 0., W = 0., sum_put = 0.;

        for (int j = 0; j < L; j++) {
            W = rnd.Gauss(0., T);
            S = S0*exp((r - sigma*sigma/2)*T + sigma*W);
            sum_put += std::max(0., K - S)*exp(-r*T);
        }
        call_put[i] = sum_put/L;
        call_put2[i] = pow(call_put[i], 2);
    }

    for (int i = 0; i < N; i++) {
        sum_prog[i] = std::accumulate(call_put.begin(), call_put.begin()+i+1, 0., std::plus<double>())/(i+1);
        su2_prog[i] = std::accumulate(call_put2.begin(), call_put2.begin()+i+1, 0., std::plus<double>())/(i+1);
        err_prog[i] = error(sum_prog, su2_prog, i);    //Statistical uncertainty
    }

    cout << "END: put direct." << endl;

//Print results to file
    print_to_file("put_direct.dat", sum_prog, err_prog);


/*-------------European call-option prices - DISCRETIZED--------------------*/
    for (int i = 0; i < N; i++) {
        double Z = 0., sum_call = 0.;

        for (int j = 0; j < L; j++) {
            double S_ti = S0;

            for (int k = 0; k < intervals; k++) {
                Z = rnd.Gauss(0., 1.);
                S_ti = S_ti*exp((r - sigma*sigma/2)*T/intervals + sigma*Z*sqrt(T/intervals));   //t_(i+1) - t_i is T/intervals (length of the interval)
            }
            sum_call += std::max(0., S_ti - K)*exp(-r*T);
        }
        call_put[i] = sum_call/L;
        call_put2[i] = pow(call_put[i], 2);
    }

    for (int i = 0; i < N; i++) {
        sum_prog[i] = std::accumulate(call_put.begin(), call_put.begin()+i+1, 0., std::plus<double>())/(i+1);
        su2_prog[i] = std::accumulate(call_put2.begin(), call_put2.begin()+i+1, 0., std::plus<double>())/(i+1);
        err_prog[i] = error(sum_prog, su2_prog, i);    //Statistical uncertainty
    }

    cout << "END: call discretized." << endl;

//Print results to file
    print_to_file("call_discretized.dat", sum_prog, err_prog);


/*---------------European put-option prices - DISCRETETIZED-------------------*/
    for (int i = 0; i < N; i++) {
        double Z = 0., sum_put = 0.;

        for (int j = 0; j < L; j++) {
            double S_ti = S0;

            for (int k = 0; k < intervals; k++) {
                Z = rnd.Gauss(0., 1.);
                S_ti = S_ti*exp((r - sigma*sigma/2)*T/intervals + sigma*Z*sqrt(T/intervals));   //t_(i+1) - t_i is T/intervals (length of the interval)
            }
            sum_put += std::max(0., K - S_ti)*exp(-r*T);
        }
        call_put[i] = sum_put/L;
        call_put2[i] = pow(call_put[i], 2);
    }

    for (int i = 0; i < N; i++) {
        sum_prog[i] = std::accumulate(call_put.begin(), call_put.begin()+i+1, 0., std::plus<double>())/(i+1);
        su2_prog[i] = std::accumulate(call_put2.begin(), call_put2.begin()+i+1, 0., std::plus<double>())/(i+1);
        err_prog[i] = error(sum_prog, su2_prog, i);    //Statistical uncertainty
    }

    cout << "END: put discretized." << endl;

//Print results to file
    print_to_file("put_discretized.dat", sum_prog, err_prog);


    return 0;
}

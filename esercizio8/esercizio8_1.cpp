//LSN_Exercise_08
//Exercise 08.1

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <numeric>
#include <iomanip>
#include "Parallel_Generator/random.h"

using namespace std;

int n_accepted = 0;        //number of accepted moves
int n_trials = 0;          //number of toval moves
double param = 0.;         //parameter for 50% acceptance (Metropolis)
double m = 1.;             //particle mass

//Trial wave function
double psi_trial (double x, double mu, double sigma) {
    double Zm = (x - mu)/sigma;
    double Zp = (x + mu)/sigma;

    return (exp(-Zm*Zm/2.) + exp(-Zp*Zp/2.));
}

//Probability density function (trial) - normalized (even though is unnecessary)
double psi_density (double x, double mu, double sigma) {
    double Zm = (x - mu)/sigma;
    double Zp = (x + mu)/sigma;
    double norm = 2.*sqrt(M_PI)*sigma*(1. + exp(-mu*mu/(sigma*sigma)));

    return ( exp(-Zm*Zm) + exp(-Zp*Zp) + 2.*exp(-(x*x + mu*mu)/(sigma*sigma)) ) / norm;
}

//Transiotion probability T(x_new | x_old) - uniform
double T_uniform (Random& rnd, double x, double delta) {

    return rnd.Rannyu(x - delta, x + delta);
}

//Metropolis algorithm
void Metropolis (std::vector<double>& x, const double& mu, const double& sigma, Random& rnd, std::function<double(double, double, double)> psi,
                 std::function<double(Random&, double, double)> T) {

    double x_end = *(x.end() - 1);

//Generate trial move
    double x_trial = T(rnd, x_end, param);

    double accept = std::min(1., psi(x_trial, mu, sigma)/psi(x_end, mu, sigma));

    if (rnd.Rannyu() <= accept) {   //accept the move
        x.push_back(x_trial);

        n_accepted++;
    }
    else {                          //reject the move
        x.push_back(x_end);
    }

    n_trials++;
}

//Kinetic part of H with hbar = 1
double H_kinetic (double x, double mu, double sigma) {
    double Zm = (x - mu)/sigma;
    double Zp = (x + mu)/sigma;

    return ( exp(-Zm*Zm/2.)*(1. - Zm*Zm) + exp(-Zp*Zp/2.)*(1. - Zp*Zp) )/(2.*m*sigma*sigma);
}

//Potential part of H
double H_potential (double x, double mu, double sigma) {

    return (x*x - 2.5)*x*x*psi_trial(x, mu, sigma);
}

//Hamiltonia expectation value
double H_expect_val (double x, double mu, double sigma) {

    return (H_kinetic(x, mu, sigma) + H_potential(x, mu, sigma))/psi_trial(x, mu, sigma);
}

//Function for statistical uncertainty estimation
double error (const std::vector<double> & AV, const std::vector<double> & AV2, const int & n) {

    return (n == 0) ? 0. : sqrt( (AV2[n] - pow(AV[n], 2))/n );
}


//Function for initializing random generator - if it fails then program is terminated
void init_random_gen (Random& gen) {
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

//Print results to file
template <typename T, typename U>
void print_to_file (const std::string & name, const std::vector<T> & v1, const std::vector<U> & v2, int prec = 5) {
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

//Print a single vector to file
template <typename T>
void print_vector (const std::string & name, const std::vector<T> & v) {
    ofstream fileout;

    fileout.open(name);
    if (fileout.is_open()) {
        for (T elem : v)
            fileout << elem << endl;
    }
    else {
        cerr << "ERROR: Unable to write file " << name << endl;
    }

    fileout.close();
}

/*************************************MAIN*************************************/
int main (int argc, char const *argv[]) {

    if (argc != 4){     //check parameters
        cerr << endl << "Wrong  call " << argv[0] << endl;
        cerr << "The program needs <param1: mu> <param2: sigma> <param3: delta_metropolis> " << endl << endl;
        return 1;
    }

    const int N = 50, L = 10000;     //number of blocks and steps per block
    double mu = atof(argv[1]), sigma = atof(argv[2]);

    std::vector<double> sample;    //trial wave function sample
    std::vector<double> ave_H(N), ave2_H(N);
    std::vector<double> sum_prog(N), su2_prog(N), err_prog(N);  //vectors for data blocking

    Random rnd;

//Initialize the simulation
    init_random_gen(rnd);
    sample.push_back(0.);
    param = atof(argv[3]);

//Use Metropolkis to generate the Sample
    for (int i = 0; i < N*L; i++)
        Metropolis(sample, mu, sigma, rnd, psi_density, T_uniform);

//Estimate <H>_T
    for (int i = 0; i < N; i++) {
        int k = 0;
        double sum_H = 0.;

        for (int j = 0; j < L; j++) {
            k = j + i*L;
            sum_H += H_expect_val(sample[k], mu, sigma);
        }
        ave_H[i] = sum_H/L;             //<H> for each block
        ave2_H[i] = pow(ave_H[i], 2);   //<H>^2
    }

//Compute mean and error on the estimation of <H>_T - data blocking
    for (int i = 0; i < N; i++) {
        sum_prog[i] = accumulate(ave_H.begin(), ave_H.begin()+i+1, 0., std::plus<double>())/(i+1);
        su2_prog[i] = accumulate(ave2_H.begin(), ave2_H.begin()+i+1, 0., std::plus<double>())/(i+1);
        err_prog[i] = error(sum_prog, su2_prog, i);    //Statistical uncertainty
    }

//Print results on screen
    cout << "Acceptance: " << (double)n_accepted/n_trials << " Mu: " << mu << " sigma: " << sigma << endl;
    cout << "E_new: " << *(sum_prog.end() - 1) << " Error: " << *(err_prog.end() - 1) << endl << endl;

//Print results to file
    print_to_file("ave_H.dat", sum_prog, err_prog);
    print_vector("sample_trial.dat", sample);

    return 0;
}

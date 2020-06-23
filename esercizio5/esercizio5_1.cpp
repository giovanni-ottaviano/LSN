//LSN_Exercise_05
//Exercise 05.1 - Hydrogen atom

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <functional>
#include <numeric>
#include <iomanip>
#include <cassert>
#include "Parallel_Generator/random.h"

using namespace std;

int n_accept = 0;       //number of accepted moves
int n_metropolis = 0;   //number of toval moves

//Function for statistical uncertainty estimation
double error(const std::vector<double> & AV, const std::vector<double> & AV2, const int & n) {

    return (n == 1) ? 0. : sqrt((AV2[n] - pow(AV[n], 2))/(n - 1.));
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

//Print results to file
template <typename T, typename U>
void print_to_file(const std::string & name, const std::vector<T> & v1, const std::vector<U> & v2, int prec = 5) {
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

//Print results to file
template <typename T>
void print_coord_to_file(const std::string & name, const std::vector<T> & x, const std::vector<T> & y, const std::vector<T> & z, int prec = 5) {
    ofstream fileout;
    unsigned int size = x.size();

//vectors with differet sizes -> don't print
    if ( (size != y.size()) || (size != z.size()) ) {
        cerr << "ERROR: function print_coord_to_file cannot print vectors of different sizes." << endl;
        return;
    }

    fileout.open(name);
    if (fileout.is_open()) {
        fileout << std::setprecision(prec) << std::fixed;
        for (unsigned int i = 0; i < size; i++) {
            fileout << x[i] << " " << y[i] << " " << z[i] << endl;
        }
    }
    else {
        cerr << "ERROR: Unable to write file " << name << endl;
    }

    fileout.close();
}


//Probability density function (Bohr radius units)
//|psi(x,y,z)|^2 (cartesian coordinates) n,l,m = 1,0,0
double psi_density100(double x, double y, double z) {
    double r = sqrt(x*x + y*y + z*z);

    return exp(-2.*r)/M_PI;
}

//Probability density function (Bohr radius units)
//|psi(x,y,z)|^2 (cartesian coordinates) n,l,m = 2,1,0
double psi_density210(double x, double y, double z) {
    double r = sqrt(x*x + y*y + z*z);

    return exp(-r)*pow(z, 2)/(32.*M_PI);
}

//Transiotion matrix T(x,y) - uniform
//delta MUST be positive
double T_uniform(Random& rnd, double end, double delta) {
    assert(delta > 0.);

    return rnd.Rannyu(end - delta, end + delta);
}

//Transiotion matrix T(x,y) - multivariate normal
//sigma MUST be positive
double T_multivar(Random& rnd, double end, double sigma) {
    assert(sigma > 0.);

    return rnd.Gauss(end, sigma);
}

//Metropolis algorithm - function perform one single step
//function takes 2 f. pointer: probability density (psi) and transition probability (T)
void Metropolis(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, Random& rnd, std::function<double(double, double, double)> psi,
                std::function<double(Random&, double, double)> T) {

    double x_end = *(x.end() - 1);
    double y_end = *(y.end() - 1);
    double z_end = *(z.end() - 1);
    double param;   //parameter for 50% acceptance

//Set param for required acceptance (different for uniform and multivariate transition)
    if ( (*T.target<double(*)(Random&,double,double)>() ) == T_uniform)
        param = ( (*psi.target<double(*)(double,double,double)>() ) == psi_density100) ? 1.2 : 2.95;
    else
        param = ( (*psi.target<double(*)(double,double,double)>() ) == psi_density100) ? 0.75 : 1.87;

//Generate trial move
    double x_trial = T(rnd, x_end, param);
    double y_trial = T(rnd, y_end, param);
    double z_trial = T(rnd, z_end, param);
    double accept = std::min(1., psi(x_trial, y_trial, z_trial)/psi(x_end, y_end, z_end));

    if (rnd.Rannyu() <= accept) {    //accept the move
        x.push_back(x_trial);
        y.push_back(y_trial);
        z.push_back(z_trial);

        n_accept++;
    }
    else {                           //reject the move
        x.push_back(x_end);
        y.push_back(y_end);
        z.push_back(z_end);
    }

    n_metropolis++;
}


/******************************************************************************/
/************************************MAIN**************************************/
/******************************************************************************/
int main() {
    const int M = 10000000;             //Number of steps
    const int N = 100;                  //Number of blocks
    const int L = M/N;                  //Number of steps per block - M multiple of N
    const int eq_steps = 100;           //Number of steps for equilibration

    std::vector<double> r(N), r2(N);                             //radius, radius^2 in each block
    std::vector<double> sum_prog(N), su2_prog(N), err_prog(N);   //Vector for data blocking

    double start;                           //Starting point
    std::string file_coord, sim_results;    //file names: coordinates, results
    auto T_p = T_uniform;                   //function pointer: transition probability
    auto psi_density = psi_density100;      //function pointer: probability density

    Random rnd;
    init_random_gen(rnd);   //Initialize Random Generator

//Run Metropolis simulation
//Each index correspond to a different setting for the simulation (1s/2p - uniform/gauss multivar)
    for (int sim_index = 0; sim_index < 4; sim_index++) {
        std::vector<double> x, y, z;    //Vectors of coordinates
        n_metropolis = 0;
        n_accept = 0;

        switch (sim_index) {
            case 0:             //1s - uniform transition probability
                start = 1.;     //Starting point (Bohr radius units)
                T_p = T_uniform;
                psi_density = psi_density100;
                file_coord = "coord_1s.dat";
                sim_results = "r_1s.dat";
                break;
            case 1:             //2p - uniform transition probability
                start = 3.;     //Starting point (Bohr radius units)
                T_p = T_uniform;
                psi_density = psi_density210;
                file_coord = "coord_2p.dat";
                sim_results = "r_2p.dat";
                break;
            case 2:             //1s - multivariate normal transition probability
                start = 1.;  //Starting point (Bohr radius units)
                T_p = T_multivar;
                psi_density = psi_density100;
                file_coord = "coord_1s_multivar.dat";
                sim_results = "r_1s_multivar.dat";
                break;
            case 3:             //2p - multivariate normal transition probability
                start = 3.;     //Starting point (Bohr radius units)
                T_p = T_multivar;
                psi_density = psi_density210;
                file_coord = "coord_2p_multivar.dat";
                sim_results = "r_2p_multivar.dat";
                break;
        }

        x.push_back(start);
        y.push_back(start);
        z.push_back(start);

//Equilibration
    double appo_x, appo_y, appo_z;
    for (int i = 0; i < eq_steps; i++)
        Metropolis(x, y, z, rnd, psi_density, T_p);

//Save last element and clean vectors
    appo_x = x.back();
    appo_y = y.back();
    appo_z = z.back();

    x.clear();
    y.clear();
    z.clear();

//Reinitialize coordinates
    x.push_back(appo_x);
    y.push_back(appo_y);
    z.push_back(appo_z);


//Sample psi using Metropolis
        for (int i = 0; i < M; i++)
            Metropolis(x, y, z, rnd, psi_density, T_p);

//Print coordinates to file (!!! for N > 10^7 huge files)
//        print_coord_to_file(file_coord, x, y, z);

//Calculate <r> on each block
        for (int i = 0; i < N; i++) {
            int k = 0;
            double sum_r = 0.;

            for (int j = 0; j < L; j++) {
                k = j + i*L;
                sum_r += sqrt(pow(x[k], 2) + pow(y[k], 2) + pow(z[k], 2));
            }

            r[i] = sum_r/L;         //estimation of r
            r2[i] = pow(r[i], 2);   //r^2
        }

//Compute mean and error on the estimation of <r> - blocking method
        for (int i = 0; i < N; i++) {
            sum_prog[i] = accumulate(r.begin(), r.begin()+i+1, 0., std::plus<double>())/(i+1);
            su2_prog[i] = accumulate(r2.begin(), r2.begin()+i+1, 0., std::plus<double>())/(i+1);
            err_prog[i] = error(sum_prog, su2_prog, i);    //Statistical uncertainty
        }

//Print results to file
        print_to_file(sim_results, sum_prog, err_prog);

        cout << "Percentage of accepted steps (Metropolis): " << ((double)n_accept/n_metropolis) * 100. << endl;
        cout << "Simulation step " << sim_index << " completed." << endl << endl;
    }

    return 0;
}

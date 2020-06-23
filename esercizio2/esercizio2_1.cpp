//LSN_Exercise_02
//Exercise 02.1

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

//function to integrate g(x) = (PI/2)*cos(x*PI/2)
double g(double x) {

    return (M_PI/2.)*cos(x*M_PI/2.);
}

//new probability density (importance sampling)
//Taylor expantion to 4 order(normalized) (il polinomio si annulla in x0 = 1.013...)
double d(double x) {
    double pim = M_PI/2.;

    return (1 - 0.5*pow(x*pim, 2) + pow(x*pim, 4)/24.)/(1 - pow(pim, 2)/6. + pow(pim, 4)/120.);
}

//function to integrate g1(x) = cos(x*PI/2)/d(x) -> importance sampling
double g1(double x) {

    return (M_PI/2.)*cos(x*M_PI/2.)/d(x);
}

//Function for acceptation
bool accept(double x, double r, double pmax = 1.) {

    return (r < d(x)/pmax) ? true : false;
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

/************************************MAIN**************************************/
int main() {
    const int M = 10000000;          //#steps
    const int N = 100;               //#blocks
    const int L = M/N;               //#steps per block
    double x_min = 0., x_max = 1.;   //bounds for integration

    std::vector<double> x0(M);      //distribuzione uniforme di samples [0,1]
    std::vector<double> I(N), I2(N);
    std::vector<double> sum_prog(N), su2_prog(N), err_prog(N);

    Random rnd;
    init_random_gen(rnd);

//fill vector with r in [0,1]
    for (int i = 0; i < M; i++)
        x0[i] = rnd.Rannyu(x_min, x_max);

//run simulation
    for (int i = 0; i < N; i++) {
        int k = 0;
        double sum = 0.;

        for (int j = 0; j < L; j++) {
            k = j + i*L;
            sum += g(x0[k]);
        }

        I[i] = (x_max - x_min)*sum/L;   //integral estimation (in each block)
        I2[i] = pow(I[i], 2);           //square
    }

//sums and error
    for (int i = 0; i < N; i++) {
        sum_prog[i] = accumulate(I.begin(), I.begin()+i+1, 0., std::plus<double>())/(i+1);
        su2_prog[i] = accumulate(I2.begin(), I2.begin()+i+1, 0., std::plus<double>())/(i+1);
        err_prog[i] = error(sum_prog, su2_prog, i);    //Statistical uncertainty
    }

//Print results to file
    print_to_file("integral.dat", sum_prog, err_prog, 5);


/*----------------------------------PART 2-----------------------------------*/
/*---------------------------IMPORTANCE SAMPLING------------------------------*/

//Reset all vectors - just to be sure
    std::fill(I.begin(), I.end(), 0.);
    std::fill(I2.begin(), I2.end(), 0.);
    std::fill(sum_prog.begin(), sum_prog.end(), 0.);
    std::fill(su2_prog.begin(), su2_prog.end(), 0.);
    std::fill(err_prog.begin(), err_prog.end(), 0.);
    x0.clear();

//make a new vector  of random numbers (pdf d(x))
//reset random generator (get same random points)
    init_random_gen(rnd);

    int n_random = 0;
    double p_max = d(0.);   //x = 0. maximise prob density d(x)

    while (n_random != M) {
        double appo_x = rnd.Rannyu(x_min, x_max);
        double appo_r  = rnd.Rannyu();

        if (accept(appo_x, appo_r, p_max)) {
            x0.push_back(appo_x);
            n_random++;
        }
    }

//run simulation
    for (int i = 0; i < N; i++) {
        int k = 0;
        double sum = 0.;

        for (int j = 0; j < L; j++) {
            k = j + i*L;
            sum += g1(x0[k]);
        }
        I[i] = (x_max - x_min)*sum/L;   //integral estimator
        I2[i] = pow(I[i], 2);           //square
    }

//sums and errors
    for (int i = 0; i < N; i++) {
        sum_prog[i] = accumulate(I.begin(), I.begin()+i+1, 0., std::plus<double>())/(i+1);
        su2_prog[i] = accumulate(I2.begin(), I2.begin()+i+1, 0., std::plus<double>())/(i+1);
        err_prog[i] = error(sum_prog, su2_prog, i);    //Statistical uncertainty
    }

//Print results to file
    print_to_file("integral_IS.dat", sum_prog, err_prog, 7);


    return 0;
}

//LSN_Exercise_01
//Exercise 01.1

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <numeric>
#include <iomanip>
#include "Parallel_Generator/random.h"

using namespace std;


//Function for statistical uncertainty estimation (n = # of block)
double error(const std::vector<double>& AV, const std::vector<double>& AV2, int n) {

    return (n == 0) ? 0 : sqrt((AV2[n] - pow(AV[n], 2))/n);
}


//Function for initializing random generator - if it fails then program is terminated
void init_random_gen (Random& gen) {
    int seed[4];
    int p1, p2;
    bool success = true;

    ifstream Primes("Parallel_Generator/Primes");
    if (Primes.is_open()) {
       Primes >> p1 >> p2;
   }
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

//Print single vector
template <typename T>
void print_vector(const std::string& name, const std::vector<T>& v1, int prec = 5) {
    ofstream fileout(name);

    if (fileout.is_open()) {
        for (T elem : v1)
            fileout << std::setprecision(prec) << std::fixed << elem << endl;
    }
    else {
        cerr << "ERROR: Unable to write file " << name << endl;
    }

    fileout.close();
}

/************************************MAIN**************************************/
int main() {

    const int M = 1000000;              //Total number of throws
    const int N = 100;                  //Number of blocks
    int L = M/N;                        //Number of throws in each block

    std::vector<double> r(M);           //U[0,1) uniform distribution
    std::vector<double> ave(N), av2(N);
    std::vector<double> sum_prog(N), su2_prog(N), err_prog(N);

//Initialize random generator
    Random rnd;
    init_random_gen(rnd);

//riempio il vettore di numeri random
    for (int i = 0; i < M; i++)
        r[i] = rnd.Rannyu();

//run simulation
    for (int i = 0; i < N; i++) {
        double sum = 0.;
        int k = 0;

        for (int j = 0; j < L; j++){
            k = j + i*L;
            sum += r[k];
        }
        ave[i] = sum/L;           //r_i
        av2[i] = pow(ave[i], 2);   //(r_i)^2
    }

//compute progressive sums and MC error
    for (int i = 0; i < N; i++) {
        sum_prog[i] = accumulate(ave.begin(), ave.begin()+i+1, 0., std::plus<double>())/(i+1);  //Cumulative average
        su2_prog[i] = accumulate(av2.begin(), av2.begin()+i+1, 0., std::plus<double>())/(i+1);  //Cumulative square average
        err_prog[i] = error(sum_prog,su2_prog,i);    //Statistical uncertainty
    }

//Write results to files
    print_to_file("r.dat", sum_prog, err_prog);


/*---------------------------------PART 2------------------------------------*/

//clean all vectors before starting point 2
    std::fill(ave.begin(), ave.end(), 0.);
    std::fill(av2.begin(), av2.end(), 0.);
    std::fill(sum_prog.begin(), sum_prog.end(), 0.);
    std::fill(su2_prog.begin(), su2_prog.end(), 0.);
    std::fill(err_prog.begin(), err_prog.end(), 0.);

//run simulation
    for (int i = 0; i < N; i++) {
        double sum = 0.;
        int k = 0;

        for (int j = 0; j < L; j++) {
            k = j + i*L;
            sum += pow(r[k] - 0.5, 2);   //Accumulate measures
        }
        ave[i] = sum/L;                 //Estimation for each block
        av2[i] = pow(ave[i], 2);
    }

//compute average value and error - blocking method
    for (int i = 0; i < N; i++) {
        sum_prog[i] = accumulate(ave.begin(), ave.begin()+i+1, 0., std::plus<double>())/(i+1);
        su2_prog[i] = accumulate(av2.begin(), av2.begin()+i+1, 0., std::plus<double>())/(i+1);
        err_prog[i] = error(sum_prog, su2_prog, i);    //Statistical uncertainty
    }

//print results to file
    print_to_file("sigma.dat", sum_prog, err_prog);


//--------------------------------PART 3---------------------------------------
    const int n = 10000;            //points for each iteration (correct value for 100 iter)
    const int MM = 100;             //number of subintervals of [0,1]
    const int n_test = 100;         //numer of chi2 tests

    std::vector<double> chi2;

    for (int i = 0; i < n_test; i++) {
        std::vector<double> n_i(MM);   //random numbers per interval

        for (int k = 0; k < MM; k++) {
            double inf = (double)k/MM, sup = (double)(k + 1)/MM;   //bounds

            for (int j = 0; j < n; j++) {
                if (r[j + i*n] >= inf &&  r[j + i*n] < sup)
                    n_i[k] += 1;
            }
        }

        for (int t = 0; t < MM; t++)    //compute chi2
            n_i[t] = pow(n_i[t] - (double)n/MM, 2);

        chi2.push_back( accumulate(n_i.begin(), n_i.end(), 0., std::plus<double>())/(n/MM) );
    }

//scrivo i risultati su un file per l'analisi
//migliora il printing con iomanip e facendo una tabella
    print_vector("chi2.dat", chi2, 2);


    return 0;
}

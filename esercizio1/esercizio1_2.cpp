//LSN_Exercise_01
//Exercise 01.2

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <numeric>
#include <iomanip>
#include "Parallel_Generator/random.h"

using namespace std;

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
    const int N_SUM = 10000;
    const int N_MIN = 1000000;      //minima dimensione necessaria per poter fare i 10^4 step su tutte le scelte di N
    std::vector<double> dice_std(N_MIN), dice_exp(N_MIN), dice_lor(N_MIN);

    Random rnd;
    init_random_gen(rnd);

//fill vectors with random numbers
    for (int i = 0; i < N_MIN; i++) {
        dice_std[i] = rnd.Rannyu();
        dice_exp[i] = rnd.Exp(1.);
        dice_lor[i] = rnd.Cauchy(0, 1.);
    }

//set all 4 case in the same for loop
    int N = 0;
    for (int i = 0; i < 4; i++) {
        std::vector<double> sum_std(N_SUM), sum_exp(N_SUM), sum_lor(N_SUM);     //vectors with sums (10^4 elem)

        switch (i) {
            case 0:
                N = 1;
                break;
            case 1:
                N = 2;
                break;
            case 2:
                N = 10;
                break;
            case 3:
                N = 100;
                break;
        }

       for (int j = 0; j < N_SUM; j++) {       //loop on sums
            sum_std[j] = accumulate(dice_std.begin() + j*N, dice_std.begin() + (j+1)*N, 0., std::plus<double>());
            sum_exp[j] = accumulate(dice_exp.begin() + j*N, dice_exp.begin() + (j+1)*N, 0., std::plus<double>());
            sum_lor[j] = accumulate(dice_lor.begin() + j*N, dice_lor.begin() + (j+1)*N, 0., std::plus<double>());
        }

//Print file with results: 1 column for every distribution
        ofstream fileout;
        std::string name_out = "clt" + to_string(N) + ".dat";

        fileout.open(name_out);
        if (fileout.is_open()) {
            for (int j = 0; j < N_SUM; j++)
                fileout << sum_std[j] << " " << sum_exp[j] << " " << sum_lor[j] << endl;
        }
        else {
            cerr << "ERROR: Unable to write file " + name_out << endl;
        }

        fileout.close();
    }


    return 0;
}

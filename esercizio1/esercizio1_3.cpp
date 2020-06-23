//LSN_Exercise_01
//Exercise 01.3 - without using M_PI or other trigonometric functions

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <numeric>
#include <iomanip>
#include <algorithm>
#include "Parallel_Generator/random.h"

using namespace std;


//Function for statistical uncertainty estimation
double error(const std::vector<double>& AV, const std::vector<double>& AV2, int n) {

    return (n == 0) ? 0 : sqrt((AV2[n] - pow(AV[n], 2))/n);
}

//funzione per vedere se l'ago interseca a linea
//prende: posizione centro dell'ago, angolo con la verticale, distanza tra le linee, lunghezza ago
bool hit_line(double x0, double costheta, double d, double L) {
    bool hit = false;

    if (x0 >= d/2.) {                         //sopra metà distanza -> guardo la linea sopra
        if ((L/2.)*costheta >= d - x0)        //distanza dalla linea sopra
            hit = true;
    }
    else {                               //sotto a metà distanza -> guardo la linea sotto
        if ((L/2.)*costheta >= x0)       //distanza dalla linea sotto ... y0 - 0
            hit = true;
    }

    return hit;
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


/***************************************MAIN***********************************/
int main() {
    const double d = 1.;       //distance between lines
    const double L = 0.5;      //length needle
    const int M = 1000000;     //#throws
    const int N = 100;         //#blocks
    const int R = M/N;         //#throws per block

    std::vector<double> x0(M), costheta(M);
    std::vector<double> v_pi(N), v_pi2(N);
    std::vector<double> sum_prog(N), su2_prog(N), err_prog(N);

    Random rnd;
    init_random_gen(rnd);


//inizializzo i vettori di numeri casuali (per x0 = centro dell'ago, theta = angolo con la verticale [0,PI/2])
//uso dei numeri uniformemente distribuito sull'angolo, ma senza usare M_PI
    for (int i = 0; i < M; i++) {
        x0[i] = rnd.Rannyu(0., d);
        double x1, y1;

        do{
            x1 = rnd.Rannyu(-1., 1.);
            y1 = rnd.Rannyu(-1., 1.);
        } while(pow(x1,2) + pow(y1,2) > 1);

        costheta[i] = fabs(x1)/sqrt(pow(x1,2) + pow(y1,2));
    }

//run simulation
    for (int i = 0; i < N; i++) {
        int Nhit = 0;
        int k = 0;

        for (int j = 0; j < R; j++) {
            k = j + i*R;

            if (hit_line(x0[k], costheta[k], d, L))
                Nhit++;
        }

        v_pi[i] = 2.*L*R/(Nhit*d);      //compute pi val in each block
        v_pi2[i] = pow(v_pi[i], 2);
    }

//compute sums and MC error
    for (int i = 0; i < N; i++) {
        sum_prog[i] = accumulate(v_pi.begin(), v_pi.begin()+i+1, 0., std::plus<double>())/(i+1);
        su2_prog[i] = accumulate(v_pi2.begin(), v_pi2.begin()+i+1, 0., std::plus<double>())/(i+1);
        err_prog[i] = error(sum_prog, su2_prog, i);    //Statistical uncertainty
    }

//print results to file
    print_to_file("pi.dat", sum_prog, err_prog);


    return 0;
}

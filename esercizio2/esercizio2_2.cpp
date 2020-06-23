//LSN_Exercise_02
//Exercise 02.2

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <functional>
#include <algorithm>
#include "Parallel_Generator/random.h"

using namespace std;


//Function for statistical uncertainty estimation
//Modified version - take element i and divide always for N
double error_mean(const std::vector<double>& AV, const std::vector<double>& AV2, const int& i, const int& n) {

    return sqrt((AV2[i] - pow(AV[i], 2))/(n - 1.));
}

//RW distance^2 from the origin
double length2_RW(const std::vector<double> & v_RW) {
    double l = 0.;

    for (double elem : v_RW)
        l += pow(elem, 2);

    return l;
}

//compute sign
int random_sgn(Random& ran) {

    return (ran.Rannyu() < 0.5) ? 1 : -1;
}

//Function for computing the next step (starts from the last inserted)
void make_step(std::vector< std::vector<double> >& RW, double a, Random& ran) {

    if (RW.size() == 0) {
        RW.push_back(std::vector<double> (3, 0.));
    }
    else {
        double dir = ran.Rannyu();  //random direction
        int sgn = random_sgn(ran);
        std::vector<double> v_appo(*RW.rbegin());

        if (dir <= 1./3.)       //make step
            v_appo[0] += sgn*a;
        else if (dir >= 2./3.)
            v_appo[1] += sgn*a;
        else
            v_appo[2] += sgn*a;

        RW.push_back(v_appo);
    }
}


//Function for computing the next step (starts from the last inserted), in the continuum
//(Si potrebbe fare anche con un metodo di rigetto per estrarre i punti, campionando i punti all'interno di una cubo e accettando quelli dentro la sfera)
void make_step_continuum (std::vector< std::vector<double> > & RW, double a, Random & ran) {

    if (RW.size() == 0) {
        RW.push_back(std::vector<double> (3, 0.));
    }
    else {
        double phi = ran.Rannyu(0., 2*M_PI); //uniform on a sphere
        double theta = acos(2*ran.Rannyu() - 1);
        std::vector<double> v_appo(*RW.rbegin());

        v_appo[0] += a*sin(theta)*cos(phi);
        v_appo[1] += a*sin(theta)*sin(phi);
        v_appo[2] += a*cos(theta);

        RW.push_back(v_appo);
    }
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


/*-------------------------------------MAIN-----------------------------------*/

int main () {
    const int M = 1000000;  //STEPS
    const int N = 100;      //blocks
    const int L = M/N;      //steps per block
    const int STEPS = 100;  //steps per run
    const double a = 1.;    //length (step)

    Random rnd;
    init_random_gen(rnd);

    std::vector<double> r2(STEPS, 0.), r2_2(STEPS, 0.);     //compute only final values
    std::vector<double> err(STEPS);

//runsimulation M times
//make step, compute length^2 and sum it to respective component of r2
    for (int i = 0; i < N; i++) {
        std::vector<double> sum_l2(STEPS, 0.);   //vector with sums of suqares of legth (for each step)

        for (int j = 0; j < L; j++) {
            std::vector< std::vector<double> > RW;

            for (int k = 0; k < STEPS; k++) {
                make_step(RW, a, rnd);
                sum_l2[k] += length2_RW(RW[k]);
            }
        }

        for (int j = 0; j < STEPS; j++) {   //compute "value" of each block and its square
            r2[j] += sum_l2[j]/L;
            r2_2[j] += pow(sum_l2[j]/L, 2);
        }
    }

//ave and error
    std::transform(r2.begin(), r2.end(), r2.begin(),
            std::bind(std::divides<double>(), std::placeholders::_1, N));
    std::transform(r2_2.begin(), r2_2.end(), r2_2.begin(),
            std::bind(std::divides<double>(), std::placeholders::_1, N));

        for (int i = 0; i < STEPS; i++)
            err[i] = error_mean(r2, r2_2, i, N);


//Print results on file, computing the value of (<r^2>)^1/2
//propagated error: sigma_y = sigma_x/(2sqrt(x))
    ofstream fileout;

    fileout.open("random_walk.dat");
    if (fileout.is_open()) {
        fileout << std::setprecision(5) << std::fixed;
        fileout << 0. << " " << 0. << endl;     //staritng point - stddev doesnt exist ave = 0
        for (int i = 1; i < STEPS; i++)
            fileout << sqrt(r2[i]) << " " << err[i]/(2.*sqrt(r2[i])) << endl;
    }
    else {
        cerr << "ERROR: Unable to write file random_walk.dat" << endl;
    }

    fileout.close();


/*----------------------------------PART 2-----------------------------------*/
//Reset all vectors - just to be sure
    std::fill(r2.begin(), r2.end(), 0.);
    std::fill(r2_2.begin(), r2_2.end(), 0.);
    std::fill(err.begin(),err.end(), 0.);

//run simulation M times
    for (int i = 0; i < N; i++) {
        std::vector<double> sum_l2(STEPS, 0.);   //vector with sums of suqares of legth (for each step)

        for (int j = 0; j < L; j++) {
            std::vector< std::vector<double> > RW;

            for (int k = 0; k < STEPS; k++) {
                make_step_continuum(RW, a, rnd);
                sum_l2[k] += length2_RW(RW[k]);
            }
        }

        for (int j = 0; j < STEPS; j++) {
            r2[j] += sum_l2[j]/L;
            r2_2[j] += pow(sum_l2[j]/L, 2);
        }
    }

//ave and error
    std::transform(r2.begin(), r2.end(), r2.begin(),
            std::bind(std::divides<double>(), std::placeholders::_1, N));
    std::transform(r2_2.begin(), r2_2.end(), r2_2.begin(),
            std::bind(std::divides<double>(), std::placeholders::_1, N));

        for (int i = 0; i < STEPS; i++)
            err[i] = error_mean(r2, r2_2, i, N);


//Print results on file, computing the value of (<r^2>)^1/2
//propagated error: sigma_y = sigma_x/(2sqrt(x))
    fileout.open("random_walk_continuum.dat");
    if (fileout.is_open()) {
        fileout << std::setprecision(5) << std::fixed;
        fileout << 0. << " " << 0. << endl;     //staritng point - stddev doesnt exist ave = 0
        for (int i = 1; i < STEPS; i++)
            fileout << sqrt(r2[i]) << " " << err[i]/(2.*sqrt(r2[i])) << endl;
    }
    else {
        cerr << "ERROR: Unable to write file random_walk.dat" << endl;
    }

    fileout.close();



    return 0;
}

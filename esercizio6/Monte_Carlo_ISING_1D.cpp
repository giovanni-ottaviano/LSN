/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "random.h"
#include "Monte_Carlo_ISING_1D.h"

using namespace std;


int main() {

    Input();  //Inizialization
    int index = 0;

    for (int iblk = 1; iblk <= nblk; iblk++) { //Simulation
        Reset(iblk);   //Reset block averages
        for (int istep = 1; istep <= nstep; istep++) {
            Move(metro);
            Measure();
            Accumulate();   //Update block averages

            index++;
        }
        Averages(iblk);     //Print results for current block
    }
    ConfFinal();            //Write final configuration

    return 0;
}


void Input() {
    ifstream ReadInput;

    cout << "Classic 1D Ising model             " << endl;
    cout << "Monte Carlo simulation             " << endl << endl;
    cout << "Nearest neighbour interaction      " << endl << endl;
    cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
    cout << "The program uses k_B = 1 and mu_B = 1 units " << endl;

//Initialize random generator
    if (restart)
        reinit_random_gen(rnd, true);
    else
        reinit_random_gen(rnd, false);

//Read input informations
    ReadInput.open("input.dat");

    ReadInput >> temp;
    beta = 1./temp;
    cout << "Temperature = " << temp << endl;

    ReadInput >> nspin;
    cout << "Number of spins = " << nspin << endl;

    ReadInput >> J;
    cout << "Exchange interaction = " << J << endl;

    ReadInput >> h;
    cout << "External field = " << h << endl << endl;

    ReadInput >> metro; // if = 1 Metropolis else Gibbs
    ReadInput >> nblk;
    ReadInput >> nstep;
    ReadInput >> restart;

    if (metro == 1)
        cout << "The program performs Metropolis moves" << endl;
    else
        cout << "The program performs Gibbs moves" << endl;

    cout << "Number of blocks = " << nblk << endl;
    cout << "Number of steps in one block = " << nstep << endl;
    cout << "Restart after equilibration = " << ( restart ? "True" : "False" ) << endl << endl;

    ReadInput.close();

//Prepare arrays for measurements
    iu = 0;         //Energy
    ic = 1;         //Heat capacity
    im = 2;         //Magnetization
    ix = 3;         //Magnetic susceptibility

    n_props = 4;    //Number of observables

//Prepare initial configuration
    if (restart) {
        ifstream filein("config.0");

        if (filein.is_open()) {
            for (int i = 0; i < nspin; i++)
                filein >> s[i];
        }
        else {
            cerr << "ERROR: Unable to open file config.0" << endl << endl;
            exit (EXIT_FAILURE);
        }
        filein.close();
    }
    else {
        for (int i = 0; i < nspin; i++)
            s[i] = (rnd.Rannyu() >= 0.5) ? 1. : -1.;
    }

//Evaluate energy etc. of the initial configuration
    Measure();

//Print initial values for the potential energy and virial
    cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
    cout << "Initial heat capacity = " << (walker[ic] - pow(walker[iu], 2))/(pow(temp, 2)*(double)nspin) << endl;
    cout << "Initial magnetization = " << walker[im]/(double)nspin << endl;
    cout << "Initial magnetic sucsptibility = " << walker[ix]/(temp*(double)nspin) << endl << endl;
}


void Move(int metro) {
    int o;
    double p, energy_old, energy_new, sm;
    double energy_up, energy_down, delta_E;

    for (int i = 0 ; i < nspin; ++i) {
        o = (int)(rnd.Rannyu()*nspin);  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)

        if (metro == 1) {   //Metropolis
            energy_old = Boltzmann(s[o], o);
            energy_new = Boltzmann(-1*s[o], o);
            delta_E = energy_new - energy_old;

//flip the spin, if it matches conditions of Metropolis algorithm
            if (delta_E <= 0.) {
                s[o] *= -1;
                accepted += 1.;
            }
            else if (rnd.Rannyu() <= exp(-beta*delta_E)) {
                s[o] *= -1;
                accepted += 1.;
            }
            attempted += 1.;
        }
        else {              //Gibbs sampling
            energy_up = Boltzmann(1, o);
            energy_down = Boltzmann(-1, o);
            p = 1./(1. + exp( -beta*(energy_down - energy_up) ));

            s[o] = (rnd.Rannyu() <= p) ? 1 : -1;    //assign new config

            attempted += 1.;    //the proposed move is always accepted
            accepted += 1.;
        }
    }
}


double Boltzmann (int sm, int ip) {
    double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;

    return ene;
}


void Measure() {
    int bin;
    double u = 0., m = 0.;

//cycle over spins
    for (int i = 0; i < nspin; i++) {
        u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
        m += s[i];
    }

    walker[iu] = u;
    walker[ic] = u*u;
    walker[im] = m;
    walker[ix] = m*m;
}


//Reset block averages
void Reset (int iblk) {

    if (iblk == 1) {
        for (int i = 0; i < n_props; i++) {
            glob_av[i] = 0.;
            glob_av2[i] = 0.;
        }
    }

    for (int i = 0; i < n_props; ++i)
        blk_av[i] = 0;

    blk_norm = 0.;
    attempted = 0.;
    accepted = 0.;
}


//Update block averages
void Accumulate() {

    for (int i = 0; i < n_props; i++)
        blk_av[i] += walker[i];

    blk_norm += 1.;
}


//Print results for current block
void Averages (int iblk) {

    ofstream Ene, Heat, Mag, Chi;
    const int wd = 14;

    cout << "Block number: " << iblk << endl;
    cout << "Acceptance rate: " << accepted/attempted << endl << endl;

//Energy
    Ene.open("output.ene.0", ios::app);
    stima_u = blk_av[iu]/((double)blk_norm*nspin);    //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u = Error(glob_av[iu], glob_av2[iu], iblk);
    Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    Ene.close();

//Heat capacity
    Heat.open("output.heat.0", ios::app);
    stima_c = (blk_av[ic]/blk_norm - pow(blk_av[iu]/blk_norm, 2))/(pow(temp, 2)*(double)nspin);    //Heat capacity
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c = Error(glob_av[ic], glob_av2[ic], iblk);
    Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
    Heat.close();

//Magnetization
    Mag.open("output.mag.0", ios::app);
    stima_m = blk_av[im]/((double)blk_norm*nspin);   //Magnetization
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m = Error(glob_av[im], glob_av2[im], iblk);
    Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    Mag.close();

//Magnetic susceptibility
    Chi.open("output.chi.0", ios::app);
    stima_x = (blk_av[ix]/blk_norm)/(temp*(double)nspin);    //Magnetica susceptibility (implemented for h=0)
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x = Error(glob_av[ix], glob_av2[ix], iblk);
    Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
    Chi.close();

    cout << "----------------------------" << endl << endl;
}

//Write final configuration
void ConfFinal() {
    ofstream WriteConf;

    cout << "Print final configuration to file config.final. " << endl << endl;
    WriteConf.open("config.final");

    for (int i = 0; i < nspin; i++)
        WriteConf << s[i] << endl;

    WriteConf.close();

    rnd.SaveSeed();     //Necessary for restart
}


//Algorithm for periodic boundary conditions
int Pbc (int i) {

    if (i >= nspin)
        i -= nspin;
    else if (i < 0)
        i += nspin;

    return i;
}


double Error (double sum, double sum2, int iblk) {

    return (iblk == 1) ? 0. : sqrt((sum2/(double)iblk - pow(sum/(double)iblk, 2) )/(iblk - 1.));
}

//Function for initializing random generator - if it fails then program is terminated
void reinit_random_gen (Random& gen, bool restart) {
    int seed[4];
    int p1, p2;
    bool success = true;
    std::string seed_name = (restart) ? "seed.out" : "seed.in";

    ifstream Primes("Primes");
    if (Primes.is_open())
       Primes >> p1 >> p2;
    else {
        cerr << "PROBLEM: Unable to open Primes" << endl;
        success = false;
    }
    Primes.close();

    ifstream input(seed_name);
    if (input.is_open()) {
       while ( !input.eof() ) {
           input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
           gen.SetRandom(seed, p1, p2);
       }
       input.close();
    }
    else {
        cerr << "PROBLEM: Unable to open " + seed_name << endl;
        success = false;
    }

    if (!success) {
        cerr << "ERROR: Unable to initialize random generator." << endl;
        cerr << "---STOP---" << endl;
        exit (EXIT_FAILURE);
    }
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

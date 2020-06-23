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
#include "Monte_Carlo_NVT.h"

using namespace std;

int main() {

    Input();  //Inizialization
    int nconf = 1;

    for (int iblk = 1; iblk <= nblk; iblk++) {    //Simulation
        Reset(iblk);   //Reset block averages
        for (int istep = 1; istep <= nstep; istep++) {
            Move();
            Measure();
            Accumulate(); //Update block averages
            Insta_values();

            if (istep%10 == 0) {
//                ConfXYZ(nconf);     //Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
                nconf++;
            }
        }
        Averages(iblk);   //Print results for current block
    }
    ConfFinal(); //Write final configuration

    return 0;
}


void Input() {
    ifstream ReadInput, ReadConf;

    cout << "Classic Lennard-Jones fluid        " << endl;
    cout << "Monte Carlo simulation             " << endl << endl;
    cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
    cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl << endl;
    cout << "The program uses Lennard-Jones units " << endl;

//Read seed for random numbers
    init_random_gen(rnd);

//Read input informations
    ReadInput.open("input.dat");

    ReadInput >> temp;
    beta = 1./temp;
    cout << "Temperature = " << temp << endl;

    ReadInput >> npart;
    cout << "Number of particles = " << npart << endl;

    ReadInput >> rho;
    cout << "Density of particles = " << rho << endl;
    vol = (double)npart/rho;
    box = pow(vol, 1./3.);
    cout << "Volume of the simulation box = " << vol << endl;
    cout << "Edge of the simulation box = " << box << endl;

    ReadInput >> rcut;
    cout << "Cutoff of the interatomic potential = " << rcut << endl << endl;

//Tail corrections for potential energy and pressure
    vtail = (8.*pi*rho)/(9.*pow(rcut, 9)) - (8.*pi*rho)/(3.*pow(rcut, 3));
    ptail = (32.*pi*rho)/(9.*pow(rcut, 9)) - (16.*pi*rho)/(3.*pow(rcut, 3));
    cout << "Tail correction for the potential energy = " << vtail << endl;
    cout << "Tail correction for the virial           = " << ptail << endl;

    ReadInput >> delta;
    ReadInput >> nblk;
    ReadInput >> nstep;

    cout << "The program perform Metropolis moves with uniform translations" << endl;
    cout << "Moves parameter = " << delta << endl;
    cout << "Number of blocks = " << nblk << endl;
    cout << "Number of steps in one block = " << nstep << endl << endl;

    ReadInput.close();

//Prepare arrays for measurements
    iv = 0; //Potential energy
    iw = 1; //Virial

    n_props = 2; //Number of observables

//measurement of g(r)
    igofr = 2;
    nbins = 100;
    n_props += nbins;
    bin_size = (box/2.)/(double)nbins;

//Read initial configuration
    cout << "Read initial configuration from file config.0 " << endl << endl;

    ReadConf.open("config.0");
    for (int i = 0; i < npart; i++) {
        ReadConf >> x[i] >> y[i] >> z[i];
        x[i] = Pbc( x[i] * box );
        y[i] = Pbc( y[i] * box );
        z[i] = Pbc( z[i] * box );
    }

    ReadConf.close();

//Evaluate potential energy and virial of the initial configuration
    Measure();

//Print initial values for the potential energy and virial
    cout << "Initial potential energy (with tail corrections) = " << walker[iv]/(double)npart + vtail << endl;
    cout << "Virial                   (with tail corrections) = " << walker[iw]/(double)npart + ptail << endl;
    cout << "Pressure                 (with tail corrections) = " << rho * temp + (walker[iw] + (double)npart * ptail) / vol << endl << endl;
}


void Move() {
    int o;
    double p, energy_old, energy_new;
    double xold, yold, zold, xnew, ynew, znew;

    for (int i = 0; i < npart; i++) {
        o = (int)(rnd.Rannyu()*npart);  //Select randomly a particle (for C++ syntax, 0 <= o <= npart-1)

//Old
        xold = x[o];
        yold = y[o];
        zold = z[o];

        energy_old = Boltzmann(xold, yold, zold, o);

//New
        xnew = Pbc( x[o] + delta*(rnd.Rannyu() - 0.5) );
        ynew = Pbc( y[o] + delta*(rnd.Rannyu() - 0.5) );
        znew = Pbc( z[o] + delta*(rnd.Rannyu() - 0.5) );

        energy_new = Boltzmann(xnew, ynew, znew, o);

//Metropolis test
        p = exp(beta*(energy_old - energy_new));
        if (p >= rnd.Rannyu()) {    //Update
            x[o] = xnew;
            y[o] = ynew;
            z[o] = znew;

            accepted += 1.;
        }
        attempted += 1.;
    }
}

double Boltzmann (double xx, double yy, double zz, int ip) {
    double ene = 0.;
    double dx, dy, dz, dr;

    for (int i = 0; i < npart; i++) {
        if(i != ip) {   // distance ip-i in pbc
            dx = Pbc(xx - x[i]);
            dy = Pbc(yy - y[i]);
            dz = Pbc(zz - z[i]);
            dr = sqrt(dx*dx + dy*dy + dz*dz);

            if(dr < rcut) {
                ene += 1./pow(dr, 12) - 1./pow(dr, 6);
            }
        }
    }

    return 4.*ene;
}

//Perform a measurement
void Measure() {
    int bin;
    double v = 0., w = 0.;
    double vij, wij;
    double dx, dy, dz, dr;

//reset the hystogram of g(r)
    for (int k = igofr; k < igofr + nbins; k++)
        walker[k] = 0.;

//cycle over pairs of particles
    for (int i = 0; i < npart-1; i++) {
        for (int j = i+1; j < npart; j++) {
            dx = Pbc(x[i] - x[j]);  // distance i-j in pbc
            dy = Pbc(y[i] - y[j]);
            dz = Pbc(z[i] - z[j]);
            dr = sqrt(dx*dx + dy*dy + dz*dz);

            if (dr < rcut) {
                vij = 1./pow(dr, 12) - 1./pow(dr, 6);
                wij = 1./pow(dr, 12) - 0.5/pow(dr, 6);

                v += vij;   // contribution to energy and virial
                w += wij;
            }

//update of the histogram of g(r)
            bin = igofr + (int) (dr/bin_size);  //find correct bin
            walker[bin] += 2.;                  //update histogram
        }
    }

    walker[iv] = 4.*v;
    walker[iw] = 48.*w/3.;
}

//Reset block averages
void Reset (int iblk) {

    if (iblk == 1) {
        for (int i = 0; i < n_props; i++) {
            glob_av[i] = 0;
            glob_av2[i] = 0;
        }
    }

    for (int i = 0; i < n_props; i++)
        blk_av[i] = 0;

   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


//Update block averages
void Accumulate() {

   for (int i = 0; i < n_props; i++)
       blk_av[i] += walker[i];

   blk_norm += 1.;
}


//Print results for current block
void Averages (int iblk) {
    double gofr, deltaV;
    ofstream Gofr, Gave, Epot, Pres;
    const int wd = 12;

    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << (double) accepted/attempted << endl << endl;

    Epot.open("output.epot.0", ios::app);
    Pres.open("output.pres.0", ios::app);
    Gofr.open("output.gofr.0", ios::app);

    stima_pot = blk_av[iv]/blk_norm/(double) npart + vtail;  //Potential energy
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot*stima_pot;
    err_pot = Error(glob_av[iv], glob_av2[iv], iblk);

    stima_pres = rho * temp + (blk_av[iw]/blk_norm + ptail * (double)npart) / vol;  //Pressure
    glob_av[iw] += stima_pres;
    glob_av2[iw] += stima_pres*stima_pres;
    err_press = Error(glob_av[iw], glob_av2[iw], iblk);             //CONTROLLA SE I GLOBAL VANNO DIVISI

//Potential energy per particle
    Epot << setw(wd) << iblk << setw(wd) << stima_pot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_pot << endl;
//Pressure
    Pres << setw(wd) << iblk << setw(wd) << stima_pres << setw(wd) << glob_av[iw]/(double)iblk << setw(wd) << err_press << endl;

//g(r)
    for (int ibin = igofr; ibin < igofr + nbins; ibin++) {
        deltaV = (4.*pi/3.)*(pow(bin_size*(ibin - igofr + 1), 3) - pow(bin_size*(ibin - igofr) ,3));
        gofr = blk_av[ibin]/(rho*m_part*deltaV)/blk_norm;     //compute G(r) (normalize histogram)

        glob_av[ibin] += gofr;
        glob_av2[ibin] += gofr*gofr;

        Gofr << setw(wd) << iblk << setw(wd) << gofr << setw(wd) << endl;
    }

    if (iblk == nblk) {
        Gave.open("output.gave.0");

        for (int ibin = igofr; ibin < igofr + nbins; ibin++) {
            err_gofr = Error(glob_av[ibin], glob_av2[ibin], iblk);

            Gave << setw(wd) << iblk << setw(wd) << glob_av[ibin]/(double)iblk << setw(wd) << err_gofr << endl;
        }

        Gave.close();
    }

    cout << "----------------------------" << endl << endl;

    Epot.close();
    Pres.close();
    Gofr.close();
}


void ConfFinal() {
    ofstream WriteConf;

    cout << "Print final configuration to file config.final " << endl << endl;
    WriteConf.open("config.final");

    for (int i = 0; i < npart; i++)
        WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;

    WriteConf.close();

    rnd.SaveSeed();
}


//Write configuration in .xyz format
void ConfXYZ (int nconf) {
    ofstream WriteXYZ;

    WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
    WriteXYZ << npart << endl;
    WriteXYZ << "This is only a comment!" << endl;

    for (int i = 0; i < npart; i++)
        WriteXYZ << "LJ  " << Pbc(x[i]) << "   " << Pbc(y[i]) << "   " << Pbc(z[i]) << endl;

    WriteXYZ.close();
}


//Algorithm for periodic boundary conditions with side L=box
double Pbc (double r) {

    return r - box * rint(r/box);
}


double Error (double sum, double sum2, int iblk) {

    return (iblk == 1) ? 0. : sqrt((sum2/(double)iblk - pow(sum/(double)iblk, 2))/(iblk - 1.));
}

//Function for initializing random generator - if it fails then program is terminated
void init_random_gen (Random& gen) {
    int seed[4];
    int p1, p2;
    bool success = true;

    ifstream Primes("../Parallel_Generator/Primes");
    if (Primes.is_open())
        Primes >> p1 >> p2;
    else {
        cerr << "PROBLEM: Unable to open Primes" << endl;
        success = false;
    }
    Primes.close();

    ifstream input("../Parallel_Generator/seed.in");
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

//Print instantaneous values of epot and pres
void Insta_values() {
    ofstream insta_epot, insta_P;

    insta_epot.open("output.insta.epot.0", ios::app);
    insta_P.open("output.insta.pres.0", ios::app);

    insta_epot << walker[iv]/(double)npart + vtail << endl;
    insta_P << rho*temp + (walker[iw] + (double)npart*ptail) / vol << endl;

    insta_epot.close();
    insta_P.close();
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

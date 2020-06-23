/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>     // Stream class to both read and write from/to files.
#include <cstdlib>     // srand, rand: to generate random number
#include <cmath>       // rint, pow
#include <vector>
#include <string>
#include <numeric>
#include <iomanip>
#include "MolDyn_NVE.h"

using namespace std;


int main() {

    Input();             //Initialization

    const int steps_bf_mis = 10;                        //Number of steps before a new measure
    int nconf = 1, steps_tot = 0;
    int mis_per_block = nsteps/steps_bf_mis;            //IMPORTANT: choose "steps" multiple of "steps_bf_mis"

    std::vector<double> Epot(nblocks), Epot2(nblocks);  //potential energy for data blocking
    std::vector<double> Ekin(nblocks), Ekin2(nblocks);  //kinetic energy for data blocking
    std::vector<double> Etot(nblocks), Etot2(nblocks);  //total energy for data blocking
    std::vector<double> T(nblocks), T2(nblocks);        //temperature for data blocking
    std::vector<double> P(nblocks), P2(nblocks);        //pressure for data blocking
    std::vector<std::vector<double>> Gofr(nbins, std::vector<double>(nblocks)), Gofr2(nbins, std::vector<double>(nblocks));  //G(r) for data blocking

    std::vector<double> sum_prog_Ep(nblocks), su2_prog_Ep(nblocks), err_prog_Ep(nblocks);   //progressive sums, squares and errors
    std::vector<double> sum_prog_Ek(nblocks), su2_prog_Ek(nblocks), err_prog_Ek(nblocks);
    std::vector<double> sum_prog_Et(nblocks), su2_prog_Et(nblocks), err_prog_Et(nblocks);
    std::vector<double> sum_prog_T(nblocks), su2_prog_T(nblocks), err_prog_T(nblocks);
    std::vector<double> sum_prog_P(nblocks), su2_prog_P(nblocks), err_prog_P(nblocks);
    std::vector<std::vector<double>> sum_prog_G(nbins, std::vector<double>(nblocks));
    std::vector<std::vector<double>> su2_prog_G(nbins, std::vector<double>(nblocks));
    std::vector<std::vector<double>> err_prog_G(nbins, std::vector<double>(nblocks));

//Run the MD simulation
    for (int iblock = 0; iblock < nblocks; iblock++) {
        double sum_ep = 0., sum_ek = 0., sum_et = 0., sum_T = 0., sum_p = 0.;
        std::vector<double> sum_gofr(nbins, 0.);

        for (int j = 0; j < nsteps; j++) {
            Move();           //Move particles with Verlet algorithm
            steps_tot++;      //Add one time step to counter

            if (steps_tot % iprint == 0)
                cout << "Number of time-steps: " << steps_tot << endl;

            if ( (j + 1) % steps_bf_mis == 0) {
                Measure();     //Properties measurement

                sum_ep += stima_pot;
                sum_ek += stima_kin;
                sum_et += stima_etot;
                sum_T += stima_temp;
                sum_p += stima_pres;

                for (int igr = 0; igr < nbins; igr++)
                    sum_gofr[igr] += stima_gofr[igr];

//                ConfXYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
                nconf++;
            }
        }

//compute mean in each block
        Epot[iblock] = sum_ep/(double)mis_per_block;    //each observable is measured "mis_per_block" times in every block
        Ekin[iblock] = sum_ek/(double)mis_per_block;
        Etot[iblock] = sum_et/(double)mis_per_block;
        T[iblock] = sum_T/(double)mis_per_block;
        P[iblock] = sum_p/(double)mis_per_block;
        Epot2[iblock] = pow(Epot[iblock], 2);           //squares
        Ekin2[iblock] = pow(Ekin[iblock], 2);
        Etot2[iblock] = pow(Etot[iblock], 2);
        T2[iblock] = pow(T[iblock], 2);
        P2[iblock] = pow(P[iblock], 2);

        for (int i = 0; i < nbins; i++) {
            Gofr[i][iblock] = sum_gofr[i]/(double)mis_per_block;
            Gofr2[i][iblock] = pow(Gofr[i][iblock], 2);
        }
    }
    ConfFinal();         //Write final configuration to restart


    cout << endl << "End Verlet moves." << endl;
    cout << "Computing averages and errors" << endl;

//Compute averages and uncertainties with blocking method
    for (int i = 0; i < nblocks; i++) {
        sum_prog_Ep[i] = accumulate(Epot.begin(), Epot.begin()+i+1, 0., std::plus<double>())/(i+1);    //potential energy
        su2_prog_Ep[i] = accumulate(Epot2.begin(), Epot2.begin()+i+1, 0., std::plus<double>())/(i+1);
        err_prog_Ep[i] = error(sum_prog_Ep, su2_prog_Ep, i+1);                                           //statistical uncertainty - potential energy

        sum_prog_Ek[i] = accumulate(Ekin.begin(), Ekin.begin()+i+1, 0., std::plus<double>())/(i+1);    //kinetic energy
        su2_prog_Ek[i] = accumulate(Ekin2.begin(), Ekin2.begin()+i+1, 0., std::plus<double>())/(i+1);
        err_prog_Ek[i] = error(sum_prog_Ek, su2_prog_Ek, i+1);                                           //statistical uncertainty - kinetic energy

        sum_prog_Et[i] = accumulate(Etot.begin(), Etot.begin()+i+1, 0., std::plus<double>())/(i+1);    //total energy
        su2_prog_Et[i] = accumulate(Etot2.begin(), Etot2.begin()+i+1, 0., std::plus<double>())/(i+1);
        err_prog_Et[i] = error(sum_prog_Et, su2_prog_Et, i+1);                                           //statistical uncertainty - total energy

        sum_prog_T[i] = accumulate(T.begin(), T.begin()+i+1, 0., std::plus<double>())/(i+1);            //temperature
        su2_prog_T[i] = accumulate(T2.begin(), T2.begin()+i+1, 0., std::plus<double>())/(i+1);
        err_prog_T[i] = error(sum_prog_T, su2_prog_T, i+1);                                             //statistical uncertainty - temperature

        sum_prog_P[i] = accumulate(P.begin(), P.begin()+i+1, 0., std::plus<double>())/(i+1);            //pressure
        su2_prog_P[i] = accumulate(P2.begin(), P2.begin()+i+1, 0., std::plus<double>())/(i+1);
        err_prog_P[i] = error(sum_prog_P, su2_prog_P, i+1);                                             //statistical uncertainty - pressure

        for (int ibin = 0; ibin < nbins; ibin++) {
            sum_prog_G[ibin][i] = accumulate(Gofr[ibin].begin(), Gofr[ibin].begin()+i+1, 0., std::plus<double>())/(i+1);        //G(r)
            su2_prog_G[ibin][i] = accumulate(Gofr2[ibin].begin(), Gofr2[ibin].begin()+i+1, 0., std::plus<double>())/(i+1);
            err_prog_G[ibin][i] = error(sum_prog_G[ibin], su2_prog_G[ibin], i+1);                                         //statistical uncertainty - G(r)
        }
    }

    cout << "---DONE---" << endl;

//Print results to file
    print_to_file("ave_epot.out", sum_prog_Ep, err_prog_Ep);
    print_to_file("ave_ekin.out", sum_prog_Ek, err_prog_Ek);
    print_to_file("ave_etot.out", sum_prog_Et, err_prog_Et, 6);
    print_to_file("ave_temp.out", sum_prog_T, err_prog_T);
    print_to_file("ave_pres.out", sum_prog_P, err_prog_P);

    ofstream Gaveout("ave_gofr.out");
    for (int i = 0; i < nbins; i++)
        Gaveout << i << " " << sum_prog_G[i].back() << " " << err_prog_G[i].back() << endl;

    Gaveout.close();


    return 0;
}

/*---------------------------------END-MAIN-----------------------------------*/


//Prepare all stuff for the simulation
void Input() {
    ifstream ReadInput, ReadConf;

    cout << "Classic Lennard-Jones fluid        " << endl;
    cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
    cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
    cout << "The program uses Lennard-Jones units " << endl;

    seed = 1;    //Set seed for random numbers
    srand(seed); //Initialize random number generator

    ReadInput.open("input.dat"); //Read input

    ReadInput >> temp;

    ReadInput >> npart;
    cout << "Number of particles = " << npart << endl;

    ReadInput >> rho;
    cout << "Density of particles = " << rho << endl;
    vol = (double) npart/rho;
    cout << "Volume of the simulation box = " << vol << endl;
    box = pow(vol, 1./3.);
    cout << "Edge of the simulation box = " << box << endl;

    ReadInput >> rcut;
    ReadInput >> delta;
    ReadInput >> nsteps;
    ReadInput >> iprint;
    ReadInput >> nblocks;
    ReadInput >> restart;

    cout << "The program integrates Newton equations with the Verlet method." << endl;
    cout << "Time step = " << delta << endl;
    cout << "Number of blocks = " << nblocks << endl;
    cout << "Number of steps per block = " << nsteps << endl;
    cout << "Restart from old configuration = " << ( restart ? "True" : "False" ) << endl << endl;
    ReadInput.close();

//Prepare array for measurements
    iv = 0; //Potential energy
    ik = 1; //Kinetic energy
    ie = 2; //Total energy
    it = 3; //Temperature
    ip = 4; //Pressure
    n_props = 5; //Number of observables

//measurement of g(r)
    igofr = 5;
    n_props += nbins;
    bin_size = (box/2.)/(double)nbins;

//Read initial configuration
    if (restart) {
        cout << "Read initial configuration from file old.0 and config.0 " << endl << endl;
        ReadConf.open("config.0");
        for (int i = 0; i < npart; ++i) {
            ReadConf >> x[i] >> y[i] >> z[i];
            x[i] *= box;
            y[i] *= box;
            z[i] *= box;
        }
        ReadConf.close();

        ReadConf.open("old.0");
        if (ReadConf.is_open()) {
            for (int i = 0; i < npart; ++i) {
                ReadConf >> xold[i] >> yold[i] >> zold[i];
                xold[i] *= box;
                yold[i] *= box;
                zold[i] *= box;
            }
        }
        else {
            cerr << "PROBLEM: Unable to open old.0" << endl;
            cerr << "---STOP---" << endl;
            exit (EXIT_FAILURE);
        }
        ReadConf.close();
    }
    else {
        cout << "Read initial configuration from file config.0 " << endl << endl;
        ReadConf.open("config.0");
        for (int i = 0; i < npart; ++i) {
            ReadConf >> x[i] >> y[i] >> z[i];
            x[i] *= box;
            y[i] *= box;
            z[i] *= box;
        }
        ReadConf.close();
    }

//Prepare initial velocities (case: start) or compute them from data (case: restart)
    cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
    double sumv[3] = {0., 0., 0.};
    double sumv2 = 0., appo_temp, fs;

    if (restart) {      //make input config if RESTART -> one step of the Verlet algo and rescale r(t)
        Move();         //one step -> r(t+dt)

        for (int i = 0; i < npart; ++i) {   //compute v(t+dt/2)
            vx[i] = Pbc(x[i] - xold[i])/delta;
            vy[i] = Pbc(y[i] - yold[i])/delta;
            vz[i] = Pbc(z[i] - zold[i])/delta;

            sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
        }
        sumv2 /= (double)npart;
        appo_temp = sumv2/3.;   //T for rescaling velocities
        fs = temp/appo_temp;    //scale factor

        cout << endl << "SCALE FACTOR: " << fs << endl;
    }
    else {      //make input config if START -> extract random veloities and make r(t-dt)
        for (int i = 0; i < npart; i++) {
            vx[i] = rand()/double(RAND_MAX) - 0.5;
            vy[i] = rand()/double(RAND_MAX) - 0.5;
            vz[i] = rand()/double(RAND_MAX) - 0.5;

            sumv[0] += vx[i];
            sumv[1] += vy[i];
            sumv[2] += vz[i];
        }

        for (int idim = 0; idim < 3; idim++)
            sumv[idim] /= (double)npart;

        for (int i = 0; i < npart; i++) {
            vx[i] -= sumv[0];
            vy[i] -= sumv[1];
            vz[i] -= sumv[2];

            sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
        }

        sumv2 /= (double)npart;
        fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
    }

//rescale velocities with appropriate scale factor and estimate new OLD spatial config -> r_new(t) = r(t+dt) - dt*v_s
    for (int i = 0; i < npart; i++) {
        vx[i] *= fs;
        vy[i] *= fs;
        vz[i] *= fs;

        xold[i] = Pbc(x[i] - vx[i]*delta);
        yold[i] = Pbc(y[i] - vy[i]*delta);
        zold[i] = Pbc(z[i] - vz[i]*delta);
    }
}


//Move particles with Verlet algorithm
void Move() {
    double xnew, ynew, znew;
    double fx[m_part], fy[m_part], fz[m_part];

    for (int i = 0; i < npart; i++) { //Force acting on particle i
        fx[i] = Force(i, 0);
        fy[i] = Force(i, 1);
        fz[i] = Force(i, 2);
    }

    for (int i = 0; i < npart; i++) {   //Verlet integration scheme
        xnew = Pbc( 2. * x[i] - xold[i] + fx[i] * pow(delta, 2) );
        ynew = Pbc( 2. * y[i] - yold[i] + fy[i] * pow(delta, 2) );
        znew = Pbc( 2. * z[i] - zold[i] + fz[i] * pow(delta, 2) );

        vx[i] = Pbc(xnew - xold[i])/(2. * delta);
        vy[i] = Pbc(ynew - yold[i])/(2. * delta);
        vz[i] = Pbc(znew - zold[i])/(2. * delta);

        xold[i] = x[i];
        yold[i] = y[i];
        zold[i] = z[i];

        x[i] = xnew;
        y[i] = ynew;
        z[i] = znew;
    }
}


//Compute forces as -Grad_ip V(r)
double Force (int ip, int idir) {
    double f = 0., dr;
    double dvec[3];

    for (int i = 0; i < npart; ++i) {
        if (i != ip) {
            dvec[0] = Pbc( x[ip] - x[i] );  //distance ip-i in pbc
            dvec[1] = Pbc( y[ip] - y[i] );
            dvec[2] = Pbc( z[ip] - z[i] );

            dr = sqrt( dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2] );

            if (dr < rcut)
                f += dvec[idir] * (48./pow(dr,14) - 24./pow(dr,8));     // -Grad_ip V(r)
        }
    }

    return f;
}


//Properties measurement
void Measure() {
    int bin;
    double vij, wij;
    double v = 0., t = 0., w = 0.;
    double dx, dy, dz, dr;
    double deltaV;
    ofstream Epot, Ekin, Etot, Temp, Pres, Gofr;

    Epot.open("output_epot.dat", ios::app);
    Ekin.open("output_ekin.dat", ios::app);
    Temp.open("output_temp.dat", ios::app);
    Etot.open("output_etot.dat", ios::app);
    Pres.open("output_pres.dat", ios::app);
    Gofr.open("output_gofr.dat", ios::app);

    for (int igr = 0; igr < nbins; igr++)
        stima_gofr[igr] = 0.;


//cycle over pairs of particles
    for (int i = 0; i < npart-1; i++) {
        for (int j = i+1; j < npart; j++) {
            dx = Pbc( xold[i] - xold[j] );   // here I use old configurations [old = r(t)]
            dy = Pbc( yold[i] - yold[j] );  // to be compatible with EKin which uses v(t)
            dz = Pbc( zold[i] - zold[j] );  // => Epot should be computed with r(t)
            dr = sqrt( dx*dx + dy*dy + dz*dz );

            if (dr < rcut) {
                vij = 1./pow(dr, 12) - 1./pow(dr, 6);
                wij = 1./pow(dr, 12) - 0.5/pow(dr, 6);

                v += vij;   //Potential energy
                w += wij;   //Virial
            }

//update of the histogram of g(r)
            bin = (int) (dr/bin_size);  //find correct bin
            stima_gofr[bin] += 2.;      //update histogram
        }
    }

    v *= 4.;
    w *= 48./3.;

//Kinetic energy
    for (int i = 0; i < npart; ++i)
        t += 0.5*(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

    stima_pot = v/(double)npart;            //Potential energy per particle
    stima_kin = t/(double)npart;            //Kinetic energy per particle
    stima_temp = (2./3.) * t/(double)npart; //Temperature
    stima_etot = (t + v)/(double)npart;     //Total energy per particle
    stima_pres = rho*stima_temp + w/vol;

    for (int ibin = 0; ibin < nbins; ibin++) {
        deltaV = (4.*M_PI/3.)*(pow(bin_size*(ibin + 1), 3) - pow(bin_size*(ibin), 3));
        stima_gofr[ibin] *= 1./(rho*m_part*deltaV);
    }

    blk++;  //block number (useful for gofr)

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;
    Pres << stima_pres << endl;

    for (int ibin = 0; ibin < nbins; ibin++) {
        Gofr << blk << " " << stima_gofr[ibin] << endl;
    }

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Pres.close();
    Gofr.close();
}


//Write final configuration
void ConfFinal() {
    ofstream WriteConf;

    cout << "Print final configuration to file config.final " << endl;
    WriteConf.open("config.final");

    for (int i = 0; i < npart; ++i)
        WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;

    WriteConf.close();

    cout << "Print old final configuration to file old.final " << endl;
    WriteConf.open("old.final");

    for (int i = 0; i < npart; ++i)
        WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;

    WriteConf.close();
    cout << endl;
}


//Write configuration in .xyz format
void ConfXYZ (int nconf) {
    ofstream WriteXYZ;

    WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
    WriteXYZ << npart << endl;
    WriteXYZ << "This is only a comment!" << endl;

    for (int i = 0; i < npart; i++)
        WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;

    WriteXYZ.close();
}


//Algorithm for periodic boundary conditions with side L=box
double Pbc (double r) {

    return r - box*rint(r/box);
}


//Function for statistical uncertainty estimation
double error (const std::vector<double> & AV, const std::vector<double> & AV2, const int & n) {

    return (n == 1) ? 0. : sqrt((AV2[n-1] - pow(AV[n-1], 2))/(n - 1.));
}

//Function for printing results to file
template <typename T, typename U>
void print_to_file (const std::string & name, const std::vector<T> & v1, const std::vector<U> & v2, int prec) {
    ofstream fileout;
    int size = std::min(v1.size(), v2.size());

    if (v1.size() != v2.size())
        cout << "WARNING: function print_to_file received vectors with different sizes." << endl;

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

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

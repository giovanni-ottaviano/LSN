/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __NVT_h_
#define __NVT_h_

#include "random.h"

//Random numbers
Random rnd;

//parameters, observables
const int m_props = 1000;
int n_props, iv, iw, igofr;
double vtail, ptail, bin_size, nbins, sd;
double walker[m_props];

//averages
double blk_av[m_props], blk_norm;
int accepted, attempted;
double glob_av[m_props], glob_av2[m_props];
double stima_pot, stima_pres, err_pot, err_press, err_gofr;

//configuration
const int m_part = 108;
double x[m_part], y[m_part], z[m_part];

//thermodynamical state
int npart;
double beta, temp, vol, rho, box, rcut;

//simulation
int nstep, nblk;
double delta;

//pigreco
const double pi = 3.1415927;

//functions
void Input();
void Reset(int);
void Accumulate();
void Averages(int);
void Move();
void ConfFinal();
void ConfXYZ(int);
void Measure();
double Boltzmann(double, double, double, int);
double Pbc(double);
double Error(double, double, int);
void init_random_gen (Random&);
void Insta_values();

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

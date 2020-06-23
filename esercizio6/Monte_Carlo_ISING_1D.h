/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#ifndef _Ising_1D_h_
#define _Ising_1D_h_

#include "random.h"

//Random numbers
Random rnd;

//parameters, observables
const int m_props = 1000;
int n_props, iu, ic, im, ix, ig;
double nbins;
double walker[m_props];

//averages
double blk_av[m_props], blk_norm, accepted, attempted;
double glob_av[m_props], glob_av2[m_props];
double stima_u, stima_c, stima_m, stima_x, stima_g;
double err_u, err_c, err_m, err_x, err_g;

//configuration
const int m_spin = 50;
double s[m_spin];

//thermodynamical state
int nspin;
double beta, temp, J, h;

//simulation parameters
int nstep, nblk, metro, restart;

//functions
void Input();
void Reset(int);
void Accumulate();
void Averages(int);
void Move(int);
void ConfFinal();
void Measure();
double Boltzmann(int, int);
int Pbc(int);
double Error(double, double, int);
void reinit_random_gen(Random&, bool);

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

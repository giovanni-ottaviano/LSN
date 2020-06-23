/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#ifndef _MolDyn_h_
#define _MolDyn_h_

//parameters, observables
const int m_props = 4;
const int nbins = 100;
int n_props;
int iv, ik, it, ie, ip, igofr;
double stima_pot, stima_kin, stima_etot, stima_temp, stima_pres, stima_gofr[nbins];
double bin_size;

// averages
double acc, att;

//configuration
const int m_part = 108;
double x[m_part], y[m_part], z[m_part], xold[m_part], yold[m_part], zold[m_part];
double vx[m_part], vy[m_part], vz[m_part];

// thermodynamical state
int npart;
double energy, temp, vol, rho, box, rcut;

// simulation
int nsteps, nblocks, blk = 0;
int iprint, seed, restart;
double delta;


//functions
void Input();
void Move();
void ConfFinal();
void ConfXYZ(int);
void Measure();
double Force(int, int);
double Pbc(double);
double error(const std::vector<double>& , const std::vector<double>&, const int&);
template <typename T, typename U>
void print_to_file(const std::string&, const std::vector<T>&, const std::vector<U>&, int prec = 5);

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

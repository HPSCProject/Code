//  Edited and expanded by PR, Aug 2013

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <locale>
#include <sstream>
#include <string>
#include <vector>
#include "functions.h"
// #include <omp.h>



using namespace std;

int main(int argc, const char* argv[])
{
    double write_start, write_time(0), eq_start, eq_time(0);
    double main_start = get_walltime();
  int nx = 251;
  int ny = 251;
  int nsteps = 250;  // are these time steps??
  double sd = 50.0;   // standard deviation
  double omega = 1.0; // ??

  double dtreal = 12;                  // micro-seconds
  double wxreal = 0.000785398;         // micro-Hz
  int dtsim = round(dtreal*wxreal*sd); // round result to nearest int
  double lambda = 1.0;

  cout << "dtsim = " << dtsim << endl;

  int t_off = 10;           // time step to turn potential off
  int t_on = t_off + dtsim; // time step to turn it back on

  bool ftrue = true; // true if potential on, false if potential off

  double u_sqr, c_dot_u, force;
  int in, jn;

  double T0 = 1.0/3.0;

  vector<vector<vector<double> > > fIn (nx, vector<vector<double> >(ny, vector<double>(Q, 0.0)));

  vector<vector<vector<double> > > fOut(nx, vector<vector<double> >(ny, vector<double>(Q, 0.0)));

  vector<vector<double> > rho (nx, vector<double>(ny, 1.0));
  vector<vector<double> > ux (nx, vector<double>(ny, 0.0));
  vector<vector<double> > uy (nx, vector<double>(ny, 0.0));

  int c[Q][D] = {{0, 0}, {1, 0}, {0, 1}, {-1, 0}, {0, -1}, {1, 1}, {-1, 1}, {-1, -1} ,{1, -1}};

  int nop[Q] = {0, 3, 4, 1, 2, 7, 8, 5, 6};

  double w1 = 4.0 / 9.0;
  double w2 = 1.0 / 9.0;
  double w3 = 1.0 / 36.0;
  double wi[] = {w1, w2, w2, w2, w2, w3, w3, w3, w3};
    double gaussian_start = get_walltime();
  init_gaussian(&fIn, &fOut, &rho, &ux, &uy, c, wi, lambda, nx, ny, sd, T0, omega);
    double gaussian_time = get_walltime() - gaussian_start;
  for (int ts = 0; ts < nsteps; ++ts)
    {
      if (ts == t_on)
	{
	  ftrue = true;
	}
      else if (ts == t_off)
	{
	  ftrue = false;
	}
        eq_start = get_walltime();
      eq_and_stream(&fIn, &fOut, &rho, &ux, &uy, c, wi, nop, lambda, nx, ny, T0, omega, sd, ftrue);
        eq_time += get_walltime() - eq_start;
      fIn = fOut;
      if(ts%10==0)
          write_start = get_walltime();
      write_gaussian(&rho, &ux, &uy, nx, ny, sd, ts);
        write_time += get_walltime() - write_start;
    }

    double main_time = get_walltime() - main_start;
    
    cout << "For size " << NX << "x" << NY << " over " << N_STEPS << " Timesteps" << endl;
    cout << "\tTotal Breathing Time: \t\t" << main_time << endl;
    cout << "\tTotal Time For init_gaussian(): \t" << gaussian_time << endl;
    cout << "\tTotal Time For eq_and_stream(): \t" << eq_time << endl;
    cout << "\tTotal Time For write_gaussian(): \t" << write_time << endl;
    
  return 0;
}

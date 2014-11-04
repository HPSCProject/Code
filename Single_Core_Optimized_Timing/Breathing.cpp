//  Edited and expanded by PR, Aug 2013
#include "functions.h"

using namespace std;

int main(int argc, const char* argv[])
{
  double write_start(0), write_time (0);
    double main_start = get_walltime();
  const int c[Q][D] = {{0, 0}, {1, 0}, {0, 1}, {-1, 0}, {0, -1}, {1, 1}, {-1, 1}, {-1, -1}, {1, -1}};
  const int nop[Q] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
  const double wi[9] = {W1, W2, W2, W2, W2, W3, W3, W3, W3};

  bool ftrue = true;	// true indicates trapping potential is turned on; false indicates trapping potential is off.

  //cout << "dtsim = " << DT_SIM << endl;

  // Input & output arrays.
  qr_type ***fIn, ***fOut;

  // Fill input & output arrays w/ 0.
  fIn = new qr_type**[NX];
  fOut = new qr_type**[NX];
  for(int i = 0; i < NX; ++i)
  {
    fIn[i] = new qr_type*[NY];
    fOut[i] = new qr_type*[NY];
    for(int j = 0; j < NY; ++j)
    {
      fIn[i][j] = new qr_type[Q];
      memset(fIn[i][j], qr_type(0.0), F_SIZE);
      fOut[i][j] = new qr_type[Q];
      memset(fOut[i][j], qr_type(0.0), F_SIZE);
    }
  }

  // Declare the macroscopic variables for the fluid
  qr_type **rho;  // Fluid density
  qr_type **ux;	  // Fluid velocity in the x direction
  qr_type **uy;	  // Fluid velocity in the y direction

  // Fill rho, ux, uy.
  rho = new qr_type*[NX];
  ux = new qr_type*[NX];
  uy = new qr_type*[NX];
  for(int i = 0; i < NX; ++i)
  {
    rho[i] = new qr_type[NY];
    memset(rho[i], qr_type(0.0), TWOD_SIZE);
    ux[i] = new qr_type[NY];
    memset(ux[i], qr_type(0.0), TWOD_SIZE);
    uy[i] = new qr_type[NY];
    memset(uy[i], qr_type(0.0), TWOD_SIZE);
  }
    double gaussian_start = get_walltime();
  init_gaussian(fIn, fOut, wi);
    double gaussian_time = get_walltime() - gaussian_start;
    double eq_time = 0;
    
  for (int ts = 0; ts < N_STEPS; ++ts)
  {
    if (ts == T_ON)
    {
      ftrue = true;
    }
    else if (ts == T_OFF)
    {
      ftrue = false;
    }
      double eq_start = get_walltime();
    eq_and_stream(fIn, fOut, rho, ux, uy, c, wi, nop, ftrue);
      double eq_iter_time = get_walltime() - eq_start;
      eq_time += eq_iter_time;
      
    // fIn = fOut
    for(int i = 0; i < NX; ++i)
    {
      for(int j = 0; j < NY; ++j)
      {
        memcpy(fIn[i][j], fOut[i][j], F_SIZE);
      }
    }
    write_start = get_walltime(); 
    if (ts % 10 == 0)
    {
      write_gaussian(rho, ux, uy, ts);
    }
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

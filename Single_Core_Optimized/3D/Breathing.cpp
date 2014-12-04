#include "functions.h"

using namespace std;

int main(int argc, const char * argv[])
{
  const double c[Q][D] = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, -1.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 0.0, -1.0},
                           {1.0, 1.0, 0.0}, {-1.0, -1.0, 0.0}, {1.0, -1.0, 0.0}, {-1.0, 1.0, 0.0}, {1.0, 0.0, 1.0}, {-1.0, 0.0, -1.0}, {1.0, 0.0, -1.0},
                           {-1.0, 0.0, 1.0}, {0.0, 1.0, 1.0}, {0.0, -1.0, -1.0}, {0.0, 1.0, -1.0}, {0.0, -1.0, 1.0}, {1.0, 1.0, 1.0}, {1.0, 1.0, -1.0},
                           {1.0, -1.0, 1.0}, {-1.0, 1.0, 1.0}, {-1.0, -1.0, 1.0}, {-1.0, 1.0, -1.0}, {1.0, -1.0, -1.0}, {-1.0, -1.0, -1.0}, {3.0, 0.0, 0.0},
                           {-3.0, 0.0, 0.0}, {0.0, 3.0, 0.0}, {0.0, -3.0, 0.0}, {0.0, 0.0, 3.0}, {0.0, 0.0, -3.0}, {3.0, 3.0, 3.0}, {3.0, 3.0, -3.0},
                           {3.0, -3.0, 3.0}, {-3.0, 3.0, 3.0}, {-3.0, -3.0, 3.0}, {-3.0, 3.0, -3.0}, {3.0, -3.0, -3.0}, {-3.0, -3.0, -3.0}};
  const double wi[Q] = {W0, W1, W1, W1, W1, W1, W1, W2, W2, W2, W2, W2, W2, W2, W2, W2, W2, W2, W2, W3, W3, W3, W3, W3, W3, W3, W3, W4, W4, W4, W4, W4, W4,
                         W5, W5, W5, W5, W5, W5, W5, W5};
  bool ftrue = true;

  // Input & output arrays.
  qr_type**** fIn, ****fOut;

  // Fill arrays w/ 0.
  fIn = new qr_type***[NX];
  fOut = new qr_type***[NX];
  for(int i = 0; i < NX; ++i)
  {
    fIn[i] = new qr_type**[NY];
    fOut[i] = new qr_type**[NY];
    for(int j = 0; j < NY; ++j)
    {
      fIn[i][j] = new qr_type*[NZ];
      fOut[i][j] = new  qr_type*[NZ];
      for(int k = 0; k < NZ; ++k)
      {
        fIn[i][j][k] = new qr_type[Q];
        fOut[i][j][k] = new qr_type[Q];
        memset(fIn[i][j][k], qr_type(0.0), F_SIZE);
        memset(fOut[i][j][k], qr_type(0.0), F_SIZE);
      }
    }
  }

  // Declare the macroscopic variables for the fluid.
  qr_type*** rho, ***ux, ***uy, ***uz;

  // Fill arrays w/ 0.
  rho = new qr_type**[NX];
  ux = new qr_type**[NX];
  uy = new qr_type**[NX];
  uz = new qr_type**[NX];
  for(int i = 0; i < NX; ++i)
  {
    rho[i] = new qr_type*[NY];
    ux[i] = new qr_type*[NY];
    uy[i] = new qr_type*[NY];
    uz[i] = new qr_type*[NY];
    for(int j = 0; j < NY; ++j)
    {
      rho[i][j] = new qr_type[NZ];
      ux[i][j] = new qr_type[NZ];
      uy[i][j] = new qr_type[NZ];
      uz[i][j] = new qr_type[NZ];
      memset(rho[i][j], qr_type(0.0), D_SIZE);
      memset(ux[i][j], qr_type(0.0), D_SIZE);
      memset(uy[i][j], qr_type(0.0), D_SIZE);
      memset(uz[i][j], qr_type(0.0), D_SIZE);
    }
  }
     
  make_lattice(fIn, fOut, c, wi);

  // Time it.
  double start = 0.0;
  double end = 0.0;

  get_walltime(start);
    
  // Perform steps.
  for(int ts = 0; ts <= STEPS; ++ts)
  {
    eq(fIn, fOut, rho, ux, uy, uz, c, wi, ftrue);

    if (ts % MEAS_STEPS == 0)
    {
	    write_gaussian(rho, ux, ts);
    }
    
    stream(fIn, fOut, c);
  }

  get_walltime(end);

  printf("Problem Size: NX=%d, NY=%d, NZ=%d, SD=%f, N_STEPS=%d\n", NX, NY, NZ, SD, STEPS);
  printf("Walltime: %fs\n", end - start);

  // Free allocated memory
  for(int i = 0; i < NX; ++i)
  {
    for(int j = 0; j < NY; ++j)
    {
      delete[] rho[i][j];
      delete[] ux[i][j];
      delete[] uy[i][j];
      delete[] uz[i][j];
    }
    delete[] rho[i];
    delete[] ux[i];
    delete[] uy[i];
    delete[] uz[i];
  }
  delete[] rho;
  delete[] ux;
  delete[] uy;
  delete[] uz;

  for(int i = 0; i < NX; ++i)
  {
    for(int j = 0; j < NY; ++j)
    {
      for(int k = 0; k < NZ; ++k)
      {
        delete[] fIn[i][j][k];
        delete[] fOut[i][j][k];
      }
      delete[] fIn[i][j];
      delete[] fOut[i][j];
    }
    delete[] fIn[i];
    delete[] fOut[i];
  }
  delete[] fIn;
  delete[] fOut;

  return 0;
}

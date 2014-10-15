//  Edited and expanded by PR, Aug 2013

#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include <locale>
#include <sstream>
#include <string>
#include "functions.h"

using namespace std;

int main(int argc, const char* argv[])
{
  const int c[Q][D] = {{0, 0}, {1, 0}, {0, 1}, {-1, 0}, {0, -1}, {1, 1}, {-1, 1}, {-1, -1}, {1, -1}};
  const int nop[Q] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
  const double wi[9] = {W1, W2, W2, W2, W2, W3, W3, W3, W3};

  bool ftrue = true;	// Description

  cout << "dtsim = " << DT_SIM << endl;

  // Input & output arrays.
  qr_type fIn[NX*NY*Q], fOut[NX*NY*Q];

  // Fill input & output arrays w/ 0.
  memcpy(fIn, 0.0, F_SIZE);
  memcpy(fOut, 0.0, F_SIZE);

  // What is this stuff?
  qr_type rho[NX*NY];	// Description
  qr_type ux[NX*NY];	// Description
  qr_type uy[NX*NY];	// Description

  // Fill rho, ux, uy.
  memcpy(rho, 1.0, 2D_SIZE);
  memcpy(ux, 0.0, 2D_SIZE);
  memcpy(uy, 0.0, 2D_SIZE);

  init_gaussian(fIn, fOut, rho, ux, uy, wi);

  for (size_t ts = 0; ts < N_STEPS; ++ts)
  {
    if (ts == T_ON)
    {
      ftrue = true;
    }
    else if (ts == T_OFF)
    {
      ftrue = false;
    }

    eq_and_stream(fIn, rho, ux, uy, c, wi, ftrue);

    fIn = fOut;

    write_gaussian(rho, ux, uy, ts);
  }

  return 0;
}

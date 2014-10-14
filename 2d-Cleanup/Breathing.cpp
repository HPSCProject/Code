//  Edited and expanded by PR, Aug 2013

#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
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

  // Input & output array.
  double fIn[NX*NY*Q], fOut[NX*NY*Q];

  // Fill input & output arrays w/ 0.
  memcpy(fIn, 0.0, sizeof(double)*NX*NY*Q);
  memcpy(fOut, 0.0, sizeof(double)*NX*NY*Q);

  // What is this stuff?
  double rho[NX*NY];	// Description
  double ux[NX*NY];		// Description
  double uy[NX*NY];		// Description

  // Fill rho, ux, uy.
  memcpy(rho, 1.0, sizeof(double)*NX*NY);
  memcpy(ux, 0.0, sizeof(double)*NX*NY);
  memcpy(uy, 0.0, sizeof(double)*NX*NY);

  init_gaussian(fIn, fOut, rho, ux, uy, c, wi, LAMBDA, NX, NY, SD, T0, OMEGA);

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

    eq_and_stream(fIn, fOut, rho, ux, uy, c, wi, nop, LAMBDA, NX, NY, T0, OMEGA, SD, ftrue);

    fIn = fOut;

    write_gaussian(rho, ux, uy, NX, NY, SD, ts);
  }

  return 0;
}

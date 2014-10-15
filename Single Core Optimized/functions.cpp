#include "functions.h"
using namespace std;

void init_gaussian(qr_type fIn[NX*NY*Q], qr_type fOut[NX*NY*Q], qr_type rho[NX*NY], qr_type ux[NX*NY], qr_type uy[NX*NY], const double wi[Q])
{
  qr_type fEq;

  for (size_t i = 0; i < NX; ++i)
  {
    for (size_t j = 0; j < NY; ++j)
    {
      rho[i*NY+j] = exp(-(0.5 / T0) * ((qr_type(i) - MIDDLE_X) / SD) * ((qr_type(i) - MIDDLE_X) / SD) - (0.5 / T0) * LAMBDA
                    * LAMBDA * ((qr_type(j) - MIDDLE_Y) / SD) * ((qr_type(j) - MIDDLE_Y) / SD));

      for (size_t n = 0; n < Q; ++n)
      {
      	fEq = wi[n] * rho[i*NY+j];

      	fIn[i*NY*Q+j*Q+n] = fEq;
      	fOut[i*NY*Q+j*Q+n] = fEq;
      }
    }
  }
}

void eq_and_stream(qr_type fIn[NX*NY*Q], qr_type rho[NX*NY], qr_type ux[NX*NY], qr_type uy[NX*NY], const int c[Q][D], const double wi[Q], const bool& ftrue)
{
  bool notify = true;
  qr_type u_sqr, c_dot_u, force;
  qr_type x, y, temp, fEq;

  for (size_t i = 0; i < NX; ++i)
  {
    for (size_t j = 0; j < NY; ++j)
    {
      x = (qr_type(i) - MIDDLE_X) / SD;
      y = LAMBDA * LAMBDA * (qr_type(j) - MIDDLE_Y) / SD;

      // Set to zero before summing.
      rho[i*NY+j] = 0.0;
      ux[i*NY+j] = 0.0;
      uy[i*NY+j] = 0.0;

      for (size_t n = 0; n < Q; ++n)
	  {
	  	temp = fIn[i*NY*Q+j*Q+n];
        rho[i*NY+j] += temp;
	    ux[i*NY+j] += c[n][0] * temp;
	    uy[i*NY+j] += c[n][1] * temp;
	  }

      ux[i*NY+j] /= rho[i*NY+j];
      uy[i*NY+j] /= rho[i*NY+j];

      u_sqr = (ux[i*NY+j] * ux[i*NY+j]) + (uy[i*NY+j] * uy[i*NY+j]);

      for (size_t n = 0; n < Q; ++n)
      {
      	c_dot_u = (c[n][0] * ux[i*NY+j]) + (c[n][1] * uy[i*NY+j]);

      	fEq = wi[n] * rho[i*NY+j] * (1.0 + (1.0 / T0) * c_dot_u + (1.0 / (2.0 * T0 * T0)) * c_dot_u * c_dot_u - (1.0 / (2.0 * T0)) * u_sqr);

      	if (ftrue)
        {
      	  force = -(1.0 / T0) * ((c[n][0] * x) + (c[n][1] * y) - (ux[i*NY+j] * x) + (uy[i*NY+j] * y)) * fEq / SD;
      	  if (notify)
          {
      	    cout << "potential is on" << endl;
      	    notify = false;
      	  }
      	}
      	else
        {
      	  force = 0.0;
      	  if (notify)
          {
      	    cout << "potential is off" << endl;
      	    notify = false;
      	  }
      	}

	    fIn[i*NY*Q+j*Q+n] *= (1.0 - OMEGA)) + OMEGA * fEq + force;
      }
    }
  }
}

void write_gaussian(qr_type rho[NX*NY], qr_type ux[NX*NY], qr_type uy[NX*NY], const size_t& ts)
{
  fstream out;
  char fname[255];

  sprintf(fname, "data/Xrho_t%i.dat", ts);
  out.open(fname, ios::out);
  //out.open("data.txt");
  for(size_t i = 0; i < NX; ++i)
  {    
    out << (i - MIDDLE_X) * SIN_V << "\t";
    out << rho[i*NY+FLOOR_Y] << "\n";
  }

  out.close();

  sprintf(fname, "data/Yrho_t%i.dat", ts);
  out.open(fname, ios::out);
  for(size_t j = 0; j < NY; ++j)
  {
    out << (j - MIDDLE_Y) * SIN_V << "\t";
    out << rho[FLOOR_X*NY+j] << "\n";
  }

  out.close();

  sprintf(fname, "data/Xux_t%i.dat", ts);
  out.open(fname, ios::out);
  for(size_t i = 0; i < NX; ++i)
  {
    out << (i - MIDDLE_X) * SIN_V << "\t";
    out << ux[i*NY+FLOOR_Y] << "\n";
  }

  out.close();

  sprintf(fname, "data/Yuy_t%i.dat", ts);
  out.open(fname, ios::out);
  for(size_t j = 0; j < NY; ++j)
  {
    out << (j - MIDDLE_Y) * SIN_V << "\t";
    out << uy[FLOOR_X*NY+j] << "\n";
  }

  out.close();
}

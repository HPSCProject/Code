#include "functions.h"
using namespace std;

void init_gaussian(qr_type ***fIn, qr_type ***fOut, const double wi[Q])
{
  qr_type fEq;

  for (int i = 0; i < NX; ++i)
  {
    for (int j = 0; j < NY; ++j)
    {
      for (int n = 0; n < Q; ++n)
      {
        fEq = wi[n] * exp(-(0.5 / T0) * ((i - MIDDLE_X) / SD) * ((i - MIDDLE_X) / SD) - (0.5 / T0) * LAMBDA
                    * LAMBDA * ((j - MIDDLE_Y) / SD) * ((j - MIDDLE_Y) / SD));
      	
        fIn[i][j][n] = fEq;
      	fOut[i][j][n] = fEq;
      }
    }
  }
}

void eq_and_stream(qr_type ***fIn, qr_type ***fOut, qr_type **rho, qr_type **ux, qr_type **uy, const int c[Q][D], const double wi[Q], const int nop[Q], const bool& ftrue)
{
  bool notify = true;
  qr_type u_sqr, c_dot_u, force;
  qr_type x, y, temp, fEq;

  for (int i = 0; i < NX; ++i)
  {
    for (int j = 0; j < NY; ++j)
    {
      // Set to zero before summing
      rho[i][j] = qr_type(0.0);
      ux[i][j] = qr_type(0.0);
      uy[i][j] = qr_type(0.0);

      x = (i - MIDDLE_X) / SD;
      y = LAMBDA * LAMBDA * (j - MIDDLE_Y) / SD;

      for (int n = 0; n < Q; ++n)
	    {
	  	  temp = fIn[i][j][n];
        rho[i][j] += temp;
	      ux[i][j] += c[n][0] * temp;
	      uy[i][j] += c[n][1] * temp;
	    }

      ux[i][j] /= rho[i][j];
      uy[i][j] /= rho[i][j];

      u_sqr = (ux[i][j] * ux[i][j]) + (uy[i][j] * uy[i][j]);

      for (int n = 0; n < Q; ++n)
      {
      	c_dot_u = (c[n][0] * ux[i][j]) + (c[n][1] * uy[i][j]);

      	fEq = wi[n] * rho[i][j] * (1.0 + (1.0 / T0) * c_dot_u + (1.0 / (2.0 * T0 * T0)) * c_dot_u * c_dot_u - (1.0 / (2.0 * T0)) * u_sqr);

      	if (ftrue)
        {
      	  force = -(1.0 / T0) * (((c[n][0] * x) + (c[n][1] * y)) - ((ux[i][j] * x) + (uy[i][j] * y))) * fEq / SD;
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

	    fIn[i][j][n] = fIn[i][j][n] * (1.0 - OMEGA) + OMEGA * fEq + force;
      }
    }
  }

  int in, jn;

  for (int i = 0; i < NX; ++i)
  {
    for (int j = 0; j < NY; ++j)
    {
      for (int n = 0; n < Q; ++n)
      {
        in = i + int(c[n][0]);
        jn = j + int(c[n][1]);

        if (in > NX - 1 || in < 0)
        {
          in = (in + NX) % NX;                                                                                                                                                                                        
        }
        if (jn > NY-1 || jn < 0)
        {
          jn = (jn + NY) % NY;                                                                                                                                                                                        
        }

        fOut[i][j][nop[n]] = fIn[in][jn][nop[n]];
      }
    }
  }
}

void write_gaussian(qr_type **rho, qr_type **ux, qr_type **uy, const int& ts)
{
  fstream out;
  char fname[255];
  
  sprintf(fname, "data/Xrho_t%i.dat", ts);
  out.open(fname, ios::out);
  for(int i = 0; i < NX; ++i)
  {    
    out << (i - MIDDLE_X) * SIN_V << "\t";
    out << rho[i][MIDDLE_Y] << "\n";
  }

  out.close();

  sprintf(fname, "data/Yrho_t%i.dat", ts);
  out.open(fname, ios::out);
  for(int j = 0; j < NY; ++j)
  {
    out << (j - MIDDLE_Y) * SIN_V << "\t";
    out << rho[MIDDLE_X][j] << "\n";
  }

  out.close();

  sprintf(fname, "data/Xux_t%i.dat", ts);
  out.open(fname, ios::out);
  for(int i = 0; i < NX; ++i)
  {
    out << (i - MIDDLE_X) * SIN_V << "\t";
    out << ux[i][MIDDLE_Y] << "\n";
  }

  out.close();

  sprintf(fname, "data/Yuy_t%i.dat", ts);
  out.open(fname, ios::out);
  for(int j = 0; j < NY; ++j)
  {
    out << (j - MIDDLE_Y) * SIN_V << "\t";
    out << uy[MIDDLE_X][j] << "\n";
  }

  out.close();
}

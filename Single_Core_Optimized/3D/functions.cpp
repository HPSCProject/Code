#include "functions.h"
using namespace std;

void make_lattice(qr_type**** fIn, qr_type**** fOut, const double c[Q][D], const double wi[Q])
{
  double c_sqr, rho;

  for (int i = 0; i < NX ; ++i)
  {
    for (int j = 0; j < NY ; ++j)
	  {
  	  for (int k = 0; k < NZ; ++k)
      {	    
	      rho = exp(-T0_H * ((i - MIDDLE_X) / SD) * ((i - MIDDLE_X) / SD) * DEF_X - T0_H * ((j - MIDDLE_Y) / SD) * ((j - MIDDLE_Y) / SD) * DEF_Y
              - T0_H * LAMBDA_SQR * ((k - MIDDLE_Z) / SD) * ((k - MIDDLE_Z) / SD) * DEF_Z) * T_EFF;

        #pragma simd
	      for (int n = 0; n < Q; ++n)
	      {
		      c_sqr = (c[n][0] * c[n][0]) + (c[n][1] * c[n][1]) + (c[n][2] * c[n][2]);
      		fIn[i][j][k][n] = wi[n] * rho * (1.0 + ((T_EFF - 1.0) / (2.0 * T0)) * (c_sqr - 3.0 * T0));
		      fOut[i][j][k][n] = fIn[i][j][k][n];	  
	      }
	    }
	  }     
  }    
}

void eq(qr_type**** fIn, qr_type**** fOut, qr_type*** rho, qr_type*** ux, qr_type*** uy, qr_type*** uz, const double c[Q][D], const double wi[Q], const bool& ftrue)
{
  qr_type x, y, z, temp, Teff, fEq, force;
  qr_type u_sqr, c_dot_u, c_sqr;
  qr_type Ptot = 0.0;
  qr_type Ntot = 0.0;

  for (int i = 0; i < NX; ++i)
  {
    for (int j = 0; j < NY; ++j)
    {
      for (int k = 0; k < NZ; ++k)
      {           
      	x = (i - MIDDLE_X) / SD;
      	y = (j - MIDDLE_Y) / SD;
      	z = LAMBDA_SQR * (k - MIDDLE_Z) / SD;

      	// Set to zero before summing.
      	rho[i][j][k] = qr_type(0.0);
      	ux[i][j][k] = qr_type(0.0);
      	uy[i][j][k] = qr_type(0.0);
      	ux[i][j][k] = qr_type(0.0);
      	            
        //#pragma simd
      	for (int n = 0; n < Q; ++n)
      	{
          temp = fIn[i][j][k][n];

      	  rho[i][j][k] += temp;
      	  ux[i][j][k] += (c[n][0] * temp);
      	  uy[i][j][k] += (c[n][1] * temp);
      	  uz[i][j][k] += (c[n][2] * temp);

      	  Ptot += temp * ((c[n][0] * c[n][0]) + (c[n][1] * c[n][1]) + (c[n][2] * c[n][2]));
      	}

      	// Divide Speeds by rho.
        temp = rho[i][j][k];
      	ux[i][j][k] /= temp;
      	uy[i][j][k] /= temp;
      	uz[i][j][k] /= temp;

      	if (ftrue)
        {
      	  ux[i][j][k] -= x / SD2;
      	  uy[i][j][k] -= y / SD2;
      	  uz[i][j][k] -= z / SD2;
      	}

      	Ptot -= rho[i][j][k] * ((ux[i][j][k] * ux[i][j][k]) + (uy[i][j][k] * uy[i][j][k]) + (uz[i][j][k] * uz[i][j][k]));
      	Ntot += rho[i][j][k];
      }
    }
  }

  Teff = Ptot / Ntot / (3.0 * T0);

  for (int i = 0; i < NX; ++i)
  {
    for (int j = 0; j < NY; ++j)
    {
      for (int k = 0; k < NZ; ++k)
      {
        x = (i - MIDDLE_X) / SD;
        y = (j - MIDDLE_Y) / SD;
        z = LAMBDA_SQR * (k - MIDDLE_Z) / SD;

      	u_sqr = (ux[i][j][k] * ux[i][j][k]) + (uy[i][j][k] * uy[i][j][k]) + (uz[i][j][k] * uz[i][j][k]);

        //#pragma simd
      	for (int n = 0; n < Q; ++n)
        {
      	  c_dot_u = (c[n][0] * ux[i][j][k]) + (c[n][1] * uy[i][j][k]) + (c[n][2] * uz[i][j][k]);
      	  c_sqr = (c[n][0] * c[n][0]) + (c[n][1] * c[n][1]) + (c[n][2] * c[n][2]);

      	  fEq = wi[n] * rho[i][j][k] * (1.0 + (c_dot_u/T0) * (1.0 + ((Teff - 1.0) / (2.0 * T0)) * (c_sqr - 5.0 * T0)) + (c_dot_u * c_dot_u) / (2.0 * T0 * T0) - u_sqr / (2.0 * T0) + ((Teff - 1.0) / (2.0 * T0)) * (c_sqr - 3.0 * T0)
                + (c_dot_u * c_dot_u * c_dot_u) / (6.0 * T0 * T0 * T0) - (u_sqr * c_dot_u) / (2.0 * T0 * T0));

      	  if (ftrue)
          {
      	    force = -(DT / SD) * CA * CA * (1.0 - 0.5 * OMEGA / (CA * CA * DT)) * wi[n] * rho[i][j][k] * (((c[n][0] * x) + (c[n][1] * y) + (c[n][2] * z)) / T0 - ((x * ux[i][j][k]) + (y * uy[i][j][k])
                    + (z * uz[i][j][k])) / T0) * (1.0 + ((Teff - 1.0) / (2.0 * T0)) * (c_sqr - 4.0 * T0) + c_dot_u / T0 + (c_dot_u * c_dot_u) / (2.0 * T0 * T0) - u_sqr / (2.0 * T0));
      	  }
      	  else
          {
      	    force = 0.0;
      	  }

      	  fOut[i][j][k][n] = ((1.0 - OMEGA) * fIn[i][j][k][n]) + (OMEGA * fEq) + force;
      	}
      }
    }
  }
}

void stream(qr_type**** fIn, qr_type**** fOut, const double c[Q][D])
{
  int in, jn, kn;

  for (int i = 0; i < NX; ++i)
  {
    for (int j = 0; j < NY; ++j)
    {
      for (int k = 0; k < NZ; ++k)
      {
        #pragma simd
	      for(int n = 0; n < Q; ++n)
        {	  
      	  in = i + int(c[n][0]);
      	  jn = j + int(c[n][1]);
      	  kn = k + int(c[n][2]);

      	  if (in > NX - 1 || in < 0)
          {
            //in %= NX;
            in = (in + NX) % NX;
      	  }
      	  if (jn > NY - 1 || jn < 0)
          {
      	    //jn %= NY;
            jn = (jn + NY) % NY;
      	  } 
      	  if (kn > NZ - 1 || kn < 0)
          {
      	    //kn %= NZ;
            kn = (kn + NZ) % NZ;
      	  }

      	  fIn[in][jn][kn][n] = fOut[i][j][k][n];
	      }
      }
    }                    
  }
}

void write_gaussian(qr_type*** rho, qr_type*** ux, const int& ts)
{
  fstream out;
  char fname[BUFF_SIZE];

  snprintf(fname, BUFF_SIZE, "%s/Xrho_t%i.dat", WRITE_DIR, ts);
  out.open(fname, ios::out);
  for(int i = 0; i < NX; ++i)
  {
    out << (i - MIDDLE_X) * SIN_V << "\t";
    out << rho[i][MIDDLE_Y][MIDDLE_Z] << "\n";
  }        
  out.close(); 

  snprintf(fname, BUFF_SIZE, "%s/Yrho_t%i.dat", WRITE_DIR, ts);
  out.open(fname, ios::out);
  for(int j = 0; j < NY; ++j)
  {
    out << (j - MIDDLE_Y) * SIN_V << "\t";
    out << rho[MIDDLE_X][j][MIDDLE_Z] << "\n";
  }        
  out.close(); 

  snprintf(fname, BUFF_SIZE, "%s/Zrho_t%i.dat", WRITE_DIR, ts);
  out.open(fname, ios::out);
  for(int k = 0; k < NZ; ++k)
  {
    out << (k - MIDDLE_Z) * SIN_V << "\t";
    out << rho[MIDDLE_X][MIDDLE_Y][k] << "\n";
  }        
  out.close(); 

  snprintf(fname, BUFF_SIZE, "%s/Xux_t%i.dat", WRITE_DIR, ts);
  out.open(fname, ios::out);
  for(int i = 0; i < NX; ++i)
  {
    out << (i - MIDDLE_X) * SIN_V << "\t";
    out << setprecision(15) << ux[i][MIDDLE_Y][MIDDLE_Z] << "\n";
  }        
  out.close();  
}

void get_walltime(double& wcTime)
{
  struct timeval tp;
  gettimeofday(&tp, NULL);
  wcTime = (double)(tp.tv_sec + tp.tv_usec/1000000.0);
}

#include "functions.h"
#include <iostream>
using namespace std;

void init_gaussian(double* fIn, vector<vector<vector<double>>>* fOut, vector<vector<double>>* rho,
                   vector<vector<double>>* ux, vector<vector<double>>* uy, int c[Q][D], double wi[Q], double lambda, int nx,
                   int ny, double sd, double T0, double omega)
{

  double u_sqr, c_dot_u, fEq;
  double x, y;

  double middlex = nx / 2;
  double middley = ny / 2;

  for (int i = 0; i < nx; ++i)
  {
    for (int j = 0; j < ny; ++j)
    {

      u_sqr = ((*ux)[i][j] * (*ux)[i][j]) + ((*uy)[i][j] * (*uy)[i][j]);

      // Data type mixing possible precision loss.
      (*rho)[i][j] = exp(-(0.5 / T0) * ((i - middlex) / sd) * ((i - middlex) / sd) - (0.5 / T0) * lambda
                     * lambda * ((j - middley) / sd) * ((j - middley) / sd));

      for (int n = 0; n < Q; ++n)
	    {
        c_dot_u = (c[n][0] * (*ux)[i][j]) + (c[n][1] * (*uy)[i][j]);

	      fEq = wi[n] * (*rho)[i][j];

	      (*fIn)[i][j][n] = fEq;
	      (*fOut)[i][j][n] = (*fIn)[i][j][n];
	    }
    }
  }
}

void eq_and_stream(vector<vector<vector<double>>>* fIn, vector<vector<vector<double>>>* fOut, vector<vector<double>>* rho,
                   vector<vector<double>>* ux, vector<vector<double>>* uy, int c[Q][D], double wi[Q], int nop[Q], double lambda,
                   int nx, int ny, double T0, double omega, double sd, bool ftrue)
{

  double u_sqr, c_dot_u, force;
  int in, jn;
  double x, y;
  double fEq[Q];

  double check_rho, check_ux, check_uy;

  double middlex = nx / 2;
  double middley = ny / 2;

  double c_sqr;

  double cdotX,udotX;

  for (int i = 0; i < nx; ++i)
  {
    for (int j = 0; j < ny; ++j)
    {
      // More possible precision loss.
      x = (i - middlex) / sd;
      y = lambda * lambda * (j - middley) / sd;

      //set to zero before summing
      (*rho)[i][j] = 0.0;
      (*ux)[i][j] = 0.0;
      (*uy)[i][j] = 0.0;

      for (int n = 0; n < Q; ++n)
	    {
        (*rho)[i][j] = (*rho)[i][j] + (*fIn)[i][j][n]; // rho
	      (*ux)[i][j] = (*ux)[i][j] + (c[n][0] * (*fIn)[i][j][n]); // ux
	      (*uy)[i][j] = (*uy)[i][j] + (c[n][1] * (*fIn)[i][j][n]); // uy
	    }

      (*ux)[i][j] = (*ux)[i][j] / (*rho)[i][j];
      (*uy)[i][j] = (*uy)[i][j] / (*rho)[i][j];

      u_sqr = ((*ux)[i][j] * (*ux)[i][j]) + ((*uy)[i][j] * (*uy)[i][j]);

      for (int n = 0; n < Q; ++n)
      {

      	c_dot_u = (c[n][0] * (*ux)[i][j]) + (c[n][1] * (*uy)[i][j]);

      	c_sqr = (c[n][0] * c[n][0]) + (c[n][1] * c[n][1]);

      	fEq[n] = wi[n] * (*rho)[i][j] * (1.0 + (1.0 / T0) * c_dot_u + (1.0 / (2.0 * T0 * T0)) * c_dot_u * c_dot_u - (1.0 / (2.0 * T0)) * u_sqr);

      	cdotX = (c[n][0] * x) + (c[n][1] * y);
      	udotX = ((*ux)[i][j] * x) + ((*uy)[i][j] * y);

      	if (ftrue)
        {
      	  force = -(1.0 / T0) * (cdotX - udotX) * fEq[n] / sd;
      	  if (i == 0 && j == 0 && n == 0)
          {
      	    cout << "potential is on" << endl;
      	  }
      	}
      	else
        {
      	  force = 0.0;
      	  if (i == 0 && j == 0 && n == 0)
          {
      	    cout << "potential is off" << endl;
      	  }
      	}

	      (*fIn)[i][j][n] = ((*fIn)[i][j][n] * (1.0 - omega)) + omega * fEq[n] + force;
      }
    }
  }
}

void write_gaussian(vector<vector<double>>* rho, vector<vector<double>>* ux, vector<vector<double>>* uy, int nx, int ny, double sd, int ts)
{
  fstream out;
  char fname[255];
  float sinv=1.0 / sd;
  int middlex = nx / 2;
  int middley = ny / 2;

  sprintf(fname,"data/Xrho_t%i.dat", ts);
  //out.open(fname, ios::out);
  out.open("data.txt");
  for(int i = 0; i < nx; ++i)
  {
    int j = middley;
    out << (i - middlex) * sinv << "\t";
    out << (*rho)[i][j] << "\n";
  }

  out.close();

  sprintf(fname,"data/Yrho_t%i.dat", ts);
  out.open(fname, ios::out);
  for(int j = 0; j < ny; ++j)
  {
    int i = middlex;
    out << (j - middley) * sinv << "\t";
    out << (*rho)[i][j] << "\n";
  }

  out.close();

  sprintf(fname,"data/Xux_t%i.dat", ts);
  out.open(fname, ios::out);
  for(int i = 0; i < nx; ++i)
  {
    int j = middley;
    out << (i - middlex) * sinv << "\t";
    out << (*ux)[i][j] << "\n";
  }

  out.close();

  sprintf(fname,"data/Yuy_t%i.dat", ts);
  out.open(fname, ios::out);
  for(int j = 0; j < ny; ++j)
  {
    int i = middlex;
    out << (j - middley) * sinv << "\t";
    out << (*uy)[i][j] << "\n";
  }

  out.close();
}

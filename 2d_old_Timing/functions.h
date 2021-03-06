#ifndef __FUNCTIONS_H_
#define __FUNCTIONS_H_
#include <vector>
#include <math.h>
#include <fstream>
#include <sys/time.h>

using namespace std;

#define Q 9
#define D 2
//Function description,input parameter,output parameters??
void init_gaussian(vector<vector<vector<double> > >* fIn, vector<vector<vector<double> > >* fOut, vector<vector<double> >* rho,
		   vector<vector<double> >* ux, vector<vector<double> >* uy, int c[Q][D], double wi[Q], double lambda, int nx,
		   int ny, double sd, double T0, double omega);


//Function description,input parameter,output parameters??
void eq_and_stream(vector<vector<vector<double> > >* fIn, vector<vector<vector<double> > >* fOut, vector<vector<double> >* rho,
		   vector<vector<double> >* ux, vector<vector<double> >* uy, int c[Q][D], double wi[Q], int nop[Q], double lambda,
		   int nx, int ny, double T0, double omega, double sd, bool ftrue);


//Function description,input parameter,output parameters??
void write_gaussian(vector<vector<double> >* rho, vector<vector<double> >* ux, vector<vector<double> >* uy, int nx, int ny, double sd, int ts);

double get_walltime();

#endif

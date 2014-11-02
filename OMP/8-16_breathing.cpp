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
#include <sys/time.h>
#include <omp.h>

#define Q 9
#define D 2

using namespace std;

void init_gaussian(vector<vector<vector <double> > > * fIn,vector<vector<vector <double> > > * fOut,vector<vector <double> > * rho, vector<vector <double> > * ux, vector<vector <double> > * uy, int c[Q][D], double wi[Q],double lambda, int nx, int ny, double sd, double T0, double omega) {
  
  double u_sqr, c_dot_u, fEq;
  double x, y;

  double middlex = nx/2;
  double middley = ny/2;
  int i,j,n;
#pragma omp parallel for schedule(static) collapse(2) private(i,j,n) shared(u_sqr, rho, c_dot_u,fEq,x,y,middlex,middley)

  for (i = 0; i < nx ; i++) {
    for (j = 0; j < ny ; j++) {

      u_sqr = (*ux)[i][j]*(*ux)[i][j] + (*uy)[i][j]*(*uy)[i][j]; 

      (*rho)[i][j] = exp( -(0.5/T0)*((i-middlex)/sd)*((i-middlex)/sd) - (0.5/T0)*lambda*lambda*((j-middley)/sd)*((j-middley)/sd) );
      #pragma simd
      for (n = 0; n < Q; n++)
	{
	  c_dot_u = c[n][0]*(*ux)[i][j] + c[n][1]*(*uy)[i][j]; 
	 
	  fEq = wi[n]*((*rho)[i][j]);
	  //#pragma omp critical
	  (*fIn)[i][j][n] = fEq;
	  (*fOut)[i][j][n] = (*fIn)[i][j][n];
	}
    }
  }
}

void eq_and_stream(vector<vector<vector <double> > > * fIn,vector<vector<vector <double> > > * fOut,vector<vector <double> > * rho, vector<vector <double> > * ux, vector<vector <double> > * uy, int c[Q][D], double wi[Q], int nop[Q], double lambda, int nx, int ny, double T0, double omega, double sd, int ftrue) {

  double u_sqr, c_dot_u, force;
  int in, jn;
  double x, y;
  double fEq[Q];

  double check_rho, check_ux, check_uy;

  double middlex = nx/2;
  double middley = ny/2;

  double c_sqr;

  double cdotX,udotX;
  int i,j,n;
  #pragma omp parallel for schedule(static) private(i,j,n,c_dot_u,c_sqr) collapse(2) shared(rho,ux,uy,u_sqr,force,cdotX,udotX,x,y,middlex,middley) 
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {

      x = (i-middlex)/sd;
      y = lambda*lambda*(j-middley)/sd;

      //set to zero before summing
      (*rho)[i][j]=0.0;
      (*ux)[i][j]=0.0;
      (*uy)[i][j]=0.0;
  	            
      for (n = 0; n < Q; n++)
	{
	  (*rho)[i][j] = (*rho)[i][j]+(*fIn)[i][j][n]; // rho
	  (*ux)[i][j] = (*ux)[i][j]+(c[n][0])*((*fIn)[i][j][n]); // ux
	  (*uy)[i][j] = (*uy)[i][j]+(c[n][1])*((*fIn)[i][j][n]); // uy
	}

      (*ux)[i][j] = (*ux)[i][j]/(*rho)[i][j];
      (*uy)[i][j] = (*uy)[i][j]/(*rho)[i][j];

      u_sqr = ((*ux)[i][j])*((*ux)[i][j]) + ((*uy)[i][j])*((*uy)[i][j]);
     //#pragma omp parallel for private(c_dot_u,c_sqr)
      for (n = 0; n < Q; n++) {

	c_dot_u = (c[n][0])*((*ux)[i][j]) + (c[n][1])*((*uy)[i][j]);

	c_sqr = (c[n][0])*(c[n][0]) + (c[n][1])*(c[n][1]);

	fEq[n] = wi[n]*((*rho)[i][j])*(1.0 + (1.0/T0)*c_dot_u + (1.0/(2.0*T0*T0))*c_dot_u*c_dot_u - (1.0/(2.0*T0))*u_sqr);

	cdotX = (c[n][0])*x + (c[n][1])*y;
	udotX = x*((*ux)[i][j]) + y*((*uy)[i][j]);

	if (ftrue == 1) {
	  force = -(1.0/T0)*(cdotX-udotX)*fEq[n]/sd;
	  // if (i==0 && j==0 && n==0) {
	  //   cout << "potential is on" << endl;
	  // }
	}
	else if (ftrue == 0) {
	  force = 0.0;
	  // if (i==0 && j==0 && n==0) {
	  //   cout << "potential is off" << endl;
	  // }
	}
	else {
	  cout << "ftrue should be either zero or one, please and thanks." << endl;
	}

	(*fIn)[i][j][n] = ((*fIn)[i][j][n])*(1.0-omega) + omega*fEq[n] + force;
      }
    }
  }

      // Check The moments:
  //     check_rho = 0.0;
  //     check_ux = 0.0;
  //     check_uy = 0.0;
  //     for (int n = 0; n < Q; n++) {

  //     	check_rho = check_rho + fEq[n];
  //     	check_ux = check_ux + (c[n][0])*fEq[n];
  //     	check_uy = check_uy + (c[n][1])*fEq[n];
  //     }

  //     check_ux = check_ux/check_rho;
  //     check_uy = check_uy/check_rho;

  //     if (i==nx/2 && j==ny/5) {
  //     	cout << "after_rho = " << check_rho << "\t after_ux = " << check_ux << "\t after_uy = " << check_uy << endl;
  //     }
  //   }
  // }
  
  #pragma omp parallel for collapse(2) shared(in,jn)
  for (i=0; i < nx; i++) {
    for (j=0; j < ny; j++) {
  //#pragma simd
      for (n=0; n < Q; n++) {

	in = i + int(c[n][0]);
	jn = j + int(c[n][1]);
	
	if (in > nx-1 || in < 0) {
	  in = (in+nx)%nx;
	  // in = i; 
	}
	if (jn > ny-1 || jn < 0) {
	  jn = (jn+ny)%ny;
	  // jn = j;
	}
 //#pragma omp critical
	(*fOut)[i][j][nop[n]] = (*fIn)[in][jn][nop[n]];
      }
    }
  }
}

void get_walltime(double* wcTime) {
  struct timeval tp;
  gettimeofday(&tp, NULL);
  *wcTime = (double)(tp.tv_sec + tp.tv_usec/1000000.0);
}


void write_gaussian(vector<vector <double> > * rho, vector<vector <double> > * ux, vector<vector <double> > * uy, int nx, int ny, double sd, int ts) {

  fstream out;
  char fname[255];
  float sinv=1.0/sd;
  int middlex = nx/2;
  int middley = ny/2;

  sprintf(fname,"data251_omp/Xrho_t%i.dat",ts);
  out.open(fname, ios::out);
  for(int i=0; i < nx; i++){
    int j=middley;
    out << (i-middlex)*sinv << "\t";
    out << (*rho)[i][j] << "\n";
  }  
        
  out.close(); 

  sprintf(fname,"data251_omp/Yrho_t%i.dat",ts);
  out.open(fname, ios::out);
  for(int j=0; j < ny; j++){
    int i=middlex;
    out << (j-middley)*sinv << "\t";
    out << (*rho)[i][j] << "\n";
  }  
        
  out.close(); 

  sprintf(fname,"data251_omp/Xux_t%i.dat",ts);
  out.open(fname, ios::out);
  for(int i=0; i < nx; i++){
    int j=middley;
    out << (i-middlex)*sinv << "\t";
    out << (*ux)[i][j] << "\n";
  }  
        
  out.close(); 

  sprintf(fname,"data251_omp/Yuy_t%i.dat",ts);
  out.open(fname, ios::out);
  for(int j=0; j < ny; j++){
    int i=middlex;
    out << (j-middley)*sinv << "\t";
    out << (*uy)[i][j] << "\n";
  }  
        
  out.close(); 

}

int main(int argc, const char * argv[])
{

  int nx = atoi( argv[1] );
  int ny = atoi( argv[2] );
  double sd = atof( argv[3] );
  int nsteps = atoi( argv[4] ); 

  //int nx = 251;
  //int ny = 251;
  //int nsteps = 250;
  //double sd = 50.0;
  double omega = 1.0;
  int meas_steps = 10;

  double dtreal = 12; // micro-seconds
  double wxreal = 0.000785398; // micro-Hz
  int dtsim = round( dtreal*wxreal*sd );
  // double lambda = 0.0345955;

  // double lambda = double(nx)/double(ny);
  // double lambda = 0.035;
  double lambda = 1.0;
 
  // cout << "dtsim = " << dtsim << endl;

  int t_off = 10; // time step to turn potential off
  int t_on = t_off + dtsim; // time step to turn it back on

  int ftrue = 1; // 1 if potential on, 0 if potential off

  double u_sqr, c_dot_u, force;
  int in, jn;

  double T0 = 1.0/3.0;

  vector<vector<vector <double> > > fIn (nx, vector<vector <double> >(ny, vector <double>(Q,0.0)));

  vector<vector<vector <double> > > fOut (nx, vector<vector <double> >(ny, vector <double>(Q,0.0)));

  vector<vector <double> > rho (nx, vector <double> (ny,1.0));
  vector<vector <double> > ux (nx, vector <double> (ny,0.0));
  vector<vector <double> > uy (nx, vector <double> (ny,0.0));

  int c[Q][D] = {{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};

  int nop[Q] = {0,3,4,1,2,7,8,5,6};

  double w1 = 4.0/9.0;
  double w2 = 1.0/9.0;
  double w3 = 1.0/36.0;
  double wi[] = {w1,w2,w2,w2,w2,w3,w3,w3,w3}; 

  init_gaussian(&fIn,&fOut,&rho,&ux,&uy,c,wi,lambda,nx,ny,sd,T0,omega);   

  double S,E;

  get_walltime(&S);

  for (int ts=0; ts < nsteps; ts++) {

    if (ts == t_off) {
      ftrue = 0;
    }
    if (ts == t_on) {
      ftrue = 1;
    }

    eq_and_stream(&fIn,&fOut,&rho,&ux,&uy,c,wi,nop,lambda,nx,ny,T0,omega,sd,ftrue);

    fIn = fOut;
    
    if (ts%meas_steps==0) {
      write_gaussian(&rho,&ux,&uy,nx,ny,sd,ts);
    }
  }

  get_walltime(&E);

  printf("Problem Size: nx=%d, ny=%d, sd=%f, nsteps=%d",nx,ny,sd,nsteps);
  printf("Walltime %f s \n",E-S);

  return 0;
}

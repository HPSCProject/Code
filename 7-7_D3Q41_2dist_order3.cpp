
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
// #include <omp.h>

#define Q 41
#define D 3

using namespace std;

//calculate velocity fields and density
void calc_marco(int nx, int ny, int nz, vector<vector<vector<vector <double> > > > * f, vector<vector<vector<vector <double> > > > * gIn, double c[Q][D],int* ts, double T0){

  // double u_sqr, chi_sqr;
 
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {  
           
	//set to zero before summing
	(*f)[i][j][k][Q]=0.0;
	(*f)[i][j][k][Q+1]=0.0;
	(*f)[i][j][k][Q+2]=0.0;
	(*f)[i][j][k][Q+3]=0.0;
	(*f)[i][j][k][Q+4]=0.0;
	            
	for (int n = 0; n < Q; n++)
	  {

	    (*f)[i][j][k][Q]=(*f)[i][j][k][Q]+(*f)[i][j][k][n]; // rho
	    (*f)[i][j][k][Q+1]=(*f)[i][j][k][Q+1]+((c[n][0])*(*f)[i][j][k][n]); // ux
	    (*f)[i][j][k][Q+2]=(*f)[i][j][k][Q+2]+((c[n][1])*(*f)[i][j][k][n]); // uy
	    (*f)[i][j][k][Q+3]=(*f)[i][j][k][Q+3]+((c[n][2])*(*f)[i][j][k][n]); // uz

	  }

	//Divide Speeds by rho
	(*f)[i][j][k][Q+1]=((*f)[i][j][k][Q+1])/(*f)[i][j][k][Q]; // ux/rho
	(*f)[i][j][k][Q+2]=((*f)[i][j][k][Q+2])/(*f)[i][j][k][Q]; // uy/rho
	(*f)[i][j][k][Q+3]=((*f)[i][j][k][Q+3])/(*f)[i][j][k][Q]; // uz/rho

	// calculate the energy
	for (int n=0; n<Q; n++) {

	  (*f)[i][j][k][Q+4] = (*f)[i][j][k][Q+4] + (*gIn)[i][j][k][n];
	}


	// if (chi_sqr < -pow(10,-15)) {
	//   cout << "(" << i << ", " << j << ", " << k << "): chi_sqr: " << chi_sqr << "\t ux: " << (*f)[i][j][k][Q+1] << "\t uy: " << (*f)[i][j][k][Q+2] << "\t uz: " << (*f)[i][j][k][Q+3] << "\t rho: " << (*f)[i][j][k][Q] << "\t theta_before: " << theta_before << endl;
	//   // for (int n =0; n < Q; n++) {
	//   //   cout << n << "\t " << (*f)[i][j][k][n] << endl;
	//   // }
	// }
  
	// (*f)[i][j][k][Q+5] = (*f)[i][j][k][Q]/pow(chi_sqr,3.0/2.0); // n

      }
    }
  }
}

void eq(int nx, int ny, int nz, vector<vector<vector<vector <double> > > > *fIn,vector<vector<vector<vector <double> > > > *fOut,vector<vector<vector<vector <double> > > > *gIn,vector<vector<vector<vector <double> > > > *gOut,double c[Q][D], double wi[Q], double sd, double omega_m, double omega_E, double lambda,int* ts, double T0){

   // double begin_ux, begin_uy, begin_uz, begin_rho, begin_E;
   // double after_ux, after_uy, after_uz, after_rho, after_E;

   // double precision = pow(10,-10);

   double c_dot_u, u_sqr, chi_sqr, c_dot_c, c_sqr;

   // double max_diff = 0.0;

   // double after_usqr;

   double T;
   // double T0 = 1.0 - sqrt(2.0/5.0);

   double fEq[Q];
   double gEq[Q];
   double force[Q];

   double cdotX, udotX;
   double ee;

   // double after_chi_sqr;

   double x, y, z;
   double middlex = nx/2;
   double middley = ny/2;
   double middlez = nz/2;

   // give the central value of temperature:
   // cout << 2.0*((*fIn)[nx/2][ny/2][nz/2][23]) << endl;

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {

	x = (i-middlex)/sd;
	y = (j-middley)/sd;
	z = lambda*lambda*(k-middlez)/sd;

	// ensure that E = 3T0 + usqr:
	u_sqr = ((*fIn)[i][j][k][Q+1])*((*fIn)[i][j][k][Q+1]) + ((*fIn)[i][j][k][Q+2])*((*fIn)[i][j][k][Q+2]) + ((*fIn)[i][j][k][Q+3])*((*fIn)[i][j][k][Q+3]);

	// T = (*fIn)[i][j][k][Q+4]/(3.0*(*fIn)[i][j][k][Q]);
	ee = (*fIn)[i][j][k][Q+4];

	T = 2.0*((ee-0.5*((*fIn)[i][j][k][Q])*u_sqr)/(3.0*(*fIn)[i][j][k][Q]));

	chi_sqr = T/T0;

	// begin_rho = (*fIn)[i][j][k][Q];
	// begin_ux = (*fIn)[i][j][k][Q+1];
	// begin_uy = (*fIn)[i][j][k][Q+2];
	// begin_uz = (*fIn)[i][j][k][Q+3];
	// begin_E = (*fIn)[i][j][k][Q+4];

	// if (i==nx/2 && j==ny/5 && k==nz/2) {
	//   cout << setprecision(15) << "begin_ux: " << begin_ux << "\t begin_uy: " << begin_uy << "\t begin_uz: " << begin_uz << "\t begin u_sqr: " << u_sqr <<"\t begin_rho: " << begin_rho << "\t begin chi_sqr: " << chi_sqr << "\t begin_E: " << begin_E << endl;
	// }

	for (int n = 0; n < Q; n++) {

	  c_dot_u = c[n][0]*(*fIn)[i][j][k][Q+1] + c[n][1]*(*fIn)[i][j][k][Q+2] + c[n][2]*(*fIn)[i][j][k][Q+3]; 

	  c_sqr = (c[n][0]*c[n][0])+(c[n][1]*c[n][1])+(c[n][2]*c[n][2]);


	  // fEq[n] = wi[n]*((*fIn)[i][j][k][Q+5])*pow(chi_sqr,3.0/2.0)*(1.0 + ((5.0-sqrt(10.0)-5.0*chi_sqr*T0)/(-5.0+sqrt(10.0))) + (0.0*(chi_sqr-1.0)*((c_sqr/T0)-3.0)) + (1.0 + 0.5*(chi_sqr-1.0)*((c_sqr/T0)-5.0))*(c_dot_u/T0) + (1.0/(2.0*T0*T0))*c_dot_u*c_dot_u - (u_sqr/(2.0*T0)) + (1.0/(6.0*T0*T0*T0))*(c_dot_u*c_dot_u*c_dot_u) - (1.0/(2.0*T0*T0))*c_dot_u*u_sqr);

	  // fEq[n] = wi[n]*((*fIn)[i][j][k][Q])*(chi_sqr + c_dot_u/T0 + (1.0/(2.0*T0*T0)*c_dot_u*c_dot_u - (1.0/(2.0*T0))*u_sqr)) - (wi[n]/3.0)*((c_sqr - 3.0*T0)/(2.0*T0*T0))*(*fIn)[i][j][k][Q+5];

	  fEq[n] = wi[n]*((*fIn)[i][j][k][Q])*(1.0 - ((-5.0*T0 + 5.0*chi_sqr*T0)/(-5.0*T0)) + (1.0 + 0.5*(chi_sqr-1.0)*((c_sqr/T0)-5.0))*(c_dot_u/T0) + (1.0/(2.0*T0*T0))*c_dot_u*c_dot_u - u_sqr/(2.0*T0) + (1.0/(6.0*T0*T0*T0))*c_dot_u*c_dot_u*c_dot_u - (1.0/(2.0*T0*T0))*c_dot_u*u_sqr) - (1.0/3.0)*wi[n]*((c_sqr-3.0*T0)/(2.0*T0*T0))*(*fIn)[i][j][k][Q+5];
	  
	  gEq[n] = wi[n]*((*fIn)[i][j][k][Q])*(1.5*chi_sqr*T0 + 0.5*u_sqr + ((5.0/2.0)*chi_sqr*T0 + 0.5*u_sqr)*(c_dot_u/T0 + (1.0/(2.0*T0*T0))*c_dot_u*c_dot_u - (1.0/(2.0*T0))*u_sqr));

	  cdotX = (c[n][0])*x + (c[n][1])*y + (c[n][2])*z;
	  udotX = x*((*fIn)[i][j][k][Q+1]) + y*((*fIn)[i][j][k][Q+2]) + z*((*fIn)[i][j][k][Q+3]);

	  force[n] = -(1.0/T0)*(cdotX-udotX)*fEq[n]/sd;

	  // gEq[n] = wi[n]*((*fIn)[i][j][k][Q])*((2332800.0*(6.0*chi_sqr*T0 + u_sqr))/((-165832612.0+52295825.0*sqrt(10))*u_sqr))*((c_dot_u/T0) + (1.0/(2.0*T0*T0))*c_dot_u*c_dot_u - u_sqr/(2.0*T0) + (1.0/(6.0*T0*T0*T0))*c_dot_u*c_dot_u*c_dot_u - (1.0/(2.0*T0*T0))*c_dot_u*u_sqr);

	  // if (fEq[n] < 0.0) {

	  //   cout << n << "\t fEq: " << fEq[n] << "\t chi_sqr: " << chi_sqr << "\t ux: " << (*fIn)[i][j][k][Q+1] << "\t uy: " << (*fIn)[i][j][k][Q+2] << "\t uz: " << (*fIn)[i][j][k][Q+3] << endl;

	  // }
	}

	fEq[0] = wi[0]*((*fIn)[i][j][k][Q])*(1.0 - ((-8065.0+3014.0*sqrt(10.0))*(-5.0*T0+5.0*chi_sqr*T0))/(2.0*(-5.0*T0)*(-5045.0+1507.0*sqrt(10.0))) - u_sqr/(2.0*T0)) - (1.0/3.0)*wi[0]*(-1.5/T0)*(*fIn)[i][j][k][Q+5];

	cdotX = (c[0][0])*x + (c[0][1])*y + (c[0][2])*z;

	force[0] = -(1.0/T0)*(cdotX-udotX)*fEq[0]/sd;

	// fEq[0] = wi[0]*((*fIn)[i][j][k][Q])*(3.0 - 2.0*chi_sqr - 1.5*u_sqr);
	// fEq[0] = wi[0]*((*fIn)[i][j][k][Q])*((405.0*(-5.0+sqrt(10.0)) + (-8065.0+3014.0*sqrt(10.0))*chi_sqr*T0)/(2.0*(-8059.0+2516.0*sqrt(10.0))) - u_sqr/(2.0*T0)) - (wi[0]/3.0)*((-1.5/T0)*((*fIn)[i][j][k][Q+5]));

	// if (i==nx/2 && j==ny/5 && k==nz/2) {

	//   cout << "n: " << (*fIn)[i][j][k][Q+5] << "\t chi_sqr: " << chi_sqr << "\t ux: " << (*fIn)[i][j][k][Q+1] << "\t uy: " << (*fIn)[i][j][k][Q+2] << "\t uz: " << (*fIn)[i][j][k][Q+3] << "\t F: " << endl;

	//   for (int n = 0; n < Q; n++) {
	    
	//     // cout << n << "\t" << fEq[n] << endl;
	//     cout << n << "\t" << wi[n]*((*fIn)[i][j][k][Q+5])*pow(chi_sqr,3.0/2.0) << endl;

	//     cout << n << "\t" << wi[n] << endl;
	//   }
	// }

	// calculate the trace of the stress tensor Pi
	(*fIn)[i][j][k][Q+5] = 0.0;
	for (int n = 0; n < Q; n++) {

	  (*fIn)[i][j][k][Q+5] = (*fIn)[i][j][k][Q+5] + (c[n][0]*c[n][0] + c[n][1]*c[n][1] + c[n][2]*c[n][2])*((*fIn)[i][j][k][n] - fEq[n]);

	  (*fOut)[i][j][k][n] = (1.0-omega_m)*(*fIn)[i][j][k][n] + omega_m*(fEq[n]) + force[n];
	  (*gOut)[i][j][k][n] = (1.0-omega_E)*(*gIn)[i][j][k][n] + omega_E*(gEq[n]);


	  // if (i==nx/2 && j==ny/2 && k==nz/2) {
	  //   cout << "diff: " << (*fIn)[i][j][k][n] - fEq[n] << endl;
	  // }
	}

	// if (i==nx/2 && j==ny/2 && k==nz/2) {

	//   cout << "Trace of Pi: " << (*fIn)[i][j][k][Q+5] << endl;
	// }

	// calculate the moments with the equilibrium distribution to verify they are equal
	// after_rho = 0.0;
	// after_ux = 0.0;
	// after_uy = 0.0;
	// after_uz = 0.0;
	// after_E = 0.0;

	// for (int n = 0; n < Q; n++) {
	//   after_rho = after_rho + (fEq[n]);
	//   after_ux = after_ux + (c[n][0])*(fEq[n]);
	//   after_uy = after_uy + (c[n][1])*(fEq[n]);
	//   after_uz = after_uz + (c[n][2])*(fEq[n]);
	// }

	// after_ux = after_ux/after_rho;
	// after_uy = after_uy/after_rho;
	// after_uz = after_uz/after_rho;

	// for (int n = 0; n < Q; n++) {
	//   after_E = after_E + gEq[n];
	// }

	// after_usqr = after_ux*after_ux + after_uy*after_uy + after_uz*after_uz;

	// after_chi_sqr = 2.0*((after_E-0.5*after_rho*after_usqr)/(3.0*after_rho*T0));

	// if (i==nx/2 && j==ny/5 && k==nz/2) {
	//   cout << "after_ux: " << after_ux << "\t after_uy: " << after_uy <<  "\t after_uz: " << after_uz << "\t after_usqr: " << after_usqr << "\t after_rho: " << after_rho << "\t chi_sqr: " << after_chi_sqr << "\t after_E: " << after_E << endl;
	// }

	//omega = (1.0/((1.0/(*omega_val))+(0.5*chi_sqr/sd)))*(chi_sqr/sd);
	// omega_m = chi_sqr/(0.5+0.5*chi_sqr);
	// omega = 0.5;

	// omega_m = 1.0;
	// omega_E = 1.5;
       

	// if (i==nx/2 && j==ny/2 && k==nz/2) {
	//   cout << "omega: " << omega << endl;
	// }

	// for (int n = 0; n < Q; n++) {

	//   (*fOut)[i][j][k][n] = (1.0-omega_m)*(*fIn)[i][j][k][n] + omega_m*(fEq[n]) + force[n];
	//   (*gOut)[i][j][k][n] = (1.0-omega_E)*(*gIn)[i][j][k][n] + omega_E*(gEq[n]);

	  // if (i==nx/2 && j==ny/2 && k==nz/2) {

	  //   cout << n << "\t fEq: " << fEq[n] << "\t chi_sqr: " << chi_sqr << "\t u_sqr: " << u_sqr << endl;

	  //   // if (n>=27) { cout << "3" << endl; }
	  // }


	  // if (isnan((*fOut)[i][j][k][n])) { cout << "n: " << n << "\t fIn: " << (*fIn)[i][j][k][n] << "\t fOut: " << (*fOut)[i][j][k][n] << "\t fEq: " << fEq[n] << endl;
	  // }
	// }
      }
    }
  }

// cout << "max_diff: " << max_diff << endl;
}

void stream(int nx, int ny, int nz, vector<vector<vector<vector <double> > > > *fIn,vector<vector<vector<vector <double> > > > *fOut, vector<vector<vector<vector <double> > > > *gIn,vector<vector<vector<vector <double> > > > *gOut, double c[Q][D], double sd){

  int in,jn,kn;

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) { 
	for(int n = 0; n< Q; n++){
	  
	  in = i + int(c[n][0]);
	  jn = j + int(c[n][1]);
	  kn = k + int(c[n][2]);

	  if (in > nx-1 || in < 0) {
	    in = (in+nx)%nx;
	    // in = i;
	  }
	  if (jn > ny-1 || jn < 0) {
	    jn =(jn+ny)%ny;
	    // jn = j;
	  } 
	  if (kn > nz-1 || kn < 0) {
	    kn = (kn+nz)%nz;
	    // kn = k;
	  }

	  (*fIn)[in][jn][kn][n] = (*fOut)[i][j][k][n];

	  (*gIn)[in][jn][kn][n] = (*gOut)[i][j][k][n];
	
	  // (*fIn)[i][j][k][n] = (*fOut)[in][jn][kn][n];
	  // (*gIn)[i][j][k][n] = (*gOut)[in][jn][kn][n];  
	}
      }
    }                    
  }
}//end eq and stream


//make initial lattice
void make_lattice(vector<vector<vector<vector <double> > > > *fIn,vector<vector<vector<vector <double> > > > *fOut,vector<vector<vector<vector <double> > > > *gIn,vector<vector<vector<vector <double> > > > *gOut, double c[Q][D], double wi[Q], int nx, int ny,int nz,double sd,double lambda,double chi_sqr_init, double T0){
    
  //argument of the expontial. 
  double arg;
  double chi_sqr = chi_sqr_init;
  double c_sqr;
  double ux = 0.0;
  double uy = 0.0;
  double uz = 0.0;
  double c_dot_u, u_sqr;
  double ee, rho;

  // #pragma omp parallel for
  for (int i = 0; i < nx ; i++)
    {
      for (int j = 0; j < ny ; j++)
	{
	  for (int k = 0; k < nz; k++) {
	    
	    // lambda = 1:
	    // arg = -(0.5/T0)*((i-nx/2)/sd)*((i-nx/2)/sd) - (0.5/T0)*((j-ny/2)/sd)*((j-ny/2)/sd) - (0.5/T0)*((k-nz/2)/sd)*((k-nz/2)/sd); 

	    // lambda = 0.045
	    arg = -(0.5/T0)*((i-nx/2)/sd)*((i-nx/2)/sd) - (0.5/T0)*((j-ny/2)/sd)*((j-ny/2)/sd) - (0.5/T0)*lambda*lambda*((k-nz/2)/sd)*((k-nz/2)/sd); 	    

	    u_sqr = ux*ux + uy*uy + uz*uz;
	   

	    for (int n = 0; n < Q; n++)
	      {
		
		c_dot_u = c[n][0]*ux + c[n][1]*uy + c[n][2]*uz; 

		c_sqr = (c[n][0]*c[n][0])+(c[n][1]*c[n][1])+(c[n][2]*c[n][2]);

		rho = exp(arg);		

		// (*fIn)[i][j][k][n] = wi[n]*exp(arg)*(chi_sqr + 3.0*c_dot_u + (9.0/2.0)*c_dot_u*c_dot_u - 1.5*u_sqr);

		// (*fIn)[i][j][k][n] = wi[n]*exp(arg)*(chi_sqr + c_dot_u/T0 + (1.0/(2.0*T0*T0))*c_dot_u*c_dot_u - (1.0/(2.0*T0))*u_sqr);

		// (*fIn)[i][j][k][n] = wi[n]*exp(arg)*(chi_sqr + c_dot_u/T0 + (1.0/(2.0*T0*T0)*c_dot_u*c_dot_u - (1.0/(2.0*T0))*u_sqr)) - (wi[n]/3.0)*((c_sqr - 3.0*T0)/(2.0*T0*T0))*(*fIn)[i][j][k][Q+5];

		(*fIn)[i][j][k][n] = wi[n]*rho*(1.0 - ((-5.0*T0 + 5.0*chi_sqr*T0)/(-5.0*T0)) + (1.0 + 0.5*(chi_sqr-1.0)*((c_sqr/T0)-5.0))*(c_dot_u/T0) + (1.0/(2.0*T0*T0))*c_dot_u*c_dot_u - u_sqr/(2.0*T0) + (1.0/(6.0*T0*T0*T0))*c_dot_u*c_dot_u*c_dot_u - (1.0/(2.0*T0*T0))*c_dot_u*u_sqr) - (1.0/3.0)*wi[n]*((c_sqr-3.0*T0)/(2.0*T0*T0))*(*fIn)[i][j][k][Q+5];
	  
	      }

		  // (*fIn)[i][j][k][0] = wi[0]*exp(arg)*((405.0*(-5.0+sqrt(10.0)) + (-8065.0+3014.0*sqrt(10.0))*chi_sqr*T0)/(2.0*(-8059.0+2516.0*sqrt(10.0))) - u_sqr/(2.0*T0));


	    (*fIn)[i][j][k][0] = wi[0]*rho*(1.0 - ((-8065.0+3014.0*sqrt(10.0))*(-5.0*T0+5.0*chi_sqr*T0))/(2.0*(-5.0*T0)*(-5045.0+1507.0*sqrt(10.0))) - u_sqr/(2.0*T0)) - (1.0/3.0)*wi[0]*(-1.5/T0)*(*fIn)[i][j][k][Q+5];

	    // (*fIn)[i][j][k][0] = wi[0]*exp(arg)*((405.0*(-5.0+sqrt(10.0)) + (-8065.0+3014.0*sqrt(10.0))*chi_sqr*T0)/(2.0*(-8059.0+2516.0*sqrt(10.0))) - u_sqr/(2.0*T0)) - (wi[0]/3.0)*((-1.5/T0)*((*fIn)[i][j][k][Q+5]));

	    rho = 0.0;
	    for (int n = 0; n < Q; n++) {
	      rho = rho + (*fIn)[i][j][k][n];
	    }
	    
	    for (int n = 0; n < Q; n++) {

	      c_dot_u = c[n][0]*ux + c[n][1]*uy + c[n][2]*uz; 

	      // (*gIn)[i][j][k][n] = wi[n]*rho*(1.5*chi_sqr*T0 + 0.5*u_sqr + ((5.0/2.0)*chi_sqr*T0 + 0.5*u_sqr)*(3.0*c_dot_u + (9.0/2.0)*c_dot_u*c_dot_u - 1.5*u_sqr));       

	      (*gIn)[i][j][k][n] = wi[n]*rho*(1.5*chi_sqr*T0 + 0.5*u_sqr + ((5.0/2.0)*chi_sqr*T0 + 0.5*u_sqr)*(c_dot_u/T0 + (1.0/(2.0*T0*T0))*c_dot_u*c_dot_u - (1.0/(2.0*T0))*u_sqr));

	      (*fOut)[i][j][k][n]=(*fIn)[i][j][k][n];
	      (*gOut)[i][j][k][n] = (*gIn)[i][j][k][n];
	    }
				  
	  }
	}     
    }
    
}//end make lattice

void write_macro(int ts, int meas_steps, vector<vector<vector<vector <double> > > > * fIn, int nx, int ny, int nz, double sd, double T0) {

  double middlex=nx/2;
  double middley=ny/2;
  double middlez=nz/2;
  float sinv=1.0/sd;

   fstream out;
  char fname[255];

  // DENSITY :

 sprintf(fname,"data/Xrho_t%i.dat",ts);
  out.open(fname, ios::out);
  for(int i=0; i < nx; i++){
    int j=middley;
    int k=middlez;
    out << (i-middlex)*sinv << "\t";
    out << ((*fIn)[i][j][k][Q]) << "\n";
  }  

  out.close();

  sprintf(fname,"data/Zrho_t%i.dat",ts);
  out.open(fname, ios::out);
  for(int k=0; k < nz; k++){
    int i=middlex;
    int j=middley;
    out << (k-middlez)*sinv << "\t";
    out << (*fIn)[i][j][k][Q] << "\n";
  }        
    
  out.close();

  // TEMPERATURE:

  // double T0 = 1.0-sqrt(2.0/5.0);
  double chi_sqr,u_sqr;

  sprintf(fname,"data/Xchi_t%i.dat",ts);
  out.open(fname, ios::out);
  for(int i=0; i < nx; i++){
    int j=middley;
    int k=middlez;
    out << (i-middlex)*sinv << "\t";
    u_sqr = ((*fIn)[i][j][k][Q+1])*((*fIn)[i][j][k][Q+1]) + ((*fIn)[i][j][k][Q+2])*((*fIn)[i][j][k][Q+2]) + ((*fIn)[i][j][k][Q+3])*((*fIn)[i][j][k][Q+3]);

    out << 2.0*((*fIn)[i][j][k][Q+4]-0.5*((*fIn)[i][j][k][Q])*u_sqr)/(3.0*T0*(*fIn)[i][j][k][Q]) << "\n";
  }      
    
  out.close();

  sprintf(fname,"data/Zchi_t%i.dat",ts);
  out.open(fname, ios::out);
  for(int k=0; k < nz; k++){
    int i=middlex;
    int j=middley;
    out << (k-middlez)*sinv << "\t";
    u_sqr = ((*fIn)[i][j][k][Q+1])*((*fIn)[i][j][k][Q+1]) + ((*fIn)[i][j][k][Q+2])*((*fIn)[i][j][k][Q+2]) + ((*fIn)[i][j][k][Q+3])*((*fIn)[i][j][k][Q+3]);

    out << 2.0*((*fIn)[i][j][k][Q+4]-0.5*((*fIn)[i][j][k][Q])*u_sqr)/(3.0*T0*(*fIn)[i][j][k][Q]) << "\n";
  }      
    
  out.close();

  // VELOCITY:

  sprintf(fname,"data/Xux_t%i.dat",ts);
  out.open(fname, ios::out);
  for(int i=0; i < nx; i++){
    int j=middley;
    int k=middlez;
    out << (i-middlex)*sinv << "\t";
    out << ((*fIn)[i][j][k][Q+1]) << "\n";
  }      
    
  out.close(); 

  sprintf(fname,"data/Zuz_t%i.dat",ts);
  out.open(fname, ios::out);
  for(int k=0; k < nz; k++){
    int i=middlex;
    int j=middley;
    out << (k-middlez)*sinv << "\t";
    out << ((*fIn)[i][j][k][Q+3]) << "\n";
  }      
    
  out.close(); 
}

int main(int argc, const char * argv[])
{
    
  // index and size variables
  int nx ,ny, nz,steps,meas_steps, ts;
  double sd, omega_m, omega_E,chi_sqr_init;

    //read in var// 
    // if(argc!=10){
    //     printf("arguments: lattice size x, lattice size y,sd,number of steps, measure every # of steps, omega \n");
    //     exit(1);
    // }
    // nx=atoi(argv[1]);
    // ny=atoi(argv[2]);
    // nz=atoi(argv[3]);
    // sd=atof(argv[4]);
    // steps=atof(argv[5]);
    // meas_steps=atoi(argv[6]);
    // omega_m = atof(argv[7]);
    // omega_E = atof(argv[8]);
    // chi_sqr_init = atof(argv[9]);

  nx = 45;
  ny = 45;
  nz = 1067;

  sd = 10.0;
  steps = 10000;
  meas_steps = 1;
  omega_m = 0.8;
  omega_E = 1.5;

  // double lambda = double(nx)/double(nz);
  // double lambda = 0.0345955;
  // double lambda = 1.0;
  double lambda = 0.035;

    chi_sqr_init = 1.0;
    // sd = sd/sqrt(2.0);
    // omega_m = omega_m*sqrt(2.0);
 
    //time step counter
    ts=0;

    //set number of threads.  If not explicity set will use all avalible.
    //omp_set_num_threads(3);

    vector<vector<vector<vector <double> > > > fIn (nx,vector<vector<vector <double> > > (ny, vector<vector <double> >(nz, vector<double>(Q+6,0.0))));

    vector<vector<vector<vector <double> > > > fOut (nx, vector<vector<vector <double> > >(ny, vector<vector <double> >(nz, vector<double>(Q,0.0))));

    vector<vector<vector<vector <double> > > > gIn (nx,vector<vector<vector <double> > > (ny, vector<vector <double> >(nz, vector<double>(Q,0.0))));

    vector<vector<vector<vector <double> > > > gOut (nx, vector<vector<vector <double> > >(ny, vector<vector <double> >(nz, vector<double>(Q,0.0))));

    double c[Q][3] = {{0.0,0.0,0.0},{1.0,0.0,0.0},{-1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,-1.0,0.0},{0.0,0.0,1.0},{0.0,0.0,-1.0},{1.0,1.0,0.0},{-1.0,-1.0,0.0},{1.0,-1.0,0.0},{-1.0,1.0,0.0},{1.0,0.0,1.0},{-1.0,0.0,-1.0},{1.0,0.0,-1.0},{-1.0,0.0,1.0},{0.0,1.0,1.0},{0.0,-1.0,-1.0},{0.0,1.0,-1.0},{0.0,-1.0,1.0},{1.0,1.0,1.0},{1.0,1.0,-1.0},{1.0,-1.0,1.0},{-1.0,1.0,1.0},{-1.0,-1.0,1.0},{-1.0,1.0,-1.0},{1.0,-1.0,-1.0},{-1.0,-1.0,-1.0},{3.0,0.0,0.0},{-3.0,0.0,0.0},{0.0,3.0,0.0},{0.0,-3.0,0.0},{0.0,0.0,3.0},{0.0,0.0,-3.0},{3.0,3.0,3.0},{3.0,3.0,-3.0},{3.0,-3.0,3.0},{-3.0,3.0,3.0},{-3.0,-3.0,3.0},{-3.0,3.0,-3.0},{3.0,-3.0,-3.0},{-3.0,-3.0,-3.0}};

    // weights 
    double w0 = (2.0/2025.0)*(5045.0-1507.0*sqrt(10.0));
    double w1 = (37.0/(5.0*sqrt(10.0)))-(91.0/40.0);
    double w2 = (1.0/50.0)*(55.0-17.0*sqrt(10.0));
    double w3 = (233.0*sqrt(10.0)-730.0)/1600.0;
    double w4 = (295.0-92.0*sqrt(10.0))/16200.0;
    double w5 = (130.0-41.0*sqrt(10.0))/129600.0;

    double wi[] = {w0,w1,w1,w1,w1,w1,w1,w2,w2,w2,w2,w2,w2,w2,w2,w2,w2,w2,w2,w3,w3,w3,w3,w3,w3,w3,w3,w4,w4,w4,w4,w4,w4,w5,w5,w5,w5,w5,w5,w5,w5};

    double T0 = 1.0 - sqrt(2.0/5.0);

    // double sum1 = 0.0;
    // double sum2 = 0.0;
    // double sum3 = 0.0;
    // double sum4 = 0.0;
    // double sum5 = 0.0;

    // for (int n = 0; n < Q; n++) {
    //   sum1 = sum1 + c[n][0];
    //   sum2 = sum2 + c[n][1];
    //   sum3 = sum3 + c[n][2];
    //   sum4 = sum4 + wi[n];
    //   sum5 = sum5 + wi[n]*(c[n][0]*c[n][0] + c[n][1]*c[n][1] + c[n][2]*c[n][2]);
    // }
    // cout << "1: " << sum1 << "2: " << sum2 << "3: " << sum3 << "4: " << sum4 << "5: " << sum5 << endl;

    // scale omega by the standard deviation
    // omega = omega/sd;
    // cout << "omega: " << omega << endl;
    
    //make initial lattice 
    make_lattice(&fIn,&fOut,&gIn,&gOut,c,wi,nx,ny,nz,sd,lambda,chi_sqr_init,T0);
    
    //do steps
    while ( ts<steps+1) {

      calc_marco(nx,ny,nz,&fIn,&gIn,c,&ts,T0);
	if (ts%meas_steps==0) {
	  write_macro(ts, meas_steps, &fIn, nx, ny, nz, sd, T0);
	  // cout << "writing step " << ts << endl;
	}
	eq(nx,ny,nz,&fIn,&fOut,&gIn,&gOut,c,wi,sd,omega_m,omega_E,lambda,&ts,T0);

        stream(nx,ny,nz,&fIn,&fOut,&gIn,&gOut,c,sd);
	ts = ts + 1;     
    }

    return 0;
}


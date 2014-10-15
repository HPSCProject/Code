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
#include <hdf5.h>
// #include <omp.h>

#define Q 9
#define D 2

using namespace std;

void init_gaussian(vector<vector<vector <double> > > * fIn,vector<vector<vector <double> > > * fOut,vector<vector <double> > * rho, vector<vector <double> > * ux, vector<vector <double> > * uy, int c[Q][D], double wi[Q],double lambda, int nx, int ny, double sd, double T0, double omega) {
  
  double u_sqr, c_dot_u, fEq;
  double x, y;

  double middlex = nx/2;
  double middley = ny/2;

  for (int i = 0; i < nx ; i++) {
    for (int j = 0; j < ny ; j++) {

      u_sqr = (*ux)[i][j]*(*ux)[i][j] + (*uy)[i][j]*(*uy)[i][j]; 

      (*rho)[i][j] = exp( -(0.5/T0)*((i-middlex)/sd)*((i-middlex)/sd) - (0.5/T0)*lambda*lambda*((j-middley)/sd)*((j-middley)/sd) );

      for (int n = 0; n < Q; n++)
	{
	  c_dot_u = c[n][0]*(*ux)[i][j] + c[n][1]*(*uy)[i][j]; 
	 
	  fEq = wi[n]*((*rho)[i][j]);

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

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {

      x = (i-middlex)/sd;
      y = lambda*lambda*(j-middley)/sd;

      //set to zero before summing
      (*rho)[i][j]=0.0;
      (*ux)[i][j]=0.0;
      (*uy)[i][j]=0.0;
	            
      for (int n = 0; n < Q; n++)
	{
	  (*rho)[i][j] = (*rho)[i][j]+(*fIn)[i][j][n]; // rho
	  (*ux)[i][j] = (*ux)[i][j]+(c[n][0])*((*fIn)[i][j][n]); // ux
	  (*uy)[i][j] = (*uy)[i][j]+(c[n][1])*((*fIn)[i][j][n]); // uy
	}

      (*ux)[i][j] = (*ux)[i][j]/(*rho)[i][j];
      (*uy)[i][j] = (*uy)[i][j]/(*rho)[i][j];

      u_sqr = ((*ux)[i][j])*((*ux)[i][j]) + ((*uy)[i][j])*((*uy)[i][j]);

      for (int n = 0; n < Q; n++) {

	c_dot_u = (c[n][0])*((*ux)[i][j]) + (c[n][1])*((*uy)[i][j]);

	c_sqr = (c[n][0])*(c[n][0]) + (c[n][1])*(c[n][1]);

	fEq[n] = wi[n]*((*rho)[i][j])*(1.0 + (1.0/T0)*c_dot_u + (1.0/(2.0*T0*T0))*c_dot_u*c_dot_u - (1.0/(2.0*T0))*u_sqr);

	cdotX = (c[n][0])*x + (c[n][1])*y;
	udotX = x*((*ux)[i][j]) + y*((*uy)[i][j]);

	if (ftrue == 1) {
	  force = -(1.0/T0)*(cdotX-udotX)*fEq[n]/sd;
	  if (i==0 && j==0 && n==0) {
	    cout << "potential is on" << endl;
	  }
	}
	else if (ftrue == 0) {
	  force = 0.0;
	  if (i==0 && j==0 && n==0) {
	    cout << "potential is off" << endl;
	  }
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

  for (int i=0; i < nx; i++) {
    for (int j=0; j < ny; j++) {
      for (int n=0; n < Q; n++) {

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
 
	(*fOut)[i][j][nop[n]] = (*fIn)[in][jn][nop[n]];
      }
    }
  }
}

void write_file(vector<vector <double> > * vec, int dim, int nx, int ny, double sd) {

  double sinv = 1.0/sd;

  hid_t file_id, dataset_id, space_id, property_id;
  herr_t status;

  hsize_t dims[2] = { dim,2 };

  double * towrite = new double[2*dim];

  // write Xrho. THIS ISN'T RIGHT!!!
  int middlex = nx/2;
  int middley = ny/2;
  int j = middley;

  for (int i=0; i<nx; i++) {
    towrite[2*i] = (i-middlex)*sinv;
    towrite[2*i+1] = (*vec)[i][j];
  }

  const string filename = "Xrho_t0.hdf5"; 
  file_id = H5Fcreate( filename.c_str(),H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

  space_id = H5Screate_simple( 2,dims,NULL );

  property_id = H5Pcreate( H5P_DATASET_CREATE );
  status = H5Pset_layout( property_id,H5D_CONTIGUOUS );

  dataset_id = H5Dcreate( file_id,"DATASET",H5T_STD_I32LE,space_id,H5P_DEFAULT,property_id,H5P_DEFAULT );

  status = H5Dwrite( dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,towrite );

  // close up the files
  status = H5Sclose( space_id );
  status = H5Dclose( dataset_id );
  status = H5Fclose( file_id );
  status = H5Pclose( property_id );
}

void write_hdf5(vector<vector <double> > * rho, int nx, int ny, double sd) {

  write_file( rho,nx,nx,ny,sd );
}

void write_gaussian(vector<vector <double> > * rho, vector<vector <double> > * ux, vector<vector <double> > * uy, int nx, int ny, double sd, int ts) {

  fstream out;
  char fname[255];
  float sinv=1.0/sd;
  int middlex = nx/2;
  int middley = ny/2;

  sprintf(fname,"data/Xrho_t%i.dat",ts);
  out.open(fname, ios::out);
  for(int i=0; i < nx; i++){
    int j=middley;
    out << (i-middlex)*sinv << "\t";
    out << (*rho)[i][j] << "\n";
  }  
        
  out.close(); 

  sprintf(fname,"data/Yrho_t%i.dat",ts);
  out.open(fname, ios::out);
  for(int j=0; j < ny; j++){
    int i=middlex;
    out << (j-middley)*sinv << "\t";
    out << (*rho)[i][j] << "\n";
  }  
        
  out.close(); 

  sprintf(fname,"data/Xux_t%i.dat",ts);
  out.open(fname, ios::out);
  for(int i=0; i < nx; i++){
    int j=middley;
    out << (i-middlex)*sinv << "\t";
    out << (*ux)[i][j] << "\n";
  }  
        
  out.close(); 

  sprintf(fname,"data/Yuy_t%i.dat",ts);
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
  int nx = 501;
  int ny = 501;
  int nsteps = 1000;
  double sd = 50.0;
  double omega = 1.0;

  double dtreal = 12; // micro-seconds
  double wxreal = 0.000785398; // micro-Hz
  int dtsim = round( dtreal*wxreal*sd );
  // double lambda = 0.0345955;

  // double lambda = double(nx)/double(ny);
  // double lambda = 0.035;
  double lambda = 1.0;
 
  cout << "dtsim = " << dtsim << endl;

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

  for (int ts=0; ts < nsteps; ts++) {

    if (ts == t_off) {
      ftrue = 0;
    }
    if (ts == t_on) {
      ftrue = 1;
    }

    eq_and_stream(&fIn,&fOut,&rho,&ux,&uy,c,wi,nop,lambda,nx,ny,T0,omega,sd,ftrue);

    fIn = fOut;

    // write_gaussian(&rho,&ux,&uy,nx,ny,sd,ts);
    write_hdf5( &rho,nx,ny,sd );
  }

  return 0;
}

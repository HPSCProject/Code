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

#define Q 9 // number of velocities that connect one lattice node to another. Shown below, labeled with their index in array "c" from left to right:
// *****   
//    \  |  /   c[6] c[2] c[5]
//   --  .  --  c[3] c[0] c[1]
//    /  |  \   c[7] c[4] c[8]

#define D 2 // number of spatial dimensions

// *********************
// THE BASIC IDEA:
// In Init_Gaussian:
// We decide that we want our initial density distribution (rho should be a gaussian with standard deviation sd. The density distribution gives the total number of particles at a point in space. We then use this density distribution to calculate the particle distribution function (fIn), which gives the number of of particles at a particular point in space moving at a particular speed (the "c" array, as defined above). Then we copy the particle distribution function (fIn) into the temporary memory storage (fOut).
// In Eq_And_Stream:
// We now use the particle distribution function (fIn) to calculate the macroscopic variable of fluid density and velocity in each spatial dimension (rho, ux, and uy). Once we have calculated these, we can compute the "equilibrium distribution" (fEq) which approximates the collisions between particles in our fluid and defines how the fluid relaxes to some equilibrium state in some characteristic relaxation time (omega).
// We then calculate how we will move particles to the next timestep. We calculate the "next neighbors" (in, jn) for each particle distribution (fIn[i][j][n]) by calculating which node we would be at next if we moved one lattice unit in the direction given by our velocity (c[n][0] and c[n][1] are the x and y components). The modular arithmetic performs periodic boundaries, so that a particle that leaves one side of the simulation enters the other side.
// Now just write your data, and continue on in time calculating density and fluid velocity from the particles distribution function and propogating the particles around the simulation like we just did!

using namespace std;

void init_gaussian(vector<vector<vector <double> > > * fIn,vector<vector<vector <double> > > * fOut,vector<vector <double> > * rho, vector<vector <double> > * ux, vector<vector <double> > * uy, int c[Q][D], double wi[Q],double lambda, int nx, int ny, double sd, double T0, double omega) {
  
  double u_sqr, c_dot_u, fEq;
  double x, y;

  // index values corresponding to the middle of the simulation
  double middlex = nx/2;
  double middley = ny/2;

  for (int i = 0; i < nx ; i++) {
    for (int j = 0; j < ny ; j++) { 

      // define the initial density to be a gaussian centered at zero
      (*rho)[i][j] = exp( -(0.5/T0)*((i-middlex)/sd)*((i-middlex)/sd) - (0.5/T0)*lambda*lambda*((j-middley)/sd)*((j-middley)/sd) );

      for (int n = 0; n < Q; n++)
	{ 
	  // calculate what the particle distribution function should be to produce the initial density function rho that we want.
	  fEq = wi[n]*((*rho)[i][j]);

	  // copy this into both current and temporary memory storages
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

  // indices for the middle of the spatial simulation volume
  double middlex = nx/2;
  double middley = ny/2;

  double c_sqr;

  double cdotX,udotX;

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {

      // used to calculate the forcing term, which is proportional to the distance of a particle from the center of the simulation
      x = (i-middlex)/sd;
      y = lambda*lambda*(j-middley)/sd;

      //set to zero before summing
      (*rho)[i][j]=0.0;
      (*ux)[i][j]=0.0;
      (*uy)[i][j]=0.0;
	            
      for (int n = 0; n < Q; n++)
	{
	  // these are the sums that give us the density and fluid velocity from the particle distribution function
	  (*rho)[i][j] = (*rho)[i][j]+(*fIn)[i][j][n]; // rho
	  (*ux)[i][j] = (*ux)[i][j]+(c[n][0])*((*fIn)[i][j][n]); // ux
	  (*uy)[i][j] = (*uy)[i][j]+(c[n][1])*((*fIn)[i][j][n]); // uy
	}

      (*ux)[i][j] = (*ux)[i][j]/(*rho)[i][j];
      (*uy)[i][j] = (*uy)[i][j]/(*rho)[i][j];

      // the square of the fluid velocity
      u_sqr = ((*ux)[i][j])*((*ux)[i][j]) + ((*uy)[i][j])*((*uy)[i][j]);

      for (int n = 0; n < Q; n++) {

	// the component of the macroscopic fluid velocity (ux and uy) along the current lattice velocity (c[n][0] and c[n][1])
	c_dot_u = (c[n][0])*((*ux)[i][j]) + (c[n][1])*((*uy)[i][j]);

	c_sqr = (c[n][0])*(c[n][0]) + (c[n][1])*(c[n][1]);

	// equilibrium distribution function
	fEq[n] = wi[n]*((*rho)[i][j])*(1.0 + (1.0/T0)*c_dot_u + (1.0/(2.0*T0*T0))*c_dot_u*c_dot_u - (1.0/(2.0*T0))*u_sqr);

	// components of your displacement x and y from the center of the simulation along the current lattice velocity
	cdotX = (c[n][0])*x + (c[n][1])*y;
	// components of your displacement x and y from the center of the simulation along the fluid velocity (ux and uy)
	udotX = x*((*ux)[i][j]) + y*((*uy)[i][j]);

	if (ftrue == 1) { // if the trapping potential is turned on, then we need a force term acting on the fluid
	  force = -(1.0/T0)*(cdotX-udotX)*fEq[n]/sd;
	  if (i==0 && j==0 && n==0) {
	    cout << "potential is on" << endl;
	  }
	}
	else if (ftrue == 0) { // if the trapping potential is off, we should have no force term
	  force = 0.0;
	  if (i==0 && j==0 && n==0) {
	    cout << "potential is off" << endl;
	  }
	}
	else {
	  cout << "ftrue should be either zero or one, please and thanks." << endl;
	}

	// compute the new particle distribution function to be partially the current particle distribution, and the other part the equilibrium distribution, depending on the relaxation time omega. Add a force if applicable.
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

	// compute the indices of the next lattice node along your lattice velocity (c[n][0] and c[n][1])
	in = i + int(c[n][0]);
	jn = j + int(c[n][1]);
	
	// implement periodic boundary conditions
	if (in > nx-1 || in < 0) {
	  in = (in+nx)%nx;
	  // in = i; 
	}
	if (jn > ny-1 || jn < 0) {
	  jn = (jn+ny)%ny;
	  // jn = j;
	}
 
	// the particle distribution on the next time step (fOut) results from moving particles all to their "neighbor" node (in,jn)
	(*fOut)[i][j][nop[n]] = (*fIn)[in][jn][nop[n]];
      }
    }
  }
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
  // problem size
  int nx = 501; // lattice size in x direction
  int ny = 501; // lattice size in y direction
  int nsteps = 1000; // number of time steps
  double sd = 50.0; // standard deviation of initial gaussian distribution
  double omega = 1.0; // relaxation time

  // physical constants used to compute whether the potential is on or off
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

  // a numerical constant for our algorithm
  double T0 = 1.0/3.0;

  vector<vector<vector <double> > > fIn (nx, vector<vector <double> >(ny, vector <double>(Q,0.0))); // current memory

  vector<vector<vector <double> > > fOut (nx, vector<vector <double> >(ny, vector <double>(Q,0.0))); // temporary memory

  vector<vector <double> > rho (nx, vector <double> (ny,1.0)); // density
  vector<vector <double> > ux (nx, vector <double> (ny,0.0)); // fluid velocity in x direction
  vector<vector <double> > uy (nx, vector <double> (ny,0.0)); // fluid velocity in y direction

  // lattice velocities, which define how a particle is allowed to propogate from one node to the next
  int c[Q][D] = {{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};

  // "opposite of n" gives the indices of the c[nop] that is opposite from c[n]. Only used in the propogation step.
  int nop[Q] = {0,3,4,1,2,7,8,5,6};

  // numerical weighting factors that (approximately) make up for the fact that the diagonal lattice velocities are longer than the straight lattice velocities.
  double w1 = 4.0/9.0;
  double w2 = 1.0/9.0;
  double w3 = 1.0/36.0;
  double wi[] = {w1,w2,w2,w2,w2,w3,w3,w3,w3}; 

  // initialize with our initial conditions we want to simulate
  init_gaussian(&fIn,&fOut,&rho,&ux,&uy,c,wi,lambda,nx,ny,sd,T0,omega);   

  for (int ts=0; ts < nsteps; ts++) {

    // compute whether the potential should be on or off
    if (ts == t_off) {
      ftrue = 0;
    }
    if (ts == t_on) {
      ftrue = 1;
    }

    // calculate the relaxation of the fluid to equilibrium and propogate particles to the nearest nodes.
    eq_and_stream(&fIn,&fOut,&rho,&ux,&uy,c,wi,nop,lambda,nx,ny,T0,omega,sd,ftrue);

    // copy temporary memory into current memory.
    fIn = fOut;
    
    // write the data
    write_gaussian(&rho,&ux,&uy,nx,ny,sd,ts);
  }

  return 0;
}


//  Edited and expanded by PR, Aug 2013

// An isothermal derivative of the file "7-7_D3Q41_2dist_order3.cpp"

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
// #include <omp.h>

#define Q 41
#define D 3

char * outdir;
// initial gaussian deformations
double defx = 1.1;
double defy = 1.1;
double defz = 1.0;

using namespace std;

void make_lattice(vector<vector<vector<vector <double> > > > *fIn,vector<vector<vector<vector <double> > > > *fOut, double c[Q][D], double wi[Q], int nx, int ny,int nz,double sd,double lambda, double T0){
    
  double c_sqr;
  double c_dot_u=0, u_sqr=0;
  double rho,Teff;

  double middlex = nx/2;
  double middley = ny/2;
  double middlez = nz/2;

  // #pragma omp parallel for
  for (int i = 0; i < nx ; i++)
    {
      for (int j = 0; j < ny ; j++)
	{
	  for (int k = 0; k < nz; k++) {
	    
	    rho = exp( -(0.5/T0)*((i-nx/2)/sd)*((i-nx/2)/sd)*defx - (0.5/T0)*((j-ny/2)/sd)*((j-ny/2)/sd)*defy - (0.5/T0)*lambda*lambda*((k-nz/2)/sd)*((k-nz/2)/sd)*defz )*sqrt(defx*defy*defz); 

	    Teff = sqrt(defx*defy*defz);

	    for (int n = 0; n < Q; n++)
	      {

		c_sqr = (c[n][0]*c[n][0])+(c[n][1]*c[n][1])+(c[n][2]*c[n][2]);	

		(*fIn)[i][j][k][n] =  wi[n]*rho*( 1.0 + (c_dot_u/T0)*(1.0 + ((Teff-1.0)/(2.0*T0))*(c_sqr - 5.0*T0)) + (c_dot_u*c_dot_u)/(2.0*T0*T0) - u_sqr/(2.0*T0) + ((Teff-1.0)/(2.0*T0))*(c_sqr-3.0*T0) + (c_dot_u*c_dot_u*c_dot_u)/(6.0*T0*T0*T0) - (u_sqr*c_dot_u)/(2.0*T0*T0) ); // newer

		(*fOut)[i][j][k][n] = (*fIn)[i][j][k][n];
	  
	      }
	  }
	}     
    }    
}

void eq(int nx, int ny, int nz, vector<vector<vector<vector <double> > > > *fIn,vector<vector<vector<vector <double> > > > *fOut,vector<vector<vector <double> > > *rho, vector<vector<vector <double> > > *ux,vector<vector<vector <double> > > *uy,vector<vector<vector <double> > > *uz,double c[Q][D], double wi[Q], double sd, double omega, double lambda,int* ts, double T0, double dt, double ca,int ftrue){
   double c_dot_u, u_sqr, c_sqr;

   double fEq[Q];
   double force[Q];

   double cdotX, udotX;

   double x, y, z;
   double middlex = nx/2;
   double middley = ny/2;
   double middlez = nz/2;

   double Ptot=0,Ntot=0;
   double Teff=1.;

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {  
           
	x = (i-middlex)/sd;
	y = (j-middley)/sd;
	z = lambda*lambda*(k-middlez)/sd;

	//set to zero before summing
	(*rho)[i][j][k]=0.0;
	(*ux)[i][j][k]=0.0;
	(*uy)[i][j][k]=0.0;
	(*uz)[i][j][k]=0.0;
	            
	for (int n = 0; n < Q; n++)
	  {

	    (*rho)[i][j][k]=(*rho)[i][j][k]+(*fIn)[i][j][k][n]; // rho
	    (*ux)[i][j][k]=(*ux)[i][j][k]+((c[n][0])*(*fIn)[i][j][k][n]); // ux
	    (*uy)[i][j][k]=(*uy)[i][j][k]+((c[n][1])*(*fIn)[i][j][k][n]); // uy
	    (*uz)[i][j][k]=(*uz)[i][j][k]+((c[n][2])*(*fIn)[i][j][k][n]); // uz

	    Ptot+=(*fIn)[i][j][k][n]*( (c[n][0])*(c[n][0])+(c[n][1])*(c[n][1])+(c[n][2])*(c[n][2]) );
	  }

	//Divide Speeds by rho
	(*ux)[i][j][k]=((*ux)[i][j][k])/(*rho)[i][j][k]; // ux/rho
	(*uy)[i][j][k]=((*uy)[i][j][k])/(*rho)[i][j][k]; // uy/rho
	(*uz)[i][j][k]=((*uz)[i][j][k])/(*rho)[i][j][k]; // uz/rho

	if (ftrue == 1) {

	  (*ux)[i][j][k] = (*ux)[i][j][k] - 0.5*x/sd;
	  (*uy)[i][j][k] = (*uy)[i][j][k] - 0.5*y/sd;
	  (*uz)[i][j][k] = (*uz)[i][j][k] - 0.5*z/sd;
	}

	u_sqr = ((*ux)[i][j][k])*((*ux)[i][j][k]) + ((*uy)[i][j][k])*((*uy)[i][j][k]) + ((*uz)[i][j][k])*((*uz)[i][j][k]);

	Ptot -= ((*rho)[i][j][k])*u_sqr;
	Ntot += (*rho)[i][j][k];
      }
    }
  }

  Teff = Ptot/Ntot/(3.0*T0);

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {

	x = (i-middlex)/sd;
	y = (j-middley)/sd;
	z = lambda*lambda*(k-middlez)/sd;

	u_sqr = ((*ux)[i][j][k])*((*ux)[i][j][k]) + ((*uy)[i][j][k])*((*uy)[i][j][k]) + ((*uz)[i][j][k])*((*uz)[i][j][k]);


	for (int n = 0; n < Q; n++) {

	  c_dot_u = c[n][0]*(*ux)[i][j][k] + c[n][1]*(*uy)[i][j][k] + c[n][2]*(*uz)[i][j][k]; 

	  c_sqr = (c[n][0]*c[n][0])+(c[n][1]*c[n][1])+(c[n][2]*c[n][2]);

	  fEq[n] = wi[n]*((*rho)[i][j][k])*( 1.0 + (c_dot_u/T0)*(1.0 + ((Teff-1.0)/(2.0*T0))*(c_sqr - 5.0*T0)) + (c_dot_u*c_dot_u)/(2.0*T0*T0) - u_sqr/(2.0*T0) + ((Teff-1.0)/(2.0*T0))*(c_sqr-3.0*T0) + (c_dot_u*c_dot_u*c_dot_u)/(6.0*T0*T0*T0) - (u_sqr*c_dot_u)/(2.0*T0*T0) ); // newer

	  cdotX = (c[n][0])*x + (c[n][1])*y + (c[n][2])*z;
	  udotX = x*((*ux)[i][j][k]) + y*((*uy)[i][j][k]) + z*((*uz)[i][j][k]);

	  if (ftrue == 1) {

	    force[n] = -(dt/sd)*ca*ca*(1.0-0.5*omega/(ca*ca*dt))*wi[n]*((*rho)[i][j][k])*(cdotX/T0-udotX/T0)*(1.0 + ((Teff-1.0)/(2.0*T0))*(c_sqr - 4.0*T0)+c_dot_u/T0+ (c_dot_u*c_dot_u)/(2.0*T0*T0) - u_sqr/(2.0*T0));//paul's form
	  }
	  else if (ftrue == 0) {
	    force[n] = 0.0;
	  }
	  else {
	    cout << "ftrue should be either zero or one, please and thanks." << endl;
	  }

	(*fOut)[i][j][k][n] = (1.0-omega)*(*fIn)[i][j][k][n] + omega*(fEq[n]) + force[n];

	}
      }
    }
  }
}

void stream(int nx, int ny, int nz, vector<vector<vector<vector <double> > > > *fIn,vector<vector<vector<vector <double> > > > *fOut, double c[Q][D], double sd){

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
	  }
	  if (jn > ny-1 || jn < 0) {
	    jn =(jn+ny)%ny;
	  } 
	  if (kn > nz-1 || kn < 0) {
	    kn = (kn+nz)%nz;
	  }

	  (*fIn)[in][jn][kn][n] = (*fOut)[i][j][k][n];
	}
      }
    }                    
  }
}//end eq and stream

void get_walltime(double* wcTime) {
  struct timeval tp;
  gettimeofday(&tp, NULL);
  *wcTime = (double)(tp.tv_sec + tp.tv_usec/1000000.0);
}

void write_gaussian(vector<vector<vector <double> > > * rho, vector<vector<vector <double> > > * ux, int nx, int ny, int nz, double sd, int ts) {

  fstream out;
  char fname[255];
  float sinv=1.0/sd;
  int middlex = nx/2;
  int middley = ny/2;
  int middlez = nz/2;

  sprintf(fname,"%s/Xrho_t%i.dat",outdir,ts);
  out.open(fname, ios::out);
  for(int i=0; i < nx; i++){
    int j=middley;
    int k=middlez;
    out << (i-middlex)*sinv << "\t";
    out << (*rho)[i][j][k] << "\n";
  }  
        
  out.close(); 

  sprintf(fname,"%s/Yrho_t%i.dat",outdir,ts);
  out.open(fname, ios::out);
  for(int j=0; j < ny; j++){
    int i=middlex;
    int k=middlez;
    out << (j-middley)*sinv << "\t";
    out << (*rho)[i][j][k] << "\n";
  }  
        
  out.close(); 

  sprintf(fname,"%s/Zrho_t%i.dat",outdir,ts);
  out.open(fname, ios::out);
  for(int k=0; k < nz; k++){
    int i=middlex;
    int j=middley;
    out << (k-middlez)*sinv << "\t";
    out << (*rho)[i][j][k] << "\n";
  }  
        
  out.close(); 

  sprintf(fname,"%s/Xux_t%i.dat",outdir,ts);
  out.open(fname, ios::out);
  for(int i=0; i < nx; i++){
    int j=middley;
    int k=middlez;
    out << (i-middlex)*sinv << "\t";
    out << setprecision(15) << (*ux)[i][j][k] << "\n";
  }  
        
  out.close();  
}

int main(int argc, const char * argv[])
{
    
  // index and size variables
  int meas_steps;
  double omega;

  int nx = atoi( argv[1] );
  int ny = atoi( argv[2] );
  int nz = atoi( argv[3] );
  double sd = atof( argv[4] );
  int steps = atoi( argv[5] ); 

  // nx = 251;
  // ny = 251;
  // nz = 251;

  // sd = 50.0;
  // steps = 10000;
  meas_steps = 10;
  omega = 1.0;

  // STANDARD LATTICE BOLTZMANN
  double dt = 1.0;
  double ca = 1.0; // MUST CHANGE CODE TO MULTIPLY VECTORS BY CA IF CA IS NO LONGER UNITY!

  double lambda = 1.0;

  char cbuffer[1000];
  sprintf(cbuffer,"data");
  outdir = cbuffer;

  vector<vector<vector<vector <double> > > > fIn (nx,vector<vector<vector <double> > > (ny, vector<vector <double> >(nz, vector<double>(Q,0.0))));
  
  vector<vector<vector<vector <double> > > > fOut (nx, vector<vector<vector <double> > >(ny, vector<vector <double> >(nz, vector<double>(Q,0.0))));

  vector<vector<vector <double> > > rho (nx,vector<vector <double> >(ny, vector<double>(nz,0.0)));

  vector<vector<vector <double> > > ux (nx,vector<vector <double> >(ny, vector<double>(nz,0.0)));

  vector<vector<vector <double> > > uy (nx,vector<vector <double> >(ny, vector<double>(nz,0.0)));

  vector<vector<vector <double> > > uz (nx,vector<vector <double> >(ny, vector<double>(nz,0.0)));
  
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
  
  int ftrue = 1;
     
  double S,E;

  get_walltime(&S);

  make_lattice(&fIn,&fOut,c,wi,nx,ny,nz,sd,lambda,T0);
    
  //do steps
  for (int ts=0; ts<steps+1; ts++) {

    eq(nx,ny,nz,&fIn,&fOut,&rho,&ux,&uy,&uz,c,wi,sd,omega,lambda,&ts,T0,dt,ca,ftrue);

    if (ts%meas_steps==0) {
      write_gaussian(&rho,&ux,nx,ny,nz,sd,ts);
    } 
      
    stream(nx,ny,nz,&fIn,&fOut,c,sd);
  }

  get_walltime(&E);

  printf("Problem Size: nx=%d, ny=%d, nz=%d, sd=%f, nsteps=%d",nx,ny,nz,sd,steps);
  printf("Walltime %f s \n",E-S);

  return 0;
}


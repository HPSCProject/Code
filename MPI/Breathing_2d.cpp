#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include <locale>
#include <sstream>
#include <cstring>
#include "functions.h"

using namespace std;

int main(int argc, char* argv[])
{
  const int c[Q][D] = {{0, 0}, {1, 0}, {0, 1}, {-1, 0}, {0, -1}, {1, 1}, {-1, 1}, {-1, -1}, {1, -1}};
  const int nop[Q] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
  const double wi[9] = {W1, W2, W2, W2, W2, W3, W3, W3, W3};

  bool ftrue = true;	// true indicates trapping potential is turned on; false indicates trapping potential is off.

  cout << "dtsim = " << DT_SIM << endl;

  // Declare macroscopic variables for the fluid
  qr_type **rho;  // Fluid density
  qr_type **ux;  // Fluid velocity in the x direction
  qr_type **uy;  // Fluid velocity in the y direction

  // Input & output arrays.
  qr_type ***fIn, ***fOut;

  // MPI Process control variables
  int my_rank,num_procs,n,ierr;
  int myid_grid,nprocs_grid;
  int coord_grid[2];
  int total_points[2]; //3 element array tot keep track of total points
  int proc_dim[2];     //defining number of dimensions per processor
  int periodicity[2];
  int l_st_I,l_st_J,l_en_I,l_en_J,l_sz_I,l_sz_J;
  MPI_Comm grid_comm_world;

  //restrict the number of processors to be either 1,2 or 3 at the most for k dir since
  //lattice size will be bigger for other dimensions.need to employ more resources there

  ierr=MPI_Init(&argc,&argv);
  ierr=MPI_Comm_rank(MPI_COMM_WORLD,&my_rank); //get processor rank
  ierr=MPI_Comm_size(MPI_COMM_WORLD,&num_procs);//get total processor num
  
  if(my_rank == 0){
    total_points[IDIR]=NX; //points in i dir
    total_points[JDIR]=NY; //points in j dir
    //total_points[KDIR]=Q;  //points in k dir

    proc_dim[IDIR] = 0; // to be populated locally by each processor
    proc_dim[JDIR] = 0; 
    //proc_dim[KDIR] = 0;

    periodicity[IDIR] = 1;
    periodicity[JDIR] = 1;
    //cout <<"Num Procs" <<num_procs<<endl;
    // periodicity[KDIR] = 1;
  }

  //Bcast IVC..............Send call for proc 0,recieve for the others 
  ierr= MPI_Bcast(&total_points,2,MPI_INT,0,MPI_COMM_WORLD);
  //cout<<"Total Points"<<total_points[0]<<" "<<total_points[1]<<" "<<total_points[2]<<endl;
  ierr= MPI_Bcast(&proc_dim,2,MPI_INT,0,MPI_COMM_WORLD);
  //cout<<"Total Points"<<total_points[0]<<" "<<total_points[1]<<" "<<total_points[2]<<endl;
  ierr= MPI_Bcast(&periodicity,2,MPI_INT,0,MPI_COMM_WORLD);
  ierr= MPI_Dims_create(num_procs,2,proc_dim);
  // cout<<"Proc in each dimension"<<proc_dim[0]<<" "<<proc_dim[1]<<endl;


  //Determine 2d cartaesian topology for procs 
  ierr = MPI_Cart_create(MPI_COMM_WORLD,2,proc_dim,periodicity,1,&grid_comm_world);
  ierr = MPI_Comm_rank(grid_comm_world,&myid_grid);//gives processor id in grid
  ierr = MPI_Comm_size(grid_comm_world,&nprocs_grid);
  ierr = MPI_Cart_coords(grid_comm_world,myid_grid,2,coord_grid);
  cout<<"Proc "<<myid_grid<<" Coords " <<coord_grid[IDIR]<<" "<<coord_grid[JDIR]<<" "<<endl; 


   //Determine local coordinates
  l_st_I=local_start(coord_grid[IDIR],proc_dim[IDIR],total_points[IDIR]);
  l_en_I=local_end(coord_grid[IDIR],proc_dim[IDIR],total_points[IDIR]);
  l_sz_I=l_en_I-l_st_I;
  l_st_J=local_start(coord_grid[JDIR],proc_dim[JDIR],total_points[JDIR]);
  l_en_J=local_end(coord_grid[JDIR],proc_dim[JDIR],total_points[JDIR]);
  l_sz_J=l_en_J-l_st_J;
  cout<<"Proc "<<myid_grid<<" Local x " <<l_st_I<<" "<<l_en_I<<" Local y "<<l_st_J<<" "<<l_en_J<<endl; 

  // Fill input & output arrays w/ 0.
  fIn = new qr_type**[l_sz_I];
  fOut = new qr_type**[l_sz_I];

  //cout<<"got here 0"<<endl;
  // FIXME.these can be omp'd
  for(int i = 0; i < l_sz_I; ++i)
    {
      fIn[i] = new qr_type*[l_sz_J];
      fOut[i] = new qr_type*[l_sz_J];
      for(int j = 0; j < l_sz_J; ++j)
	{
	  fIn[i][j] = new qr_type[Q];
	  memset(fIn[i][j], qr_type(0.0), F_SIZE);
	  fOut[i][j] = new qr_type[Q];
	  memset(fOut[i][j], qr_type(0.0), F_SIZE);
	}
    }
  //cout<<"got here 1"<<endl;

  // Fill rho, ux, uy.
  rho = new qr_type*[l_sz_I];
  ux = new qr_type*[l_sz_I];
  uy = new qr_type*[l_sz_I];
  // FIXME.these can be omp'd
  for(int i = 0; i < NX; ++i)
    {
      rho[i] = new qr_type[l_sz_J];
      memset(rho[i], qr_type(0.0), TWOD_SIZE);
      ux[i] = new qr_type[l_sz_J];
      memset(ux[i], qr_type(0.0), TWOD_SIZE);
      uy[i] = new qr_type[l_sz_J];
      memset(uy[i], qr_type(0.0), TWOD_SIZE);
    }
  //cout<<"got here 2"<<endl;

  init_gaussian(fIn, fOut, wi,coord_grid[IDIR],proc_dim[IDIR],total_points[IDIR],coord_grid[JDIR],proc_dim[JDIR],total_points[JDIR]);

  fstream out;
  char fname[255];

  if(my_rank == 2){  
    // sprintf(fname, "data/Xrho_t%i.dat", ts);
  out.open("test_f", ios::out);
  for (int i = 0; i < l_sz_I; ++i)
  {
    for (int j = 0; j < l_sz_J; ++j)
    {
      for (int n = 0; n < Q; ++n)
      {
        out << fIn[i][j][n] << endl;
      }
    }
  }
}
 
  /*

  if(my_rank == 0) //condition for master,I don't think that is necessary
  else if (my_rank%2==1) //processors for init conditions
    {
      init_gaussian(fIn, fOut, rho, ux, uy, wi,l_st_idx,l_ed_idx);
      // init gaussian will initially initialize fin and rho per row,
      // to be changed later to larger array size
        for (size_t ts = 0; ts < N_STEPS; ++ts)
     {
    if (ts == T_ON)
    {
      ftrue = true;
    }
    else if (ts == T_OFF)
    {
      ftrue = false;
    }
    eq_and_stream(fIn, rho, ux, uy, c, wi, ftrue);
    //can we possiibly eliminate this?
    memcpy(fIn,fOut,sizeof(fIn));
    }
  else if (my_rank%2== 0) //write processors
    {
      //place MPI wait statement
      //use blocking receive.each processor will receive 
      //chunk per timesteop,or every step we wish to write
      write_gaussian(rho, ux, uy, ts);
    }
    } */
  ierr = MPI_Finalize();
  return 0;
}

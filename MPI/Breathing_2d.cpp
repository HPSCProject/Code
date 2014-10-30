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

  // Declare the macroscopic variables for the fluid
  qr_type rho[NX][NY];	// Fluid density
  qr_type ux[NX][NY];	// Fluid velocity in the x direction
  qr_type uy[NX][NY];	// Fluid velocity in the y direction

  // Input & output arrays.
  qr_type fIn[NX][NY][Q], fOut[NX][NY][Q];

  // MPI Process control variables
  int my_rank,num_procs,n,ierr;
  int l_arr_sz_q,l_st_idx_q,l_ed_idx_q;
  int l_arr_sz,l_st_idx,l_ed_idx;

  int myid_grid,nprocs_grid;
  int coord_grid[3];

  int total_points[3]; //3 element array tot keep track of total points
  int proc_dim[3];     //defining number of dimensions per processor
  int periodicity[3];

  MPI_Comm grid_comm_world;

  ierr=MPI_Init(&argc,&argv);
  ierr=MPI_Comm_rank(MPI_COMM_WORLD,&my_rank); //get processor rank
  ierr=MPI_Comm_size(MPI_COMM_WORLD,&num_procs);//get total processor num

  if(my_rank == 0){
    total_points[IDIR]=NX; //points in i dir
    total_points[JDIR]=NY; //points in j dir
    total_points[KDIR]=Q;  //points in k dir

    proc_dim[IDIR] = 0; // to be populated locally by each processor
    proc_dim[JDIR] = 0; 
    proc_dim[KDIR] = 0;

    periodicity[0] = 1;
    periodicity[1] = 1;
    periodicity[2] = 1;
  }
  //Send call for proc 0,recieve for the others 
  ierr= MPI_Bcast(&total_points,3,MPI_INT,0,MPI_COMM_WORLD);
  //cout<<"Total Points"<<total_points[0]<<" "<<total_points[1]<<" "<<total_points[2]<<endl;

  ierr= MPI_Bcast(&proc_dim,3,MPI_INT,0,MPI_COMM_WORLD);
  //cout<<"Total Points"<<total_points[0]<<" "<<total_points[1]<<" "<<total_points[2]<<endl;

  //if(my_rank==0){
  ierr= MPI_Dims_create(num_procs,3,proc_dim);
  // cout<<"Proc in each dimension"<<proc_dim[0]<<" "<<proc_dim[1]<<" "<<proc_dim[2]<<endl;

  ierr = MPI_Cart_create(MPI_COMM_WORLD,3,proc_dim,periodicity,1,&grid_comm_world);
  ierr = MPI_Comm_rank(grid_comm_world,&myid_grid);
  ierr = MPI_Comm_size(grid_comm_world,&nprocs_grid);
  ierr = MPI_Cart_coords(grid_comm_world,myid_grid,3,coord_grid);

  cout<<"Proc "<<myid_grid<<" Coords " <<coord_grid[0]<<" "<<coord_grid[1]<<" "<<coord_grid[2]<<endl; 
 
  /*
  //these may not be neccessary
  l_arr_sz_q = (NX*NY*Q)/num_procs;
  l_st_idx_q = my_rank*l_arr_sz_q;
  l_ed_idx_q = l_st_idx_q+l_arr_sz_q;

  l_arr_sz = (NX*NY)/num_procs;
  l_st_idx = my_rank*l_arr_sz;
  l_ed_idx = l_st_idx+l_arr_sz;

  // Fill input & output arrays w/ 0
  // ASHAHID: I think these could be parallelized with OMP
  memset(fIn,0.0,sizeof(qr_type)*l_arr_sz);
  memset(fOut,0.0,sizeof(qr_type)*l_arr_sz);

  // Fill rho, ux, uy.
  // ASHAHID: I think these could be parallelized with OMP
  memset(rho, 1.0, sizeof(qr_type)*l_arr_sz);
  memset(ux, 0.0, sizeof(qr_type)*l_arr_sz);
  memset(uy, 0.0, sizeof(qr_type)*l_arr_sz);

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

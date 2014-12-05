#include <mpi.h>
#include <omp.h>
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
    double time1, time2;
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
    int src,dst;     //For neighbor identification
    int bufsize[2];  //For 2D
    int l_st_I,l_st_J,l_en_I,l_en_J,l_sz_I,l_sz_J;
    MPI_Comm grid_comm_world;
    qr_type *bufTx; //buffer for sending halo points
    qr_type *bufRx; //buffer for receiveing halo points
    qr_type *temp; //buffer for storing corner points
    double send[2];
    send[0] = 2;
    send[1] = 4;
    double receive[2];
    
    //restrict the number of processors to be either 1,2 or 3 at the most for k dir since
    //lattice size will be bigger for other dimensions.need to employ more resources there
    MPI_Request request;
    MPI_Status status;
    int provided;
    
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
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
        time1 = get_walltime();
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
    int local_dim[2];
    for (int i = 0;i<2;i++){
        local_dim[i]=total_points[i]/proc_dim[i];//no.of pts in 1 dir divided by number of procs in that dir
        // if (coord_grid[i] < (spat_dim[i]%proc_dim[i]))
        //local_dim[i] = local_dim[i]+1;
        // cout<<"Proc "<<myid_grid<<"local " <<local_dim[i]<<endl;
    }
    
    l_st_I=local_start(coord_grid[IDIR],proc_dim[IDIR],total_points[IDIR]);
    l_en_I=local_end(coord_grid[IDIR],proc_dim[IDIR],total_points[IDIR]);
    // l_sz_I=(l_en_I-l_st_I)+1;
    l_sz_I = local_dim[IDIR];
    l_st_J=local_start(coord_grid[JDIR],proc_dim[JDIR],total_points[JDIR]);
    l_en_J=local_end(coord_grid[JDIR],proc_dim[JDIR],total_points[JDIR]);
    // l_sz_J=(l_en_J-l_st_J)+1;
    l_sz_J = local_dim[JDIR];
    //  cout<<"Proc "<<myid_grid<<" Local x " <<l_st_I<<" "<<l_en_I<<" Local y "<<l_st_J<<" "<<l_en_J<<endl;
    // cout<<"Proc "<<myid_grid<<" Local x " <<l_sz_I<<" Local y "<<l_sz_J<<endl;
    // Fill input & output arrays w/ 0.
    fIn = new qr_type**[l_sz_I+2];
    fOut = new qr_type**[l_sz_I+2];
    
    // FIXME.these can be omp'd
    for(int i = 0; i <= l_sz_I+1; ++i)
    {
        fIn[i] = new qr_type*[l_sz_J+2];
        fOut[i] = new qr_type*[l_sz_J+2];
        for(int j = 0; j <= l_sz_J+1; ++j)
        {
            fIn[i][j] = new qr_type[Q];
            memset(fIn[i][j], qr_type(0.0), F_SIZE);
            fOut[i][j] = new qr_type[Q];
            memset(fOut[i][j], qr_type(0.0), F_SIZE);
        }
    }
    // Fill rho, ux, uy.
    rho = new qr_type*[l_sz_I+2];
    ux = new qr_type*[l_sz_I+2];
    uy = new qr_type*[l_sz_I+2];
    // FIXME.these can be omp'd
    
    bufTx = new qr_type[(2+2+l_sz_I)*Q]; //assumption:lattice is perfect square/cube
    bufRx = new qr_type[(2+2+l_sz_I)*Q];
    temp  = new qr_type[4*Q]; //for storing cornern points
    //setting buffer length
    bufsize[IDIR]=(2+2+l_sz_J)*Q; //x-plane another 2 for storing corner points
    bufsize[JDIR]=(2+2+l_sz_I)*Q; //y-plane
    //temp =
    
    for(int i = 0; i <= l_sz_I+1; ++i)
    {
        rho[i] = new qr_type[l_sz_J+2];
        memset(rho[i], qr_type(0.0), sizeof(qr_type)*(l_sz_J+2));
        ux[i] = new qr_type[l_sz_J+2];
        memset(ux[i], qr_type(0.0), sizeof(qr_type)*(l_sz_J+2));
        uy[i] = new qr_type[l_sz_J+2];
        memset(uy[i], qr_type(0.0),sizeof(qr_type)*(l_sz_J+2));
    }
    
    // init_gaussian(fIn, fOut, wi,coord_grid[IDIR],proc_dim[IDIR],total_points[IDIR],coord_grid[JDIR],proc_dim[JDIR],total_points[JDIR],my_rank,l_sz_I,l_sz_J);
    init_gaussian(fIn, fOut, wi,coord_grid[JDIR],proc_dim[JDIR],total_points[JDIR],coord_grid[IDIR],proc_dim[IDIR],total_points[IDIR],my_rank,l_sz_I,l_sz_J);
    fstream out;
    if(myid_grid == 0) out.open("data/pinit", ios::out);
    for (int i = 1; i <= l_sz_I; ++i)
    {
        for (int j = 1; j <= l_sz_J; ++j)
        {
            for (int n = 0; n < Q; ++n)
            {
                if(myid_grid == 0) out << fIn[i][j][n] << endl;
            }
        }
    }
    out.close();
    for (int ts = 0; ts < N_STEPS; ++ts) ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////NSTEPS
    {
        if (ts == T_ON)
        {
            ftrue = true;
            // cout<<" Potential is on" <<endl;
        }
        else if (ts == T_OFF)
        {
            ftrue = false;
            //cout<<" Potential is off" <<endl;
        }
        eq_and_stream(fIn, fOut, rho, ux, uy, c, wi, nop, ftrue,coord_grid[JDIR],proc_dim[JDIR],total_points[JDIR],coord_grid[IDIR],proc_dim[IDIR],total_points[IDIR],my_rank,l_sz_I,l_sz_J);
        for(int i = 0; i <=l_sz_I+1; ++i)
        {
            for(int j = 0; j <= l_sz_J+1; ++j)
                memcpy(fIn[i][j], fOut[i][j], F_SIZE);
        }
        
        if(myid_grid == 0)  out.close();
        
        //copy the points to vertical neighbours
        
        for (int disp=-1;disp <2; disp=disp+2)
        {
            //Determine neighbouring processors
            ierr = MPI_Cart_shift(grid_comm_world,1,disp,&src,&dst);
            //if(myid_grid==1) cout<<"rgt nbr"<<src<<" lft  nbr"<<dst<<endl;
            if(src !=MPI_PROC_NULL){
                ierr = MPI_Irecv(&bufRx[0],2*Q,MPI_DOUBLE_PRECISION,src,0,grid_comm_world,&request);
            }
            if(dst !=MPI_PROC_NULL){
                cpy_send_cr_pts(fIn,l_sz_I,bufTx,disp,myid_grid);
                ierr = MPI_Send(&bufTx[0],2*Q,MPI_DOUBLE_PRECISION,dst,0,grid_comm_world);
            }
            if(src !=MPI_PROC_NULL){
                ierr = MPI_Wait(&request,&status);
                cpy_receive_cr_pts(temp,bufRx,disp);
                
            }
        } //end for disp
        
//        if(my_rank==0){
//            for(int i=0;i<Q;i++)
//                //cout << "Tx Pt" << fIn[5][5][i]<<endl;
//        }
//        
//        if(my_rank==1){
//            for(int i=0;i<Q;i++)
//               // cout << "Rx early Pt" << temp[(Q)+i]<<endl;
//        }
        
        for (int disp=-1;disp <2; disp=disp+2)
        {
            for (int dir=0;dir<2;dir++)
            {
                //Determine neighbouring processors
                ierr = MPI_Cart_shift(grid_comm_world,dir,disp,&src,&dst);
                // if(myid_grid==0) cout<<"rgt nbr"<<src<<" lft  nbr"<<dst<<endl;
                if(src !=MPI_PROC_NULL){
                    ierr = MPI_Irecv(&bufRx[0],bufsize[dir],MPI_DOUBLE_PRECISION,src,0,grid_comm_world,&request);
                }
                if(dst !=MPI_PROC_NULL){
                    cpy_send_buf(fIn,l_sz_I,bufTx,dir,disp,temp,myid_grid);
                    ierr = MPI_Send(&bufTx[0],bufsize[dir],MPI_DOUBLE_PRECISION,dst,0,grid_comm_world);
                    //cout << "Proc: " << myid_grid <<" sent to " <<dst<<endl;
                }
                if(src !=MPI_PROC_NULL){
                    ierr = MPI_Wait(&request,&status);
                    cpy_receive_buf(fIn,l_sz_I,bufRx,dir,disp,temp,myid_grid);
                }
            }// end for dir
        }//end for disp
        
        if(my_rank==3){
            for(int i=0;i<Q;i++)
                cout << "End Rx Pt" << fIn[0][0][i]<<endl;
        }
        
        if(myid_grid == 0) out.open("data/aeq1", ios::out);
        for (int i = 1; i <= l_sz_I; ++i)
        {
            for (int j = 1; j <= l_sz_J; ++j)
            {
                for (int n = 0; n < Q; ++n)
                {
                    if(myid_grid == 0) out << fIn[i][j][n] << endl;
                }
            }
        }
        out.close();
        
        
        /*
         if(myid_grid == 1) out.open("data/bdry6row1", ios::out);
         for (int i = 1; i <= l_sz_I; ++i)
         {
         for (int n = 0; n < Q; ++n)
         {
         if(myid_grid == 1) out << fIn[1][i][n] << endl;
         }
         }
         if(myid_grid == 1) out.close();
         
         if(myid_grid == 2) out.open("data/bdry6col1", ios::out);
         for (int i = 1; i <= l_sz_I; ++i)
         {
         for (int n = 0; n < Q; ++n)
         {
         if(myid_grid == 2) out << fIn[i][1][n] << endl;
         }
         }
         if(myid_grid == 2) out.close();
         
         
         if(myid_grid == 0) out.open("data/bdry6row", ios::out);
         for (int i = 1; i <= l_sz_I; ++i)
         {
         for (int n = 0; n < Q; ++n)
         {
         if(myid_grid == 0) out << fIn[6][i][n] << endl;
         }
         }
         if(myid_grid == 0) out.close();
         
         
         if(myid_grid == 0) out.open("data/bdry6col", ios::out);
         for (int i = 1; i <= l_sz_I; ++i)
         {
         for (int n = 0; n < Q; ++n)
         {
         if(myid_grid == 0) out << fIn[i][6][n] << endl;
         }
         }
         if(myid_grid == 0) out.close(); */
        
        
        eq_and_stream(fIn, fOut, rho, ux, uy, c, wi, nop, ftrue,coord_grid[JDIR],proc_dim[JDIR],total_points[JDIR],coord_grid[IDIR],proc_dim[IDIR],total_points[IDIR],my_rank,l_sz_I,l_sz_J);
        for(int i = 0; i <=l_sz_I+1; ++i)
        {
            for(int j = 0; j <= l_sz_J+1; ++j)
            {
                memcpy(fIn[i][j], fOut[i][j], F_SIZE);
            }
        }
        
        
        if(myid_grid == 0) out.open("data/aeq", ios::out);
        for (int i = 1; i <= l_sz_I; ++i)
        {
            for (int j = 1; j <= l_sz_J; ++j)
            {
                for (int n = 0; n < Q; ++n)
                {
                    if(myid_grid == 0) out << fIn[i][j][n] << endl;
                }
            }
        }
        if(myid_grid == 0) out.close();
        
    }//end for ts
    if(rank == 0) {
        time2 = get_walltime() - time1;
        cout << "Total Time: " << time2 << endl;
    }
    
    ierr = MPI_Finalize();
    return 0;
}



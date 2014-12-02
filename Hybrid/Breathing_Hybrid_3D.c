#include <stdlib.h>
#include <omp.h>
#include <mpi.h>
#include <math.h>

void init(double* rho, double* ux, double* uy, double* df, int* local_dims, int* spacial_dims, int* local_index, double lambda, double sd, double* weight);
//T0

int main(int argc, char** argv)
{
    //init mpiworld
    int rank, size, i, j, k;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // reading command line args
    int nx, ny, tsteps, write_steps;
    double sd, lambda;
    
    for(i = 0; i < argc; ++i)
    {
        if(strcmp(argv[i],"-nx") == 0)
            nx = atoi(argv[i+1]);
        if(strcmp(argv[i],"-ny") == 0)
            ny = atoi(argv[i+1]);
        if(strcmp(argv[i],"-tsteps") == 0)
            tsteps = atoi(argv[i+1]);
        if(strcmp(argv[i],"-write_steps") == 0)
            write_steps = atoi(argv[i+1]);
        if(strcmp(argv[i],"-sd") == 0)
            sd = atof(argv[i+1]);
        if(strcmp(argv[i],"-lambda") == 0)
            lambda = atof(argv[i+1]);
    }
    
    // setup cartesian topology
    int spacial_dims[] = {nx, ny};
    int block_dims[] = {0,0};
    int periodicity[] = {1,1};
    MPI_Comm MPI_COMM_GRID;
    MPI_Dims_create(size, 2, block_dims); // determines best topology
    MPI_Cart_create(MPI_COMM_WORLD, 2, block_dims, periodicity, 0, &MPI_COMM_GRID); // new communicator with grid info attached
    int size_grid, rank_grid;
    MPI_Comm_size(MPI_COMM_GRID, &size_grid);
    MPI_Comm_rank(MPI_COMM_GRID, &rank_grid);
    
    int proc_coord[] = {0,0};
    MPI_Cart_coords(MPI_COMM_GRID, rank_grid, 2, proc_coord);
    
    // grid sizes
    int local_dims[2], local_index[2];
    for(i = 0; i < 2; ++i)
    {
        local_dims[i] = spacial_dims[i]/block_dims[i];
        local_index[i] = local_dims[i]*proc_coord[i];
        if(proc_coord[i] < (spacial_dims[i] % block_dims[i])) { //deals with division incompatibility
            local_dims[i] += 1;
            local_index[i] +=(spacial_dims[i] % block_dims[i]) - proc_coord[i];
        }
    }
    
    
    
    
    // setup initial conditions
    // aligned to 64 bit helps with vectorization
    double* rho = (double *) aligned_alloc(64, sizeof(double)*(local_dims[0]+2)*(local_dims[1]+2));
    double* ux = (double *) aligned_alloc(64, sizeof(double)*(local_dims[0]+2)*(local_dims[1]+2));
    double* uy = (double *) aligned_alloc(64, sizeof(double)*(local_dims[0]+2)*(local_dims[1]+2));
    double* df = (double *) aligned_alloc(64, sizeof(double)*(local_dims[0]+2)*(local_dims[1]+2)*41);
    
    double stencil[] = {{0,0},{1,0},{-1,0},{0,1},{0,-1},{1,1},{-1,1},{1,-1},{-1,-1}};
    double weight[] = {4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36};
    
    init(rho, ux, uy, df, local_dims, spacial_dims, local_index, lambda, sd, weight);
    
    printf("Value at index (%i,%i) is %f \n",local_index[0],local_index[1], rho[0] );

    
}


void init(double* rho, double* ux, double* uy, double* df, int* local_dims, int* spacial_dims, int* local_index, double lambda, double sd, double* weight)
{
    int imax = local_dims[0], jmax = local_dims[1];
    int nx = spacial_dims[0], ny = spacial_dims[1];
    double xbar = nx/2, ybar = ny/2;
    double c0, c1, c2, c3, val1, val2;
    c0 = -1/(2*T0*sd*sd);
    c1 = lambda*lambda;
    
    for(int j = 0; j < jmax; ++j) {
        c2 = (local_index[1]+j-jbar) * (local_index[1]+j-jbar);
        for(int i = 0; i < imax; ++i) {
            c3 = (local_index[0] + i - ibar) * (local_index[0] + i - ibar);
            val = exp(c0 * ( c3 + (c1 * c2)));
            rho[(i+1)+((j+1)*(imax+2))] = val;
            for(k = 0; k < 9; ++k) {
                df[((i+1)+((j+1)*(imax+2)))*9 + k] = weight[k] * val;
            }
        }
    }
    
    // ux uy do not need to be initialized since all values are reset at beginning of eq
    
}














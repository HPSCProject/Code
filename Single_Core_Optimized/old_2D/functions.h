#ifndef __FUNCTIONS_H_
#define __FUNCTIONS_H_
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <sys/time.h>

// PREPROCESSOR MACROS
#define ROUND(a) 	(round(a))
#define SIZE_OF(b) 	(sizeof(b))

// CONSTANT, PREPROCESSOR DEFINTIONS
typedef double qr_type;

#define NX 			int(251)						// Lattice size in x direction.
#define NY 			int(251) 						// Lattice size in y direction.
#define MIDDLE_X 	int((NX / 2))					// Middle x index value of the stream, rounded down.
#define MIDDLE_Y 	int((NY / 2))					// Middle y index value of the stream, rounded down.

#define N_STEPS		int(250) 						// Number of time steps.
#define SD 			qr_type(50.0)					// Standard deviation of initial guassian distibution.
#define OMEGA 		qr_type(1.0)					// Relaxation time.
#define SIN_V		qr_type((1.0 / SD))				// Reciprocal of standard deviation.

#define DT_REAL		int(12)							// Description (in microseconds)
#define WX_REAL		qr_type(0.000785398)			// Description (micro-Hz)
#define DT_SIM		int(ROUND(DT_REAL*WX_REAL*SD))	// Description
#define LAMBDA 		qr_type(1.0)					// Description

#define T_OFF 		int(10)							// Time step to turn potential off.
#define T_ON 		int((T_OFF + DT_SIM))			// Time step to turn potential back on.

#define T0 			qr_type((1.0 / 3.0))			// Description

#define W1 			qr_type((4.0 / 9.0))			// Description
#define W2 			qr_type((1.0 / 9.0))			// Description
#define W3 			qr_type((1.0 / 36.0))			// Description

#define Q 			int(9)							// Number of velocities that connect one lattice node to another. 
													// Shown below, labeled with their index in array "c" from left to right:
													//    \  |  /   c[6] c[2] c[5]
													//   --  .  --  c[3] c[0] c[1]
													//    /  |  \   c[7] c[4] c[8]
#define D 			int(2) 							// Number of spatial dimensions.
#define F_SIZE		int((SIZE_OF(qr_type)*Q))		// Size of the innermost fIn/fOut array.
#define TWOD_SIZE	int((SIZE_OF(qr_type)*NY))		// Size of the innermost rho/ux/uy array.
#define TWOD_TOTAL	int((SIZE_OF(qr_type)*NX*NY))	// Size of the entire rho/ux/uy array.

using namespace std;

/* THE BASIC IDEA:
* In Init_Gaussian:
* We decide that we want our initial density distribution (rho should be a gaussian with standard deviation sd. The density distribution gives the total number of particles at a point in space. We then use this density distribution to calculate the particle distribution function (fIn), which gives the number of of particles at a particular point in space moving at a particular speed (the "c" array, as defined above). Then we copy the particle distribution function (fIn) into the temporary memory storage (fOut).
* In Eq_And_Stream:
* We now use the particle distribution function (fIn) to calculate the macroscopic variable of fluid density and velocity in each spatial dimension (rho, ux, and uy). Once we have calculated these, we can compute the "equilibrium distribution" (fEq) which approximates the collisions between particles in our fluid and defines how the fluid relaxes to some equilibrium state in some characteristic relaxation time (omega).
* We then calculate how we will move particles to the next timestep. We calculate the "next neighbors" (in, jn) for each particle distribution (fIn[i][j][n]) by calculating which node we would be at next if we moved one lattice unit in the direction given by our velocity (c[n][0] and c[n][1] are the x and y components). The modular arithmetic performs periodic boundaries, so that a particle that leaves one side of the simulation enters the other side.
* Now just write your data, and continue on in time calculating density and fluid velocity from the particles distribution function and propogating the particles around the simulation like we just did!
*/

// Initializes the values in both fIn and fOut to step 0 values.
void init_gaussian(qr_type ***fIn, qr_type ***fOut, const double wi[Q]);

 // Performs one step. Description about what does to each input param?
void eq_and_stream(qr_type ***fIn, qr_type ***fOut, qr_type **rho, qr_type **ux, qr_type **uy, const int c[Q][D], const double wi[Q], const int nop[Q], const bool& ftrue);

 // Writes out rho, ux, and uy to data files.
void write_gaussian(qr_type **rho, qr_type **ux, qr_type **uy, const int& ts);

// Get the current wall time.
void get_walltime(double& wcTime);

#endif
#ifndef __FUNCTIONS_H_
#define __FUNCTIONS_H_
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <sys/time.h>

// PREPROCESSOR MACROS
#define SQRT(a) 	(sqrt(a))
#define SIZE_OF(b) 	(sizeof(b))

// CONSTANT, PREPROCESSOR DEFINTIONS
typedef double qr_type;

#define Q			int(41)														// Number of velocities.
#define D 			int(3)														// Number of spatial dimensions.

#define DEF_X		qr_type(1.1)												// Gaussian deformation in x-direction.
#define DEF_Y		qr_type(1.1)												// Gaussian deformation in y-direction.
#define DEF_Z		qr_type(1.0)												// Gaussian deformation in z-direction.
#define T_EFF		qr_type(SQRT(DEF_X * DEF_Y * DEF_Z))						// Description.

#define NX			int(251)													// Lattice size in x-direction.
#define NY			int(251)													// Lattice size in y-direction.
#define NZ			int(251)													// Lattice size in z-direction.
#define MIDDLE_X	int(NX / 2)													// Middle x index value of the stream, rounded down.
#define MIDDLE_Y	int(NY / 2)													// Middle y index value of the stream, rounded down.
#define MIDDLE_Z	int(NZ / 2)													// Middle z index value of the stream, rounded down.

#define SD 			qr_type(50.0)												// Standard deviation of initial guassian distibution.
#define SD2			qr_type(2 * SD)												// Double standard deviation.
#define SIN_V		qr_type(1.0 / SD)											// Reciprocal of standard deviation.
#define STEPS 		int(10000)													// Number of iterations.
#define MEAS_STEPS	int(5) 														// How many steps before a measurement is written out.

#define OMEGA 		qr_type(1.0)												// Relaxation time.
#define LAMBDA 		qr_type(1.0)												// Description.
#define LAMBDA_SQR	qr_type(LAMBDA * LAMBDA)									// Lambda squared.
#define DT 			qr_type(1.0)												// Description.
#define CA 			qr_type(1.0)												// Description.

#define W0 			qr_type((2.0 / 2025.0) * (5045.0 - 1507.0 * SQRT(10.0)))	// First weight.
#define W1 			qr_type((37.0 / (5.0 * SQRT(10.0))) - (91.0 / 40.0))		// Second weight.
#define W2 			qr_type((1.0 / 50.0) * (55.0 - 17.0 * SQRT(10.0)))			// Third weight.
#define W3 			qr_type((233.0 * SQRT(10.0) - 730.0) / 1600.0)				// Fourth weight.
#define W4 			qr_type((295.0 - 92.0 * SQRT(10.0)) / 16200.0)				// Fifth weight.
#define W5 			qr_type((130.0 - 41.0 * SQRT(10.0)) / 129600.0)				// Sixth weight.

#define T0 			qr_type(1.0 - SQRT(2.0 / 5.0))								// Description.
#define T0_H		qr_type(0.5 / T0)											// Description.

#define F_SIZE		int((SIZE_OF(qr_type) * Q))									// Size of the innermost fIn/fOut array.
#define D_SIZE		int((SIZE_OF(qr_type) * NZ))								// Size of the innermost rho/ux/uy array.

#define WRITE_DIR	"data"														// The directory to which the results are written.
#define BUFF_SIZE	255															// Filename buffer size.

using namespace std;

void make_lattice(qr_type**** fIn, qr_type**** fOut, const double c[Q][D], const double wi[Q]);

void eq(qr_type**** fIn, qr_type**** fOut, qr_type*** rho, qr_type*** ux, qr_type*** uy, qr_type*** uz, const double c[Q][D], const double wi[Q], const bool& ftrue);

void stream(qr_type**** fIn, qr_type**** fOut, const double c[Q][D]);

void write_gaussian(qr_type*** rho, qr_type*** ux, const int& ts);

// Get the current wall time.
void get_walltime(double& wcTime);

#endif

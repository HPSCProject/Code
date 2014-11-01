#ifndef __FUNCTIONS_H_
#define __FUNCTIONS_H_
#include <math.h>
#include <fstream>
#include <iostream>

// PREPROCESSOR MACROS
#define ROUND(a) 	(round(a))
#define FLOOR(b) 	(floor(b))
#define SIZE_OF(c) 	(sizeof(c))

// CONSTANT, PREPROCESSOR DEFINTIONS
typedef double qr_type;

#define NX 			10							// Description
#define NY 			10 							// Description
#define MIDDLE_X 	(NX/2.0)					// Description
#define MIDDLE_Y 	(NY/2.0)					// Description
#define FLOOR_X		FLOOR(MIDDLE_X)				// MIDDLE_X rounded down.
#define FLOOR_Y		FLOOR(MIDDLE_Y)				// MIDDLE_Y rounded down.

#define N_STEPS		20 							// Description
#define SD 			50.0						// Description
#define OMEGA 		1.0 						// Description
#define SIN_V		(1.0/SD)					// Description

#define DT_REAL		12							// Description
#define WX_REAL		0.000785398					// Description
#define DT_SIM		ROUND(DT_REAL*WX_REAL*SD)	// Description
#define LAMBDA 		1.0 						// Description

#define T_OFF 		10							// Description
#define T_ON 		(T_OFF+DT_SIM)				// Description

#define T0 			(1.0/3.0)					// Description

#define W1 			(4.0/9.0)					// Description
#define W2 			(1.0/9.0)					// Description
#define W3 			(1.0/36.0)					// Description

#define Q 			9							// Description
#define D 			2 							// Description
#define F_SIZE		(SIZE_OF(qr_type)*NX*NY*Q)	// Size of the fIn/fOut array.
#define TWOD_SIZE		(SIZE_OF(qr_type)*NX*NY)	// Size of the rho/ux/uy.

#define IDIR            0
#define JDIR            1
#define KDIR            2

//Function description,input parameter,output parameters??
void init_gaussian(qr_type fIn[NX*NY*Q], qr_type fOut[NX*NY*Q], qr_type rho[NX*NY], qr_type ux[NX*NY], qr_type uy[NX*NY], const double wi[Q],int l_start_index,int l_stop_index);

 //Function description,input parameter,output parameters??
void eq_and_stream(qr_type fIn[NX*NY*Q], qr_type rho[NX*NY], qr_type ux[NX*NY], qr_type uy[NX*NY], const int c[Q][D], const double wi[Q], const bool& ftrue);

 //Function description,input parameter,output parameters??
void write_gaussian(qr_type rho[NX*NY], qr_type ux[NX*NY], qr_type uy[NX*NY], const size_t& ts);

int local_start(int dir_rank,int dir_num_procs,int dir_t_points);
int local_end(int dir_rank,int dir_num_procs,int dir_t_points);
#endif

#ifndef __FUNCTIONS_H_
#define __FUNCTIONS_H_
#include <math.h>
#include <fstream>

// PREPROCESSOR MACROS
#define ROUND(a) (round(a))

// CONSTANT, PRE_PROCESSOR DEFINTIONS
#define NX 		10							// Description
#define NY 		10 							// Description
#define N_STEPS	20 							// Description
#define SD 		50.0						// Description
#define OMEGA 	1.0 						// Description

#define DT_REAL	12							// Description
#define WX_REAL	0.000785398					// Description
#define DT_SIM	ROUND(DT_REAL*WX_REAL*SD)	// Description
#define LAMBDA 	1.0 						// Description

#define T_OFF 	10							// Description
#define T_ON 	(T_OFF+DT_SIM)				// Description

#define T0 		(1.0/3.0)					// Description

#define W1 		(4.0/9.0)					// Description
#define W2 		(1.0/9.0)					// Description
#define W3 		(1.0/36.0)					// Description

#define Q 		9							// Description
#define D 		2 							// Description

using namespace std;

//Function description,input parameter,output parameters??
void init_gaussian(double fIn[NX*NY*Q], double fOut[NX*NY*Q], double rho[NX*NY], double ux[NX*NY], double uy[NX*NY], const int c[Q][D], const double wi[Q]);

 //Function description,input parameter,output parameters??
void eq_and_stream(double fIn[NX*NY*Q], double fOut[NX*NY*Q], double rho[NX*NY], double ux[NX*NY], double uy[NX*NY], const int c[Q][D], const double wi[Q],
				   const int nop[Q], const bool& ftrue);

 //Function description,input parameter,output parameters??
void write_gaussian(double rho[NX*NY], double ux[NX*NY], double uy[NX*NY], const size_t& ts);

#endif
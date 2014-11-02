#include "functions.h"
using namespace std;

void init_gaussian(qr_type ***fIn, qr_type ***fOut, const double wi[Q],int x_rank,int x_num_procs,int x_t_points,int y_rank,int y_num_procs,int y_t_points)
{
  qr_type fEq;
  int global_x=local_start(x_rank,x_num_procs,x_t_points);//FIXME:consider changing to long as lattic increases
  int global_y=local_start(y_rank,y_num_procs,y_t_points); 

  int l_st_I=local_start(x_rank,x_num_procs,x_t_points);
  int l_en_I=local_end(x_rank,x_num_procs,x_t_points);
  int l_sz_I=l_en_I-l_st_I;
  int l_st_J=local_start(x_rank,x_num_procs,x_t_points);
  int l_en_J=local_end(x_rank,x_num_procs,x_t_points);
  int l_sz_J=l_en_J-l_st_J;

  //cout<<"got here 3"<<endl;

  for (int i = 0; i < l_sz_I; ++i)
  {
    for (int j = 0; j < l_sz_J; ++j)
    {
      #pragma simd
      for (int n = 0; n < Q; ++n)
      {
	//cout<<"global_x "<<global_x<< " global_y " <<global_y<<endl;
        fEq = wi[n] * exp(-(0.5 / T0) * (((i+global_x) - MIDDLE_X) / SD) * (((i+global_x) - MIDDLE_X) / SD) - (0.5 / T0) * LAMBDA
			  * LAMBDA * (((j+global_y) - MIDDLE_Y) / SD) * (((j+global_y) - MIDDLE_Y) / SD));
      	
        fIn[i][j][n] = fEq;
      	fOut[i][j][n] = fEq;
      }
    }
  }
}

void eq_and_stream(qr_type ***fIn, qr_type ***fOut, qr_type **rho, qr_type **ux, qr_type **uy, const int c[Q][D], const double wi[Q], const int nop[Q], const bool& ftrue)
{
  bool notify = true;
  qr_type u_sqr, c_dot_u, force;
  qr_type x, y, temp, fEq;

  for (int i = 0; i < NX; ++i)
  {
    for (int j = 0; j < NY; ++j)
    {
      // Set to zero before summing
      rho[i][j] = qr_type(0.0);
      ux[i][j] = qr_type(0.0);
      uy[i][j] = qr_type(0.0);

      x = (i - MIDDLE_X) / SD;
      y = LAMBDA * LAMBDA * (j - MIDDLE_Y) / SD;

      for (int n = 0; n < Q; ++n)
	    {
	  	  temp = fIn[i][j][n];
        rho[i][j] += temp;
	      ux[i][j] += c[n][0] * temp;
	      uy[i][j] += c[n][1] * temp;
	    }

      ux[i][j] /= rho[i][j];
      uy[i][j] /= rho[i][j];

      u_sqr = (ux[i][j] * ux[i][j]) + (uy[i][j] * uy[i][j]);

      for (int n = 0; n < Q; ++n)
      {
      	c_dot_u = (c[n][0] * ux[i][j]) + (c[n][1] * uy[i][j]);

      	fEq = wi[n] * rho[i][j] * (1.0 + (1.0 / T0) * c_dot_u + (1.0 / (2.0 * T0 * T0)) * c_dot_u * c_dot_u - (1.0 / (2.0 * T0)) * u_sqr);

      	if (ftrue)
        {
      	  force = -(1.0 / T0) * (((c[n][0] * x) + (c[n][1] * y)) - ((ux[i][j] * x) + (uy[i][j] * y))) * fEq / SD;
      	  if (notify)
          {
      	    cout << "potential is on" << endl;
      	    notify = false;
      	  }
      	}
      	else
        {
      	  force = 0.0;
      	  if (notify)
          {
      	    cout << "potential is off" << endl;
      	    notify = false;
      	  }
      	}

	    fIn[i][j][n] = fIn[i][j][n] * (1.0 - OMEGA) + OMEGA * fEq + force;
      }
    }
  }

  int in, jn;

  for (int i = 0; i < NX; ++i)
  {
    for (int j = 0; j < NY; ++j)
    {
      #pragma simd
      for (int n = 0; n < Q; ++n)
      {
        in = i + int(c[n][0]);
        jn = j + int(c[n][1]);

        if (in > NX - 1 || in < 0)
        {
          in = (in + NX) % NX;                                                                                                                                                                                        
        }
        if (jn > NY-1 || jn < 0)
        {
          jn = (jn + NY) % NY;                                                                                                                                                                                        
        }

        fOut[i][j][nop[n]] = fIn[in][jn][nop[n]];
      }
    }
  }
}

void write_gaussian(qr_type **rho, qr_type **ux, qr_type **uy, const int& ts)
{
  fstream out;
  char fname[255];
  
  sprintf(fname, "data/Xrho_t%i.dat", ts);
  out.open(fname, ios::out);
  for(int i = 0; i < NX; ++i)
  {    
    out << (i - MIDDLE_X) * SIN_V << "\t";
    out << rho[i][MIDDLE_Y] << "\n";
  }

  out.close();

  sprintf(fname, "data/Yrho_t%i.dat", ts);
  out.open(fname, ios::out);
  for(int j = 0; j < NY; ++j)
  {
    out << (j - MIDDLE_Y) * SIN_V << "\t";
    out << rho[MIDDLE_X][j] << "\n";
  }

  out.close();

  sprintf(fname, "data/Xux_t%i.dat", ts);
  out.open(fname, ios::out);
  for(int i = 0; i < NX; ++i)
  {
    out << (i - MIDDLE_X) * SIN_V << "\t";
    out << ux[i][MIDDLE_Y] << "\n";
  }

  out.close();

  sprintf(fname, "data/Yuy_t%i.dat", ts);
  out.open(fname, ios::out);
  for(int j = 0; j < NY; ++j)
  {
    out << (j - MIDDLE_Y) * SIN_V << "\t";
    out << uy[MIDDLE_X][j] << "\n";
  }

  out.close();
}


int local_start(int dir_rank,int dir_num_procs,int dir_t_points)
{
  //return((dir_rank+1)*(dir_t_points/dir_num_procs) + dir_rank*(dir_t_points/dir_num_procs));
  return(dir_rank*(dir_t_points/dir_num_procs)); 
}

int local_end(int dir_rank,int dir_num_procs,int dir_t_points)
{
  //return((dir_rank+1)*(dir_t_points/dir_num_procs) + dir_rank*(dir_t_points/dir_num_procs));
  return((dir_rank+1)*(dir_t_points/dir_num_procs)); 
}



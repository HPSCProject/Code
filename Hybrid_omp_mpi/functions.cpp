#include "functions.h"
#include <omp.h>
using namespace std;
//x corresponds to a row,y corresponds to col
void init_gaussian(qr_type ***fIn, qr_type ***fOut, const double wi[Q],int x_rank,int x_num_procs,int x_t_points,int y_rank,int y_num_procs,int y_t_points,int my_rank,int l_sz_I,int l_sz_J)
{
  qr_type fEq;
  int bl_x=local_start(x_rank,x_num_procs,x_t_points);//FIXME:consider changing to long as lattice sz increases
  int bl_y=local_start(y_rank,y_num_procs,y_t_points); 
  int gl_x,gl_y;
  //#pragma omp parallel for schedule(static) private(fEq)
  for (int i = 0; i <= l_sz_I+1; ++i)
  {
    gl_x = i+bl_x-1;
    if(gl_x < 0) gl_x = x_t_points-1;
    else if(gl_x == x_t_points) gl_x = 0;
    for (int j = 0; j <= l_sz_I+1; ++j)
    {
      gl_y=j+bl_y-1;
      if(gl_y < 0) gl_y = y_t_points-1;//FIXME:costly to check for wrap ard condition at all the procs.identify bdry procs and only apply this to them
      else if(gl_y == y_t_points) gl_y = 0;
      // if(my_rank == 2 ) cout <<"L_sz_I " <<gl_x<<"L_sz_J "<<gl_y<<endl;
      #pragma simd
      for (int n = 0; n < Q; ++n)
      {
	fEq = wi[n] * exp(-(0.5 / T0) * ((gl_x - MIDDLE_X) / SD) * ((gl_x - MIDDLE_X) / SD) - (0.5 / T0) * LAMBDA * LAMBDA * ((gl_y - MIDDLE_Y) / SD) * ((gl_y - MIDDLE_Y) / SD)); 
         fIn[i][j][n] = fEq;
      	 fOut[i][j][n] = fEq;
	 // if(my_rank == 3 ) out <<fIn[i][j][n]<<endl;
      }
    }
  }
  //out.close();
}

void eq_and_stream(qr_type ***fIn, qr_type ***fOut, qr_type **rho, qr_type **ux, qr_type **uy, const int c[Q][D], const double wi[Q], const int nop[Q], const bool& ftrue,int x_rank,int x_num_procs,int x_t_points,int y_rank,int y_num_procs,int y_t_points,int my_rank,int l_sz_I,int l_sz_J)
{
  bool notify = true;
  qr_type u_sqr, c_dot_u, force;
  qr_type x, y, temp, fEq;
    bool introduced = true;

  int bl_x=local_start(x_rank,x_num_procs,x_t_points);//FIXME:consider changing to long as lattice sz increases
  int bl_y=local_start(y_rank,y_num_procs,y_t_points); 
  int gl_x,gl_y;
    cout << "first for bitch!" << endl;
  for (int i = 0; i <= l_sz_I+1; ++i)
  {
    gl_x = i+bl_x-1;
    if(gl_x < 0) gl_x = x_t_points-1;
    else if(gl_x == x_t_points) gl_x = 0;
      cout << "spawn threads bitch!" << endl;
    #pragma omp parallel for private(gl_x,gl_y,x,y,force,c_dot_u,u_sqr,fIn,fEq,rho,ux,uy,temp,introduced,notify)
    for (int j = 0; j <= l_sz_J+1; ++j)
    {
        cout << "spawned the threads bitch!" << endl;
        if(introduced) {
            cout << "Hello I am process #" << my_rank << " my j is " << j << " thread #" << endl;
            //<< omp_get_thread_num() << endl;
            introduced = false;
        }
        
        
      gl_y=j+bl_y-1;
      if(gl_y < 0) gl_y = y_t_points-1;//FIXME:costly to check for wrap ard condition at all the procs.identify bdry procs and only apply this to them
      else if(gl_y == y_t_points) gl_y = 0;
      // if(my_rank == 0 ) cout <<"L_sz_I " <<gl_x<<"L_sz_J "<<gl_y<<endl;
      // Set to zero before summing
      rho[i][j] = qr_type(0.0);
      ux[i][j] = qr_type(0.0);
      uy[i][j] = qr_type(0.0);
        
        cout << "What up bitch i'm here bitch!" << endl;

      x = (gl_x - MIDDLE_X) / SD;
      y = LAMBDA * LAMBDA * (gl_y - MIDDLE_Y) / SD;
      for (int n = 0; n < Q; ++n)
	    {
	      temp = fIn[i][j][n];
	      rho[i][j] += temp;
	      // if(my_rank == 0) out <<fIn[i][j][n]<<endl;
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
        //if(my_rank == 0) out << fEq<<endl;
      	if (ftrue)
        {
      	  force = -(1.0 / T0) * (((c[n][0] * x) + (c[n][1] * y)) - ((ux[i][j] * x) + (uy[i][j] * y))) * fEq / SD;
      	  if (notify)
          {
	    //cout << "potential is on" << endl;
      	    notify = false;
      	  }
      	}
      	else
        {
      	  force = 0.0;
      	  if (notify)
          {
      	    //cout << "potential is off" << endl;
      	    notify = false;
      	  }
      	}
	    fIn[i][j][n] = fIn[i][j][n] * (1.0 - OMEGA) + OMEGA * fEq + force;
            //if(my_rank == 0) out << fIn[i][j][n]<<endl;
      }//end for n
    }//end for j
  }// for i

  fstream out;
  if(my_rank==0)  out.open("data/bpop", ios::out);

  for (int i = 1; i <= l_sz_I; ++i)
  {
    for (int j = 1; j <= l_sz_J; ++j)
    {
       for (int n = 0; n < Q; ++n)
      {
	 if(my_rank == 0) out << fIn[i][j][n]<<endl;
      }
    }
  }
  if (my_rank == 0)  out.close(); 
  
  //dependancy
  int in, jn;

  for (int i = 1; i <= l_sz_I; ++i)
  {
    for (int j = 1; j <= l_sz_J; ++j)
    {
      #pragma simd
      for (int n = 0; n < Q; ++n)
      {
        in = i + int(c[n][0]);
        jn = j + int(c[n][1]);
	/* if (in > NX - 1 || in < 0)
        {in = (in + NX) % NX }
        if (jn > NY-1 || jn < 0)
        {jn = (jn + NY) % NY; }*/
	//if(my_rank==0 ) cout <<"Picking"<<in<<","<<jn<<" to fill "<<nop[n]<<endl;
        fOut[i][j][nop[n]] = fIn[in][jn][nop[n]];
      }
    }
  }

  if(my_rank==0)  out.open("data/apop", ios::out);
  for (int i = 1; i <= l_sz_I; ++i)
  {
    for (int j = 1; j <= l_sz_J; ++j)
    {
       for (int n = 0; n < Q; ++n)
      {
	 if(my_rank == 0) out << fOut[i][j][n]<<endl;
      }
    }
  }
  if (my_rank == 0)  out.close(); 
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
  return(dir_rank*(dir_t_points/dir_num_procs)); 
}

int local_end(int dir_rank,int dir_num_procs,int dir_t_points)
{
  return((dir_rank+1)*(dir_t_points/dir_num_procs)-1); 
}

void cpy_send_cr_pts(qr_type ***fIn,int stop,qr_type *buf,int end,int my_rank)
{
//For Corners procs,move corner points to corner halo points
//Transferring points to the vertical neigbours in tasks cpy_send_cr_pts/cpy_receive_cr_pts \
// These points will then be distributed to horizontall neighbours  in cpy_send/cpy_rx tasks \
// and placed into their respective locations. 
// fIn[0][0] = fIn[1][1];
// fIn[0][stop] = fIn[1][stop-1];
// fIn[stop][0] = fIn[stop-1][1];
// fIn[stop][stop]= fIn[stop-1][stop-1];

 if(end == -1){ //tx data to top
      for (int j=0;j<Q;j++){
	buf[j]=fIn[1][1][j];//1
        buf[Q+j]=fIn[1][stop][j];//2
      }
 }
 if(end == 1){ //tx data to bottom	
      for (int j=0;j<Q;j++){
	buf[j]=fIn[stop][1][j];//4
        buf[Q+j]=fIn[stop][stop][j];//3
      }
 }
}

void cpy_receive_cr_pts(qr_type *temp,qr_type *buf,int end)
//1..2
//4..3
{
 if(end == -1){ //rx data from bottom
      for (int j=0;j<Q;j++){
	temp[(3*Q)+j]=buf[j];//4
        temp[(2*Q)+j]=buf[Q+j];//3
      }
 }

 if(end == 1){ //rx data from top	
      for (int j=0;j<Q;j++){
	temp[j]=buf[j];//1
        temp[Q+j]=buf[Q+j];//2
       }
  }

}

void cpy_send_buf(qr_type ***fIn,int stop,qr_type *buf,int dir,int end,qr_type *temp,int my_rank)
{
int k = 0;
  if(dir==0 && end == -1){ //tx data to left          //1   2                                                                                                                                                                                                    //4   3    
    for(int i=1;i<=stop;i++){
      for (int j=0;j<Q;j++){
	buf[k] = fIn[i][1][j];
        k++;
      }
    }
    for (int j=0;j<Q;j++){
	buf[k]=temp[j];
        buf[Q+k]=temp[(3*Q)+j];
	k++;
	}

  } 

  else if(dir==0 && end == 1){ //tx data to right
    for(int i=1;i<=stop;i++){
      for (int j=0;j<Q;j++){
	buf[k] = fIn[i][stop][j];
	k++;
      }
    }  
    for (int j=0;j<Q;j++){
	buf[k]=temp[Q+j];
        buf[Q+k]=temp[(2*Q)+j];
        if(my_rank == 1)  cout <<"sending "<<k<<"at "<<temp[Q+j]<<endl;
	k++;
	}
 
  }

  else if(dir==1 && end == -1){ //tx data to top
    for(int i=1;i<=stop;i++){
      for (int j=0;j<Q;j++){
	//if(my_rank == 1) cout <<"Tx"<<fIn[1][i][j]<<endl;
	buf[k]=fIn[1][i][j];
	k++;
      }
    }  

  }

  else if(dir==1 && end == 1){ //tx data to bottom
    for(int i=1;i<=stop;i++){
      for (int j=0;j<Q;j++){
	buf[k]=fIn[stop][i][j];
	k++;
      }
    }  
  }

}

void cpy_receive_buf(qr_type ***fIn,int stop,qr_type *buf,int dir,int end,qr_type *temp,int my_rank)
{
 int k = 0;
  //dir = dimension                                                                                                                                                                                            //end = displacement                                                                                                                                                                                       //FIXME:change the start parameter to hardcoded value since it will not change
  if(dir==0 && end == -1){ //rx data from right                                                                                                                                                             
    for(int i=1;i<=stop;i++){
      for (int j=0;j<Q;j++){
        fIn[i][stop+1][j] = buf[k];
	k++;
      }
    }
    for (int j=0;j<Q;j++){
        fIn[0][stop+1][j] = buf[k];
        fIn[stop+1][stop+1][j] = buf[Q+k];
	k++;
	}
  }

  else if(dir==0 && end == 1 ){ //rx data from left          
    for(int i=1;i<=stop;i++){
      for (int j=0;j<Q;j++){
        fIn[i][0][j] = buf[k];
	k++;
      }
    }
    for (int j=0;j<Q;j++){
        if(my_rank == 3)  cout <<"receiving rx "<<k<<"at "<<buf[k]<<endl;
	fIn[0][0][j]=buf[k];
        fIn[stop+1][0][j]=buf[Q+k];
	k++;
	}
  }

  else if(dir==1 && end == -1 ){ //rx data from bottom          
    for(int i=1;i<=stop;i++){
      for (int j=0;j<Q;j++){
	//if(my_rank == 0) cout <<"RX"<<buf[k]<<endl;
        fIn[stop+1][i][j] = buf[k];
	k++;
      }
    }

  }

  else if(dir==1 && end == 1 ){ //rx data from top          
    for(int i=1;i<=stop;i++){
      for (int j=0;j<Q;j++){
        fIn[0][i][j] = buf[k] ;
	k++;
      }
    }
  }
}


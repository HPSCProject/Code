#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
using namespace std;

// usage: ./
int main(int argc,char* argv[])
{
    ifstream infile;
    ofstream outfile;
    //get mtx dimensions from file 1
    int row = atoi(argv[1]);
    int col = atoi(argv[2]);
    int i,j;
    double A[row][col],B[row][col];
    double number;
    int count = 0;
    string line;
    //read file 1
      infile.open(argv[3]);
    if(infile.is_open())
    {
      for(i=0;i<row;i++)
        {
          for(j=0;j<col;j++)
            {
              infile >> number;
              //cout << "The number is:"<<setprecision(5)<<number<<endl;
              A[i][j]=number;
            }
        }
      infile.close();
    }

    //read file 2
      infile.open(argv[4]);
    if(infile.is_open())
    {
      for(i=0;i<row;i++)
        {
          for(j=0;j<col;j++)
            {
              infile >> number;
              //cout << "The number is:"<<setprecision(5)<<number<<endl;
              B[i][j]=number;
            }
        }
      infile.close();
    }

    //find the difference
        for (int i=0;i<row;i++)
        {
          for(int j=0;j<col;j++)
            {
                if((abs(A[i][j] - B[i][j])) > 0.0e-06){
                    cout<<" A and B have different values at location "<<i<<"x"<<j<<"   A value "<<A[i][j]<<"  B value "<<B[i][j]<<endl;
                    count++;}
            }
        }

        if(count)
            cout<<"Comparison failed! There are "<<count <<"differences" <<endl;
        else
            cout<<"Comparison passed!!"<<endl;

     return 0;
}

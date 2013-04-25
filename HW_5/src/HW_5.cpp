#include <iostream>
#include <fstream>
#include <cmath>
#include "nr3.h"
#include "gaussj.h"
#include "ludcmp.h"

using namespace std;

int r = 1; //resistance (ohms) of all resistors

/*
int gaussJordan(int m, int n)
{

  MatDoub_IO R(m,n);
  MatDoub_IO v(m,1);

  //Define LHS
  //first row
  R[0][0]=1;
  R[0][1]=0;
  R[0][2]=4;

  //second row
  R[1][0]=0;
  R[1][1]=2;
  R[1][2]=4;

  //third row
  R[2][0]=1;
  R[2][1]=1;
  R[2][2]=-1;

  //Define RHS
  v[0][0]=1;
  v[1][0]=2;
  v[2][0]=0;

  //perform Gauss-Jordan
  gaussj(R,v);

  cout << "Solution is: \n";
  for (int i=0;i<3;i++)
  {
    cout << i << ":\t" << v[i][0]<<endl;
  }

  cout<<"Finished"<<endl;

  return 0;
}
*/

/*
int luDecomp()
{
	int size=3;

	MatDoub_IO A(size,size);
	VecDoub b(size), x(size);

	//Define LHS

	//first row
	A[0][0]=1;
	A[0][1]=0;
	A[0][2]=4;

	//second row
	A[1][0]=0;
	A[1][1]=2;
	A[1][2]=4;

	//third row
	A[2][0]=1;
	A[2][1]=1;
	A[2][2]=-1;

	//Define RHS
	b[0]=1;
	b[1]=2;
	b[2]=0;

	//Solves for currents via LU Decomposition
	LUdcmp alu(A);

	//solve problem
	alu.solve(b,x);
	cout << "LU Solution is: \n";
	for (int i=0;i<3;i++){
	cout << i << ":\t" << x[i] << endl;
	}
	cout<<"Finished"<<endl;

	return 0;
}

*/

int currentSolver(int m, int n)
{

  /////////////Makes Reference Matrix

    //dimensions of reference matrix with the 0 border
    int M = (m+2);
    int N = (n+2);

    int reference[M][N];
    int refCount;

    //initializes reference matrix;
    for(int i=0; i<M; i++)
    {
		 for(int j=0; j<N; j++)
		 {
		   reference[i][j] = 0;
		 }
    }



    //fills bottom row of reference matrix w/ 1's;
    for(int j=1; j<(N-1); j++)
    {
    	reference[M-2][j] = 1;
    }


    refCount = 2;

    //fills remainings rows of reference matrix
    for(int j=1; j<(N-1); j++)
    {
		for(int i=(M-3); i>0; i--)
		{
		 reference[i][j] = refCount;
		 refCount++;
		}
    }


  	//data output of reference matrix
    	  ofstream of;
  	  of.open("referenceMatrix.csv" );
        //int elements = 0;

  			for(int i=0; i<M ;i++)
  			 {
  				 for(int j=0; j<N ;j++)
  				 {
  					 if((j+1)%N == 0)
  					 {
  						 of << reference[i][j] << "\n";
  					 }

  					 if((j+1)%N != 0)
  					 {
  						 of << reference[i][j] <<",";
  					 }
  				 }
  			 }

  			//make_filename( "body", ".txt" ).c_str()

  			 of << "\n" << "\n";
  			 of.close();




  	/////////////Makes Resistance Matrix

  	//determines dimension of resistance matrix
  	int x = refCount-1;
  	int y = refCount-1;

    MatDoub_IO R(x,y);
    MatDoub_IO v(x,1);


  	  //initializes resistance matrix;
  	  for(int i=0; i<x; i++)
  	  {
  		 for(int j=0; j<y; j++)
  		 {
  		  R[i][j] = 0;
  		 }
  	  }

  	  //defines row-1 of resistance matrix
  	  int resIndex = 0;

  	  for(int j=1; j<(N-1); j++)
  	  {
  		 if(reference[M-3][j] != 0)
  		 {
  			 R[resIndex][ reference[M-3][j]-1 ] = -1;
  		 }
  	  }

  	  R[0][0] = n;




  		resIndex = 1;

  	  //defines remaining rows of resistance matrix
  	  for(int j=1; j<(N-1); j++)
  	  {

  		  for(int i=(M-3); i>0; i--)
  		  {

  			 //look up
  			 if(reference[i-1][j] != 0)
  			 {
  				 R[resIndex][ reference[i-1][j]-1 ] = -1;
  			 }

  			 //look down
  			 if(reference[i+1][j] != 0)
  			 {
  				 R[resIndex][ reference[i+1][j]-1 ] = -1;
  			 }

  			 //look right
  			 if(reference[i][j+1] != 0)
  			 {
  				 R[resIndex][ reference[i][j+1]-1 ] = -1;
  			 }

  			 //look left
  			 if(reference[i][j-1] != 0)
  			 {
  				 R[resIndex][ reference[i][j-1]-1 ] = -1;
  			 }
  			 //Who gets the 4
  			 		 R[resIndex][resIndex] = 4;
  			 resIndex++;

  		  }

  	  }

  	  //data output of Resistance matrix;
  	  		  of.open("RESistanceMatrix.csv" );
  	  	      //int elements = 0;

  	  				for(int i=0; i<x ;i++)
  	  				 {
  	  					 for(int j=0; j<y ;j++)
  	  					 {
  	  			//			 elements++;

  	  						 if((j+1)%y == 0)
  	  						 {
  	  							 of << R[i][j] << "\n";
  	  						 }

  	  						 if((j+1)%y != 0)
  	  						 {
  	  							 of << R[i][j] <<",";
  	  						 }
  	  					 }
  	  				 }

  	  				//make_filename( "body", ".txt" ).c_str()

  	  				 of << "\n" << "\n";
  	  				 of.close();



  	 //Define Voltage Vector (RHS)
	 for(int i=0; i<x; i++)
	 {
		 v[i][0]=0;
	 }

 	 v[0][0]=1;

  //Solves for current via Gauss-Jordan
  gaussj(R,v);

  cout << "Gauss Jordan solution is: \n";
  for (int i=0;i<x;i++)
  {
    cout << i << ":\t" << v[i][0]<<endl;
  }

  cout<<"Finished Gauss Jordan"<<endl<<endl;



  return 0;
}


int main()
{
	currentSolver(3,2);
	//resistance(3,3);
	return 0;
}
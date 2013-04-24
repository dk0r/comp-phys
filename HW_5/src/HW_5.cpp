#include <iostream>
#include <fstream>
#include <cmath>
#include "nr3.h"
#include "gaussj.h"

using namespace std;

int r = 1; //resistance (ohms) of all resistors


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


int resistance(int m, int n)
{


 //Makes Reference Matrix
 ////////////////////////////////////////////////////////////////////////////////////

  //dimension of reference matrix with 0 buffer
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
	  //cout << "reference[" << i << "][" << j << "]=" << reference[i][j] << endl;
	 }
  }



  //fills bottom loop of reference matrix w/ 1's;
  for(int j=1; j<(N-1); j++)
  {
	reference[M-2][j] = 1;
  }


  refCount = 2;

  for(int j=1; j<(N-1); j++)
  {
	for(int i=(M-3); i>0; i--)
	{
	 reference[i][j] = refCount;
	 refCount++;
	}
  }


	//data output
  	  ofstream of;
	  of.open("referenceMatrix.csv" );
      //int elements = 0;

			for(int i=0; i<M ;i++)
			 {
				 for(int j=0; j<N ;j++)
				 {
		//			 elements++;

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

/////////////////////////////////////////////////////////////////
	//Makes Resistance Matrix
	int x = refCount-1;
	int y = refCount-1;



	int resistance[x][y];
	  //initializes resistance matrix;
	  for(int i=0; i<x; i++)
	  {
		 for(int j=0; j<y; j++)
		 {
		  resistance[i][j] = 0;
		 }
	  }

	  //defines line-1 of resistance matrix

	  resistance[0][0] = n;

	  int resIndex = 0;

	  for(int j=1; j<(N-1); j++)
	  {
		 if(reference[M-3][j] != 0)
		 {
			 resistance[resIndex][ reference[M-3][j]-1 ] = -1;
		 }
	  }


		resIndex = 1;

	  for(int j=1; j<(N-1); j++)
	  {

		  for(int i=(M-3); i>0; i--)
		  {

			 //look up
			 if(reference[i-1][j] != 0)
			 {
				 resistance[resIndex][ reference[i-1][j]-1 ] = -1;
			 }

			 //look down
			 if(reference[i+1][j] != 0)
			 {
				 resistance[resIndex][ reference[i+1][j]-1 ] = -1;
			 }

			 //look right
			 if(reference[i][j+1] != 0)
			 {
				 resistance[resIndex][ reference[i][j+1]-1 ] = -1;
			 }

			 //look left
			 if(reference[i][j-1] != 0)
			 {
				 resistance[resIndex][ reference[i][j-1]-1 ] = -1;
			 }
			 //Who gets the 4
			 		  resistance[resIndex][resIndex] = 4;
			 resIndex++;

		  }

	  }


	  //data output
	  		  of.open("RESistanceMatrix.csv" );
	  	      //int elements = 0;

	  				for(int i=0; i<x ;i++)
	  				 {
	  					 for(int j=0; j<y ;j++)
	  					 {
	  			//			 elements++;

	  						 if((j+1)%y == 0)
	  						 {
	  							 of << resistance[i][j] << "\n";
	  						 }

	  						 if((j+1)%y != 0)
	  						 {
	  							 of << resistance[i][j] <<",";
	  						 }
	  					 }
	  				 }

	  				//make_filename( "body", ".txt" ).c_str()

	  				 of << "\n" << "\n";
	  				 of.close();


	  		gaussJordan(resistance);




  return 0;
}


int main()
{
	resistance(3,3);
	return 0;
}

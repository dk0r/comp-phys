#include <iostream>
#include <fstream>
#include <cmath>
#include <nr3.h>
#include <rk4.h>

using namespace std;


int test()
{
	int x = 100;
	int y = 100;
	int n = 100;
	int periodic = 0;
	double func[x][y];
	ofstream of;
	of.open("func.txt");
	of <<"test"<<endl;
	of.close();

	bool donotcompute[x][y];

	//initializes grid
	 for(int ix=0; ix<n; ix++)
	 {
		 for (int iy=0; iy<n; iy++)
		 {
			 donotcompute[ix][iy] = false;
			 func[ix][iy] = 0.0;
		 }
	 }

	 //defines: 1. upper and lower horizontal boundary conditions
	 //         2.
	 for(int iy=0; iy<n; iy++)
	 {
		 if(periodic==0)
		 {
			 donotcompute[0][iy] = true;
			 donotcompute[n-1][iy] = true;
		 }

		 func[0][iy] = 10;
		 func[n-1][iy] = 10;
	 }

	 //defines left and right boundary conditions
	 for(int ix=0; ix<n; ix++)
	 {
		 donotcompute[ix][0] = true;
		 donotcompute[ix][n-1] = true;
		 func[ix][0] = 10;
		 func[ix][n-1] = 10;
	 }




	 //jacobian relaxation
	 for(int ix=0; ix<n; x++)
	 {
		 for(int iy=0; iy<n; iy++)
		 {
			 if(donotcompute[ix][iy]==false)
			{
			 func[ix][iy] = 0.25*(func[ix+1][iy]+func[ix-1][iy]+func[ix][iy+1]+func[ix][iy-1]);
			}
		 }
	}

	 return 0;
}

int main()
{

}

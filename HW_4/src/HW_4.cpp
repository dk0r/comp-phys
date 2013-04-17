#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <nr3.h>
#include <rk4.h>
#include <time.h>

clock_t init, final;

using namespace std;
int fileNamerCounter = 0;
string make_filename( const string& basename, const string& ext )
  {
  ostringstream result;
  result << basename << fileNamerCounter << ext;
  fileNamerCounter++;
  return result.str();
  }


int jacobi()
{
	int x = 7;
	int y = 7;
	int n = 7;
	int periodic = 0;
	double func[x][y];
	double old[x][y];
	double voltage = 10;


	bool donotcompute[x][y];

	//initializes grid
	 for(int ix=0; ix<n; ix++)
	 {
		 for (int iy=0; iy<n; iy++)
		 {
			 donotcompute[ix][iy] = false;
			 func[ix][iy] = 0.0;
			 old[ix][iy] = 0.0;
		 }
	 }



	 //defines: 1. upper and lower horizontal boundary conditions
	 //         2. remaining non-boundaries
	 for(int iy=0; iy<n; iy++)
	 {
		 if(periodic==0)
		 {
			 donotcompute[0][iy] = true;
			 donotcompute[n-1][iy] = true;
		 }

		 old[0][iy] = voltage;
		 old[n-1][iy] = voltage;

		 func[0][iy] = voltage;
		 func[n-1][iy] = voltage;

	 }


	 //defines left and right boundary conditions
	 for(int ix=0; ix<n; ix++)
	 {
		 donotcompute[ix][0] = true;
		 donotcompute[ix][n-1] = true;

		 old[ix][0] = voltage;
		 old[ix][n-1] = voltage;

		 func[ix][0] = voltage;
		 func[ix][n-1] = voltage;
	 }


	 int elements = 0;
	 int jacobiCounter = 0;
	 double stop = 0;
	 ofstream of;


	while(voltage - stop > 0.0000000000001)
	{

				 for(int iy=0; iy<n ;iy++)
				 {
					 for(int ix=0; ix<n ;ix++)
					 { //cout << "Entering Bottom Loop" << endl;
							 if(donotcompute[ix][iy]==false)
							 {
								 func[ix][iy] = 0.25*(old[ix+1][iy]+old[ix-1][iy]+old[ix][iy+1]+old[ix][iy-1]);
							 }
					 }

				 }

				 //writes complete func[][] values to old[][] for next iteration
				 for(int iy=0; iy<n ;iy++)
				 {
					 for(int ix=0; ix<n ;ix++)
					 {
							 if(donotcompute[ix][iy]==false)
							 {
								 old[ix][iy] = func[ix][iy];
							 }
					 }

				 }


				//data output
				of.open(make_filename( "/home/dk0r/csv/jacobi", ".csv" ).c_str());

				for(int i=0; i<n ;i++)
				 {
						 for(int j=0; j<n ;j++)
						 {
								 elements++;

								 if(elements%n == 0)
								 {
									 of << old[i][j] << "\n";
								 }

								 if(elements%n != 0)
								 {
									 of << old[i][j] <<",";
								 }
						 }
				 }

				//make_filename( "body", ".txt" ).c_str()

				 of << "\n" << "\n";
				 //of.close();

				 stop = func[(int)(n/2)][(int)(n/2)];
				 jacobiCounter++;
				 cout << jacobiCounter << ")  " << "Voltage - " << setprecision(15) << func[(int)(n/2)][(int)(n/2)]  << "  =  " << voltage - func[(int)(n/2)][(int)(n/2)] << endl;

	}
	cout << "All DooOOOOooOOooOOooOne" << endl;
	of.close();
	 return 0;
}


int gaussSeidel()
{
	int x = 7;
	int y = 7;
	int n = 7;
	int periodic = 0;
	double func[x][y];
	double old[x][y];
	double voltage = 10;


	bool donotcompute[x][y];

	//initializes grid
	 for(int ix=0; ix<n; ix++)
	 {
		 for (int iy=0; iy<n; iy++)
		 {
			 donotcompute[ix][iy] = false;
			 func[ix][iy] = 0.0;
			 old[ix][iy] = 0.0;
		 }
	 }



	 //defines: 1. upper and lower horizontal boundary conditions
	 //         2. remaining non-boundaries
	 for(int iy=0; iy<n; iy++)
	 {
		 if(periodic==0)
		 {
			 donotcompute[0][iy] = true;
			 donotcompute[n-1][iy] = true;
		 }

		 old[0][iy] = voltage;
		 old[n-1][iy] = voltage;

		 func[0][iy] = voltage;
		 func[n-1][iy] = voltage;

	 }


	 //defines left and right boundary conditions
	 for(int ix=0; ix<n; ix++)
	 {
		 donotcompute[ix][0] = true;
		 donotcompute[ix][n-1] = true;

		 old[ix][0] = voltage;
		 old[ix][n-1] = voltage;

		 func[ix][0] = voltage;
		 func[ix][n-1] = voltage;
	 }


	 int elements = 0;
	 int gaussSeidelCounter = 0;
	 double stop = 0;
	 ofstream of;


	while(voltage - stop > 0.0000000000001)
	{

				 for(int iy=0; iy<n ;iy++)
				 {
					 for(int ix=0; ix<n ;ix++)
					 { //cout << "Entering Bottom Loop" << endl;
							 if(donotcompute[ix][iy]==false)
							 {
								 func[ix][iy] = 0.25*(func[ix+1][iy]+func[ix-1][iy]+func[ix][iy+1]+func[ix][iy-1]);
							 }
					 }

				 }


				//data output
				of.open(make_filename( "/home/dk0r/csv/gaussSeidel", ".csv" ).c_str());

				for(int i=0; i<n ;i++)
				 {
						 for(int j=0; j<n ;j++)
						 {
								 elements++;

								 if(elements%n == 0)
								 {
									 of << old[i][j] << "\n";
								 }

								 if(elements%n != 0)
								 {
									 of << old[i][j] <<",";
								 }
						 }
				 }

				//make_filename( "body", ".txt" ).c_str()

				 of << "\n" << "\n";
				 //of.close();

				 stop = func[(int)(n/2)][(int)(n/2)];
				 gaussSeidelCounter++;
				 cout << gaussSeidelCounter << ")  " << "Voltage - " << setprecision(15) << func[(int)(n/2)][(int)(n/2)]  << "  =  " << voltage - func[(int)(n/2)][(int)(n/2)] << endl;

	}
	cout << "All DooOOOOooOOooOOooOne" << endl;
	of.close();
	 return 0;
}

int ouRelaxation(double alpha)
{
	int x = 7;
	int y = 7;
	int n = 7;
	int periodic = 0;
	double func[x][y];
	double old[x][y];
	double voltage = 10;


	bool donotcompute[x][y];

	//initializes grid
	 for(int ix=0; ix<n; ix++)
	 {
		 for (int iy=0; iy<n; iy++)
		 {
			 donotcompute[ix][iy] = false;
			 func[ix][iy] = 0.0;
			 old[ix][iy] = 0.0;
		 }
	 }



	 //defines: 1. upper and lower horizontal boundary conditions
	 //         2. remaining non-boundaries
	 for(int iy=0; iy<n; iy++)
	 {
		 if(periodic==0)
		 {
			 donotcompute[0][iy] = true;
			 donotcompute[n-1][iy] = true;
		 }

		 old[0][iy] = voltage;
		 old[n-1][iy] = voltage;

		 func[0][iy] = voltage;
		 func[n-1][iy] = voltage;

	 }


	 //defines left and right boundary conditions
	 for(int ix=0; ix<n; ix++)
	 {
		 donotcompute[ix][0] = true;
		 donotcompute[ix][n-1] = true;

		 old[ix][0] = voltage;
		 old[ix][n-1] = voltage;

		 func[ix][0] = voltage;
		 func[ix][n-1] = voltage;
	 }


	 int elements = 0;
	 int gaussSeidelCounter = 0;
	 double stop = 0;
	 ofstream of;


	while(voltage - stop > 0.0000000000001)
	{

				 for(int iy=0; iy<n ;iy++)
				 {
					 for(int ix=0; ix<n ;ix++)
					 { //cout << "Entering Bottom Loop" << endl;
							 if(donotcompute[ix][iy]==false)
							 {
								 func[ix][iy] = 0.25*(func[ix+1][iy]+func[ix-1][iy]+func[ix][iy+1]+func[ix][iy-1]);
							 }
					 }

				 }


				//data output
				of.open(make_filename( "/home/dk0r/csv/gaussSeidel", ".csv" ).c_str());

				for(int i=0; i<n ;i++)
				 {
						 for(int j=0; j<n ;j++)
						 {
								 elements++;

								 if(elements%n == 0)
								 {
									 of << old[i][j] << "\n";
								 }

								 if(elements%n != 0)
								 {
									 of << old[i][j] <<",";
								 }
						 }
				 }

				//make_filename( "body", ".txt" ).c_str()

				 of << "\n" << "\n";
				 //of.close();

				 stop = func[(int)(n/2)][(int)(n/2)];
				 gaussSeidelCounter++;
				 cout << gaussSeidelCounter << ")  " << "Voltage - " << setprecision(15) << func[(int)(n/2)][(int)(n/2)]  << "  =  " << voltage - func[(int)(n/2)][(int)(n/2)] << endl;

	}
	cout << "All DooOOOOooOOooOOooOne" << endl;
	of.close();
	 return 0;
}



int main()
{

	//jacobi();


	//gaussSeidel();


	return 0;
}

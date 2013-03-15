
#include <iostream>
#include <fstream>
#include <cmath>
#include <nr3.h>
#include <rk4.h>

using namespace std;


////Physical Constants////

//Gravitational Constant (SI-units:  m^3 / kg*s^2)
double G  = 1; //6.67398e-11;

//Solar Mass (SI-units: kg)
double M  = 1; //1.9891e30;

//Heliocentric Gravitational Constant (Astro-units: au^3/days^2)
double GM = 2.95912208e-4;





void derivs(const Doub x, VecDoub_I & y, VecDoub_O & dydx)
{
	dydx[0] = y[1];
	dydx[1] =  ((-G*M) *  y[0]) / (pow( (y[0]*y[0] + y[2]*y[2]) , (1.5)));
	dydx[2] = y[3];
	dydx[3] = ((-G*M) *  y[2]) / pow( y[0]*y[0] + y[2]*y[2] , 1.5);
}




int main ()
{
	VecDoub y(4),dydx(4); //
	Doub x, xmin, kmax=5000, h=0.001;

	VecDoub yout(4);
	int k;

	xmin = 0;

	//Initial Conditions:  If using GM [NOT G*M] in derivs(), units are kg, au and days
	y[0] = 0.5;  //units: au_x
	y[1] = 0;    //units: au_x/day
	y[2] = 0;    //units au_y
	y[3] = 1.63; //units au_u/day

	derivs(xmin, y, dydx);

	//file output streams
	ofstream ofALL, ofX, ofY, ofXY;
	ofALL.open("all-data.csv");
	ofX.open("x-components.csv");
	ofY.open("y-components.csv");
	ofXY.open("x-y-position.csv");


	for(k=0; k < kmax; k++)
	{
			x=xmin+k*h;

			rk4(y, dydx,  x, h, yout, derivs);

			//display output
			cout << "k = " << k
				 << "    x = " << yout[0] << "    x'= " << yout[1]
				 << "    y = " << yout[2] << "    y'= " << yout[3] << endl;

			//file output streams
			ofALL << yout[0] << "," << yout[1] << "," << yout[2] << "," << yout[3] << endl;
			ofX << yout[0] << "," << yout[1] << endl;
			ofY << yout[2] << "," << yout[3] << endl;
			ofXY << yout[0] << "," << yout[2] << endl;

			y = yout;

			derivs(x,y,dydx);
	}

	//closes file output streams
	ofALL.close();
	ofX.close();
	ofY.close();
	ofXY.close();
}

/*

//Last Working w/ NumericalMethodsGuy

void derivs(const Doub x, VecDoub_I & y, VecDoub_O & dydx)
{
	dydx[0]= y[1];
	dydx[1]= ( 11*exp((-1)*x) - 3*y[1] - 5*y[0] ) / 2;

}




int main (int argc, char * const argv[])
{
	VecDoub y(2),dydx(2);
	Doub x, xmin, xmax, kmax=1000, h=0.25;

	VecDoub yout(2);
	int k;

	xmin=0;
	//xmax=2; //appears to do nothing

	y[0]=7;
	y[1]=13;

	derivs(xmin,y,dydx);

	ofi.open("test.csv");
	ofi << "f(x+h),f'(x+h)" << endl;

	for(k=0;k<kmax;k++)
	{
		x=xmin+k*h;
		rk4(y, dydx,  x, h, yout, derivs);
		cout << "k = " << k << "    f(x+h) = " << yout[0] << "    f'(x+h) = " << yout[1] << endl;
		ofi << yout[0] << "," << yout[1] << endl;
		y[0]=yout[0];
		y[1]=yout[1];
		derivs(x,y,dydx);
	}

}

*/


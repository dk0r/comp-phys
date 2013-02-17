#include <iostream>
#include <math.h>
#include <iomanip>
#include <fstream>
using namespace std;



double e = 2.71828182845904523536028747135266249775724709369995;
double pi = atan(1)*4;


double epsilon()
	{
		double x, epsilon;
		x=1;

		for(;1;)
		{
			x/=2;

			if(x+1==1)
				break;
		}

		epsilon=x*2;

		return epsilon;
	}


double trapezoid(int i, double a, double b)
	{
		//Known integral of x^2 between 1 and 5
		double known = 124/3;

		ofstream fo;
		fo.open("trapezoid.csv");
		fo << "n,epsilon" << endl;

		for(int n=1; n<=i; n++)
		{
			fo << n << ",";

			//Determine interval width
			double deltaX = (b-a)/n;
																		///cout << "deltaX = " << deltaX << endl;
			//Calculates areas only at limits of integration
			double area = deltaX/2 * (pow(a,2));
				   area += deltaX/2 * (pow(b, 2));				  //Stores value of initial area @ limits of integration
																		///double limitArea = area;
																		///cout << "Area @ limits = " << limitArea << endl;
			//Calculates remaining areas between limits of integration
			for(int j=1; j<n; j++)
			{
				double xSub_j = a + j*(deltaX);
				area += ( deltaX * (pow(xSub_j,2)) );
			}
																		///cout << "Area inbetween = " << (area - limitArea ) << endl;
			// cout << "n = " << n << "   area = " << setprecision(17) << area << endl;

			fo << setprecision(17) << fabs( (area-known)/known ) << endl;

			area = 0;
		}

		fo.close();
		return 0;
	}

double simpson(int n, double a, double b)
	{   //Determine interval width
		double deltaX = (b-a)/n;
																	///cout << "deltaX = " << deltaX << endl;
		//Calculates areas only at limits of integration
		double area = deltaX/3 * pow(a,2);
			   area += deltaX/3 * pow(b,2);				  //Stores value of initial area @ limits of integration
																	///double limitArea = area;
																	///cout << "Area @ limits = " << limitArea << endl;
		//Calculates remaining areas between limits of integration
		for(int i=1; i<n; i++)
		{
			double xSub_i = a + i*(deltaX);

			if(i%2 == 0)
			{
				area += ( deltaX/3 * 2 * pow(xSub_i,2) );
			}

			else
			{
				area += ( deltaX/3 * 4 * pow(xSub_i,2) );
			}
		}
																	///cout << "Area inbetween = " << (area - limitArea ) << endl;
		return area;
	}

double gauss3pt(double a, double b)
	{
		//defines 3pt gauss-quadrature weighting(c)/function(x) factors
		double c[3] = {0.555555556, 0.888888889, 0.555555556};
		double x[3] = {-0.774596669, 0.000000000, 0.774596669};

		//Determine interval width
		double deltaX1 = (b-a)/2;
		double deltaX2 = (b+a)/2;
																	///cout << "deltaX = " << deltaX << endl;
		//Holds area from gauss-quad
		double area = 0;

		for(int i=0; i<3; i++)
		{
			area += c[i] * deltaX1 * pow( (deltaX1*x[i]+deltaX2) , 2 );
		}

		return area;
	}



double forwardDiffx2()
	{
		double area = 0;
		double error = 0;
		double deltaX =0.00000000001;
		double x = 3;

		ofstream of;
		of.open("forwardDiffx2.csv");
		//of << "h,epsilon" << endl;

		for(int i=1; i<1000000; i++)
		{

				area = ( pow((x+deltaX),2) - pow(x,2) ) / deltaX;

				of << setprecision(7) << deltaX << ",";
			cout << "i=" << i << ")  h = " << setprecision(7) << deltaX;
				deltaX -= 0.00000000000000005;

				error = fabs( (area-6)/6 );

				of << setprecision(17) << error << endl;

			cout << "    area = " << setprecision(17) << area;
			cout << "    RelError = " << setprecision(8) << error*100 << "%" << endl;
			if(error>=1)
				{break;}
			error = 0;
			area = 0;

		}

		of.close();
		return 0;
	}

double centralDiffx2()
{
		double area = 0;
		double error = 0;
		double deltaX =0.00000000001;
		double x = 3;

		ofstream of;
		of.open("centralDiffx2.csv");
		//of << "h,epsilon" << endl;

		for(int i=1; i<1000000; i++)
		{

			area = ( (pow((x+(deltaX/2)),2))-(pow(x-(deltaX/2),2)) ) / deltaX;

				of << setprecision(7) << deltaX << ",";
			cout << "i=" << i << ")  h = " << setprecision(7) << deltaX;
				deltaX -= 0.00000000000000005;

				error = fabs( (area-6)/6 );

				of << setprecision(17) << error << endl;

			cout << "    area = " << setprecision(17) << area;
			cout << "    RelError = " << setprecision(8) << error*100 << "%" << endl;
			if(error>=1)
				{break;}
			error = 0;
			area = 0;

		}

		of.close();
		return 0;
	}

double fivePointx2()
{
	double area = 0;
	double error = 0;
	double deltaX =0.00000000001;
	double x = 3;

	ofstream of;
	of.open("fivePointx2.csv");
	//of << "h,epsilon" << endl;

	for(int i=1; i<1000000; i++)
	{

		area = (  (-1)*(pow((x+(2*deltaX)),2)) + (8)*(pow((x+deltaX),2)) - (8)*(pow((x-deltaX),2)) + pow((x-(2*deltaX)),2)  )/(12*deltaX);


			of << setprecision(7) << deltaX << ",";
		cout << "i=" << i << ")  h = " << setprecision(7) << deltaX;
			deltaX -= 0.00000000000000005;

			error = fabs( (area-6)/6 );

			of << setprecision(17) << error << endl;

		cout << "    area = " << setprecision(17) << area;
		cout << "    RelError = " << setprecision(8) << error*100 << "%" << endl;
		if(error>=1)
			{break;}
		error = 0;
		area = 0;

	}

	of.close();
	return 0;
}


double forwardDiffx3()
	{
		double area = 0;
		double error = 0;
		double deltaX =0.00000000001;
		double x = 3;

		ofstream of;
		of.open("forwardDiffx3.csv");
		//of << "h,epsilon" << endl;

		for(int i=1; i<1000000; i++)
		{

				area = ( pow((x+deltaX),3) - pow(x,3) ) / deltaX;

				of << setprecision(7) << deltaX << ",";
			cout << "i=" << i << ")  h = " << setprecision(7) << deltaX;
				deltaX -= 0.00000000000000005;

				error = fabs( (area-27)/27 );

				of << setprecision(17) << error << endl;

			cout << "    area = " << setprecision(17) << area;
			cout << "    RelError = " << setprecision(8) << error*100 << "%" << endl;
			if(error>=1)
				{break;}
			error = 0;
			area = 0;

		}

		of.close();
		return 0;
	}

double centralDiffx3()
{
		double area = 0;
		double error = 0;
		double deltaX =0.00000000001;
		double x = 3;

		ofstream of;
		of.open("centralDiffx3.csv");
		//of << "h,epsilon" << endl;

		for(int i=1; i<1000000; i++)
		{

			area = ( (pow((x+(deltaX/2)),3))-(pow(x-(deltaX/2),3)) ) / deltaX;

				of << setprecision(7) << deltaX << ",";
			cout << "i=" << i << ")  h = " << setprecision(7) << deltaX;
				deltaX -= 0.00000000000000005;

				error = fabs( (area-27)/27 );

				of << setprecision(17) << error << endl;

			cout << "    area = " << setprecision(17) << area;
			cout << "    RelError = " << setprecision(8) << error*100 << "%" << endl;
			if(error>=1)
				{break;}
			error = 0;
			area = 0;

		}

		of.close();
		return 0;
	}

double fivePointx3()
{
	double area = 0;
	double error = 0;
	double deltaX =0.00000000001;
	double x = 3;

	ofstream of;
	of.open("fivePointx3.csv");
	//of << "h,epsilon" << endl;

	for(int i=1; i<1000000; i++)
	{

		area = (  (-1)*(pow((x+(2*deltaX)),3)) + (8)*(pow((x+deltaX),3)) - (8)*(pow((x-deltaX),3)) + pow((x-(2*deltaX)),3)  )/(12*deltaX);


			of << setprecision(7) << deltaX << ",";
		cout << "i=" << i << ")  h = " << setprecision(7) << deltaX;
			deltaX -= 0.00000000000000005;

			error = fabs( (area-27)/27 );

			of << setprecision(17) << error << endl;

		cout << "    area = " << setprecision(17) << area;
		cout << "    RelError = " << setprecision(8) << error*100 << "%" << endl;
		if(error>=1)
			{break;}
		error = 0;
		area = 0;

	}

	of.close();
	return 0;
}



double forwardDiffE()
{
	double area = 0;
	double error = 0;
	double deltaX =0.00000000001;
	double x = 3;

	ofstream of;
	of.open("forwardDiffE.csv");
	//of << "h,epsilon" << endl;

	for(int i=1; i<1000000; i++)
	{

		area = ( pow(e,-(x+deltaX)) - pow(e,-x) ) / deltaX;

			of << setprecision(7) << deltaX << ",";
		cout << "i=" << i << ")  h = " << setprecision(7) << deltaX;
			deltaX -= 0.00000000000000005;

			error = fabs( (area-(pow(-e,-3)))/pow(-e,-3) );

			of << setprecision(17) << error << endl;

		cout << "    area = " << setprecision(17) << area;
		cout << "    RelError = " << setprecision(8) << error*100 << "%" << endl;
		if(error>=1)
			{break;}
		error = 0;
		area = 0;

	}

	of.close();
	return 0;
}

double centralDiffE()
{
		double area = 0;
		double error = 0;
		double deltaX =0.00000000001;
		double x = 3;

		ofstream of;
		of.open("centralDiffE.csv");
		//of << "h,epsilon" << endl;

		for(int i=1; i<1000000; i++)
		{

			area = ( (pow(e,(-1)*(x+(deltaX/2))))-(pow(e,(-1)*(x-(deltaX/2)))) ) / deltaX;

				of << setprecision(7) << deltaX << ",";
			cout << "i=" << i << ")  h = " << setprecision(7) << deltaX;
				deltaX -= 0.00000000000000005;

				error = fabs( (area-(pow(-e,-3)))/pow(-e,-3) );

				of << setprecision(17) << error << endl;

			cout << "    area = " << setprecision(17) << area;
			cout << "    RelError = " << setprecision(8) << error*100 << "%" << endl;
			if(error>=1)
				{break;}
			error = 0;
			area = 0;

		}

		of.close();
		return 0;
	}

double fivePointE()
{
	double area = 0;
	double error = 0;
	double deltaX =0.00000000001;
	double x = 3;

	ofstream of;
	of.open("fivePointE.csv");
	//of << "h,epsilon" << endl;

	for(int i=1; i<1000000; i++)
	{

		area = (  (-1)*(pow(e,(-1)*((x+(2*deltaX))))) + (8)*(pow(e,(-1)*(x+deltaX))) - (8)*(pow(e,(-1)*(x-deltaX))) + pow(e,(-1)*((x-(2*deltaX))))  )/(12*deltaX);


			of << setprecision(7) << deltaX << ",";
		cout << "i=" << i << ")  h = " << setprecision(7) << deltaX;
			deltaX -= 0.00000000000000005;

			error = fabs( (area-(pow(-e,-3)))/pow(-e,-3) );

			of << setprecision(17) << error << endl;

		cout << "    area = " << setprecision(17) << area;
		cout << "    RelError = " << setprecision(8) << error*100 << "%" << endl;
		if(error>=1)
			{break;}
		error = 0;
		area = 0;

	}

	of.close();
	return 0;
}


double forwardDiffsinln()
{
		double area = 0;
		double error = 0;
		double deltaX =0.00000000001;
		double x = pi;

		ofstream of;
		of.open("forwardDiffsinln.csv");
		//of << "h,epsilon" << endl;

		for(int i=1; i<1000000; i++)
		{

			area = ( (pow(sin(x+deltaX),2)+log(x+deltaX)) - (pow(sin(x),2)+log(x) ) ) / deltaX;

				of << setprecision(7) << deltaX << ",";
			cout << "i=" << i << ")  h = " << setprecision(7) << deltaX;
				deltaX -= 0.00000000000000005;

				error = fabs( (area-(1/pi))/(1/pi) );

				of << setprecision(17) << error << endl;

			cout << "    area = " << setprecision(17) << area;
			cout << "    RelError = " << setprecision(8) << error*100 << "%" << endl;
			if(error>=1)
				{break;}
			error = 0;
			area = 0;

		}

		of.close();
		return 0;
	}

double centralDiffsinln()
{
		double area = 0;
		double error = 0;
		double deltaX =0.00000000001;
		double x = pi;

		ofstream of;
		of.open("centralDiffsinln.csv");
		//of << "h,epsilon" << endl;

		for(int i=1; i<1000000; i++)
		{

			area = ( ((pow(sin(x+(deltaX/2)),2)+log(x+(deltaX/2))) - (pow(sin(x-(deltaX/2)),2)+log(x-(deltaX/2))) ) / deltaX);


				of << setprecision(7) << deltaX << ",";
			cout << "i=" << i << ")  h = " << setprecision(7) << deltaX;
				deltaX -= 0.00000000000000005;

				error = fabs( (area-(1/pi))/(1/pi) );

				of << setprecision(17) << error << endl;

			cout << "    area = " << setprecision(17) << area;
			cout << "    RelError = " << setprecision(8) << error*100 << "%" << endl;
			if(error>=1)
				{break;}
			error = 0;
			area = 0;

		}

		of.close();
		return 0;
	}

double fivePointsinln()
{
	double area = 0;
	double error = 0;
	double deltaX =0.00000000001;
	double x = pi;

	ofstream of;
	of.open("fivePointsinln.csv");
	//of << "h,epsilon" << endl;

	for(int i=1; i<1000000; i++)
	{
		area = ((-1)*(pow(sin((x+(2*deltaX))),2) + log((x+(2*deltaX)))) + ((8)*(pow(sin((x+deltaX)),2) + log((x+deltaX)) )) - ((8)*(pow(sin((x-deltaX)),2) + log((x-deltaX)))) + (pow(sin((x-(2*deltaX))),2) + log((x-(2*deltaX))))) / (12*deltaX);


			of << setprecision(7) << deltaX << ",";
		cout << "i=" << i << ")  h = " << setprecision(7) << deltaX;
			deltaX -= 0.00000000000000005;

			error = fabs( (area-(1/pi))/(1/pi) );

			of << setprecision(17) << error << endl;

		cout << "    area = " << setprecision(17) << area;
		cout << "    RelError = " << setprecision(8) << error*100 << "%" << endl;
		if(error>=1)
			{break;}
		error = 0;
		area = 0;

	}

	of.close();
	return 0;
}



double newRaph(double x_0)
{
	double x_i = x_0;
	double x_iPlus1 = 0;
	double error = 0;

	for(int i=0; i<10000000; i++)
	{
		x_iPlus1 = (  x_i - (pow(x_i,3)-4) / (3*(pow(x_i,2)))  );
		cout << "x_i = " << x_i;
		cout << "  x_iPlus1 = " << x_iPlus1;

		error = fabs(    ( (x_iPlus1)-(x_i) )/(x_iPlus1)    );
		cout << "    error = " << error*100 << "%" << endl;
		x_i = x_iPlus1;

		if(error < 0.000000001)
			break;

	}


	return error*100;

}

double bisection(double l, double u)

	{
		double m = 0;

					//Determines whether user provided bounds contains a root.
		if((pow(l,2)-1)*(pow(u,2)-1) >= 0)
			{
				cout << "Invalid Bounds: Bounds must contain root" << endl;
				return 0;
			}

					//Finds root within user provided bounds
		for(int i=0; i<1000; i++)
		{

					//Defines midpoint
			m = (l+u)/2;

					//Found Root
			if((pow(l,2)-1)*(pow(m,2)-1) == 0)
			{
				cout << "Found Zero @ " << m << endl;
				l = m;

				break;
			}

					//Excludes upper bound from root search
			if((pow(l,2)-1)*(pow(m,2)-1) < 0)
			{
				u = m;
			}

					//Excludes lower bound from root search
			if((pow(l,2)-1)*(pow(m,2)-1) > 0)
			{
				l = m;
			}

		}

	return 0;

	}


int main()
	{
		//int n;		//# of trapezoid(s)
		//int a,b;	//lower/upper limits of integration
		//int x, h;
		//cout << "Enter:          Trapezoid quantity" << endl;
		//cout << "        Lower limit of integration" << endl;
		//cout << "       Upper limit of integration" << endl;
		//cin >> n >> a >> b;

		//cout << "Enter:   point to eval (x)" << endl << "      stepsize(h)" << endl;
		//cin >> x >> h;
		//forwardDiff();

		//forwardDiffsinln();

		cout << "Machine Epsilon = " << epsilon() << endl;

		fivePointsinln();
		//cout << newRaph(3) << endl;
		//bisection(0,5);
		//forwardDiffx2();
		//cout << endl << "Trapezoidal: " << setprecision(17) << trapezoid(n, 1, 5) << endl;
		//cout << "    Simpson's: " << simpson(n, a, b) << endl;
		//cout << " Gauss-Quadrature: " << gauss3pt(a, b) << endl;
	//cout << "area=" << (   ( pow(e,-(3+0.000001)) - pow(e,-3) ) / 0.000001 ) << endl;;
		return 0;
	}

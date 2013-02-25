#include <ctime>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;

double powerResidue(double seed)
{

	ofstream of;
	of.open("powerResidue.csv");

	double r_i = seed;
	double a = 57, c = 1, M = 256;

	of << "i,defaultPRNG" << endl;
					///cout << "r_i = " << r_i/M << endl;


	for(int i=0; i<15000; i++)
	{
		r_i = fmod(( a*(r_i)+c),M);
					///cout << "r_i = " << r_i/M << endl;
		of << i << "," << r_i/M << endl;
	}

	of.close();
	return 0;
}


double dRand()
{
	ofstream of;
	of.open("drand48().csv");
	of << "i,drand48()" << endl;

	for(int i=0; i<15000; i++)
	{
		of << i << "," << drand48() << endl;
	}

	of.close();
	return 0;
}

double defaultPRNGk()
{
	double xi[20000][1];
	double xik[20000][1];
	double xRaisedK = 0;
	double k = 0.00000000001;
	double N = 10000;
	double autoCorr = 0;

	ofstream of;
	of.open("defaultPRNGk.csv");
	of << "x[i+k],x[i]" << endl;

	for(int i=0; i<N; i++)
	{
		xi[i][0] = (i);
		xi[i][1] = ((double)rand() / RAND_MAX);

		xRaisedK += pow(xi[i][1],k);

		xik[i][0] = (i+k);
	    xik[i][1] = ((double)rand() / RAND_MAX);

		autoCorr += (xi[i][1])*(xik[i][1])/N;
		cout << "xi[" << i << "][1] = " << xi[i][1] << endl;

	    of << xik[i][1] << "," << xi[i][1] << endl;
	}
	double kMoment = 1/(k+1);
	double calcKmoment = xRaisedK/N;
	cout << "k^th moment (calculated) = " << calcKmoment << endl;
	cout << "k^th moment (theoretical) = " << (1/(k+1)) << endl;
	cout << "Relative error: " << (fabs(calcKmoment-kMoment)/calcKmoment)*100 << "%" << endl;

	cout << endl << "Auto-correlation: " << autoCorr << endl;
	of.close();
	return 0;
}

double defaultPRNG()
{
	ofstream of;
	of.open("defaultPRNG.csv");
	of << "i,x[i]" << endl;

	for(double i=0; i<15000; i++)
	{

		of << i << "," << (double)rand() / RAND_MAX << endl;
				///cout << "i=" << i << "defaultPRNG="<< (double)rand() / RAND_MAX << endl;
	}
	return 0;
}

 double prngDouble(double a, double b)
 {
	 //Recall need to seed: "srand((unsigned)time(0))"

	 double c, d;

	 d = b-a;
	 c= ((double)rand()/(static_cast<double>(RAND_MAX)+1))*d+a;

	 return c;
 }

 double monteCarlo5D()
 {
	 double u1=1,u2=1,u3=1,u4=1,u5=1;
	 double l1=0,l2=0,l3=0,l4=0,l5=0;

	 double N = 4000;
	 double coeff5D = ( (u1-l1)*(u2-l2)*(u3-l3)*(u4-l4)*(u5-l5) ) / N;
	 cout << "coeff5D"
			 "=" << coeff5D << endl;
	 double doublePRNG1,doublePRNG2,doublePRNG3,doublePRNG4,doublePRNG5;
	 double summation;

	 	 for(int i=0; i<N; i++)
	 	 {
	 	doublePRNG1 = prngDouble(0,1);
	 	//cout << doublePRNG1 << endl;
	 	doublePRNG2 = prngDouble(0,1);
	 	//cout << doublePRNG2 << endl;
		doublePRNG3 = prngDouble(0,1);
		//cout << doublePRNG3 << endl;
		doublePRNG4 = prngDouble(0,1);
		//cout << doublePRNG4 << endl;
		doublePRNG5 = prngDouble(0,1);
		//cout << doublePRNG5 << endl << endl << endl;

	 	summation += pow((doublePRNG1+doublePRNG2+doublePRNG3+doublePRNG4+doublePRNG5),3);

	 	 }

	 double random5D = coeff5D * summation;
	 cout << "Summation = " << setprecision(10) << random5D << endl;
	 cout << "Relative error = " << fabs ( 100*(random5D - 18.75) / 18.75 ) << "%" << endl;
	 return 0;
 }

int main()
{
	monteCarlo5D();
	return 0;
}

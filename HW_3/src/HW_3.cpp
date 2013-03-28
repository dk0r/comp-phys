
#include <iostream>
#include <fstream>
#include <cmath>
#include <nr3.h>
#include <rk4.h>

using namespace std;


///Gravitational Constants
double Gsi  = 6.67398e-11; // (SI-units:  m^3 / kg*s^2)              [old# used: 6.67428e-11]
double Gastro = 1.48811382e-34; // (Astro-units: au^3 / kg*days^2)   [old# used: 1.48818071e-34]

//Celestial Mass  (SI-units: kg)
double Msun =    1.9891e30;
double Mmer =    0.3302e24;
double Mven =    4.8685e24;
double Mear =    5.9736e24;
double Mmar =    0.64185e24;
double Mjup = 1898.6e24; //actual: 1898.6e24
double Msat =  568.46e24;
double Mura =   86.832e24;
double Mplu =    0.0125e24;
double Mnep =  102.43e24;

//Solar distances @ Perihelion   //(Astro-units:  au)
double Pmer =  0.307491008;
double Pven =  0.718459424;
double Pear =  0.98323592;
double Pmar =  1.38116939;
double Pjup =  4.95007046;
double Psat =  9.04123831;
double Pura = 18.3244587;
double Pplu = 29.6583098;
double Pnep = 29.7093132;

//Celestial Perihelion Velocities (Astro-units:  au/days)
double Vmer =  34.0638003e-3;
double Vven =  20.364354e-3;
double Vear =  17.4939388e-3;
double Vmar =  15.3050307e-3;
double Vjup =   7.92396305e-3;
double Vsat =   5.87944197e-3;
double Vura =   4.10636861e-3;
double Vplu =   3.5230448e-3;
double Vnep =   3.1765158e-3;

//Gravitational Constants (Astro-units: au^3/days^2)
double GMsun = Gastro*Msun; //2.96014025e-4;
double GMmer = Gastro*Mmer;
double GMven = Gastro*Mven;
double GMear = Gastro*Mear; //8.88979629e-10;
double GMmar = Gastro*Mmar;
double GMjup = Gastro*Mjup; //2.82545990e-7;
double GMsat = Gastro*Msat;
double GMura = Gastro*Mura;
double GMplu = Gastro*Mplu;
double GMnep = Gastro*Mnep;

//Pi
double Pi = atan(1)*4;

double g = 9.8; //gravitational acceleration near earth's surface
//Mass1
double M1 = 1;   //mass of Mass1
double R1 = 1; //pendulum length
//Mass2
double M2 = 1;	 //mass of Mass1
double R2 = 1; //pendulum length

//Returns the radicand of the radius (x^2+y^2) between bodies
double R(double x, double y)


{
	double r;

	r = fabs(pow(x,2)+pow(y,2));

	return r;
}
//Returns the radicand raised to an exponent
double Re(double x, double y, double exp)

{
	double r;

	r = fabs(   pow( (pow(x,2)+pow(y,2)),exp)   );

	return r;
}



void pendDerivs(const Doub x, VecDoub_I & y, VecDoub_O & dydx)
{
	dydx[0] = y[1];
	dydx[1] = (  -g*(2*M1+M2)*sin(y[0])-M2*g*sin(y[0]-2*y[2])-2*sin(y[0]-y[2])*M2*(y[3]*y[3]*R2+y[1]*y[1]*R1*cos(y[0]-y[2]))  )
			/ (  R1*(2*M1+M2-M2*cos(2*y[0]-2*y[2]))  );


	dydx[2] = y[3];
	dydx[3] = (  2*sin(y[0]-y[2])*(y[1]*y[1]*R1*(M1+M2)+g*(M1+M2)*cos(y[0])+y[3]*y[3]*R2*M2*cos(y[0]-y[2]) ) )
			/ (  R2*(2*M1+M2-M2*cos(2*y[0]-2*y[2])) );


	/*
	dydx[0] = y[1];
	dydx[1] = ((-sin(y[0]-y[0]))*(y[1]*y[1]*A+y[3]*y[3])+w*(-2*sin(y[0])+sin(y[2])*A))/(2-A*A);

	dydx[2] = y[3];
	dydx[3] = ((-sin(y[0]-y[0]))*(2*y[1]*y[1]+y[3]*y[3]*A)+2*w*(-sin(y[2])+sin(y[0])*A))/(2-A*A);
	 */


	/*
	//////Theta-1
	//Angular Velocity
	dydx[0] = y[1];
	//Angular Acceleration
	dydx[1] = g*M1*R1*sin(y[0])+g*M2*R1*sin(y[0])+M2*R1*R2*sin(y[0]-y[2])*pow(y[3],2)
			+ M1*pow(R1,2)*;


	//////Theta-2
	//Angular Velocity
	dydx[2] = y[3];
	//Angular Acceleration
	dydx[3] =   (-GMsun * y[4]) / pow( pow(y[4],2) + pow(y[6],2) , 1.5 )
			  + (-GMear * y[4]) / pow( pow(y[0]-y[4],2) + pow(y[2]-y[6],2) , 1.5 );
	 */
}

void doublePend()
{
		VecDoub y(4);
		VecDoub dydx(4); //vector of positions & velocities for earth and
		VecDoub yout(4);

		double x;
		double xmin = 0;      //minimum starting position (units: au)
		double kmax=1000;  //max iterations (units: days)
		double h=0.01;       //time step size (units: days)

		///Initial Conditions:

		//Mass1
		y[0] = 0.17*atan(1)+0.008;       //Angular-Position
		y[1] = 0;        //Angular-Velocity

		//Mass2
		y[2] = 4*atan(1);        //Angular-Position
		y[3] = 0; 		 //Angular-Velocity


		pendDerivs(xmin, y, dydx);

		//file output streams
		ofstream ofPositionM1, ofPositionM2;
		ofPositionM1.open("x-y_position.M1.csv");
		ofPositionM2.open("x-y_position.M2.csv");


		for(int k=0; k < kmax; k++)
		{
			x=xmin+k*h;

			rk4(y, dydx,  x, h, yout, pendDerivs);

			/*
			//display output
			cout << "k = " << k
				 << "    Xe = " << yout[0] << "    X'e= " << yout[1]
				 << "    Ye = " << yout[2] << "    Y'e= " << yout[3] << endl

				 << "    Xj = " << yout[4] << "    X'j= " << yout[5]
				 << "    Yj = " << yout[6] << "    Y'j= " << yout[7] << endl << endl;
			 */

			//file output streams
			ofPositionM1 << R1*sin(y[0]) << "," << -R1*cos(y[0]) << endl;
			ofPositionM2 << R1*sin(y[0]) + R2*sin(y[2]) << "," << (-R1*cos(y[0]) - R2*cos(y[2])) << endl;


			y = yout;

			pendDerivs(x,y,dydx);
		}

		//closes file output streams
		ofPositionM1.close();
		ofPositionM2.close();
}




//Two-body Earth/Jupiter about Sun
void derivs(const Doub x, VecDoub_I & y, VecDoub_O & dydx)
{
	///Earth

	//X-components
	dydx[0] = y[1];
	dydx[1] = (-GMsun * y[0]) / pow( pow(y[0],2) + pow(y[2],2) , 1.5 )
			+ (-GMjup * y[0]) / pow( pow(y[0]-y[4],2) + pow(y[2]-y[6],2) , 1.5 );
	//Y-components
	dydx[2] = y[3];
	dydx[3] = (-GMsun * y[2]) / pow( pow(y[0],2) + pow(y[2],2) , 1.5 )
			+ (-GMjup * y[2]) / pow( pow(y[0]-y[4],2) + pow(y[2]-y[6],2) , 1.5 );


	///Jupiter

	//X-components
	dydx[4] = y[5];
	dydx[5] =   (-GMsun * y[4]) / pow( pow(y[4],2) + pow(y[6],2) , 1.5 )
			  + (-GMear * y[4]) / pow( pow(y[0]-y[4],2) + pow(y[2]-y[6],2) , 1.5 );
	//Y-components
	dydx[6] = y[7];
	dydx[7] =   (-GMsun * y[6]) / pow( pow(y[4],2) + pow(y[6],2) , 1.5 )
			  + (-GMear * y[6]) / pow( pow(y[0]-y[4],2) + pow(y[2]-y[6],2) , 1.5 );
}

int sunEarthJupiter()
{
	//cout << "GMs =" << GMs << ",  G_astro*M_s = " << G_astro*M_s << endl;
	//cout << "GMe =" << GMe << ",  G_astro*M_e = " << G_astro*M_e << endl;
	//cout << "GMj =" << GMj << ",  G_astro*M_j = " << G_astro*M_j << endl;

	VecDoub y(8);
	VecDoub dydx(8); //vector of positions & velocities for earth and
	VecDoub yout(8);

	double x;
	double xmin;      //minimum starting position (units: au)
	double kmax=10*2*365;  //max iterations (units: days)
	double h=0.5;       //time step size (units: days)

	xmin = 0;

	///Initial Conditions:

	//Earth
	y[0] = 0.98329;       //X-Position          (units: au)
	y[1] = 0;                      //X-Velocity (units: au/day)
	y[2] = 0;             //Y- Position         (units: au)
	y[3] = 1.74939488e-2; 		   //Y-Velocity (units: au/day)

	//Jupiter
	y[4] = 4.950429;      //X-Position          (units: au)
	y[5] = 0;                      //X-Velocity (units: au/day)
	y[6] = 0;             //Y- Position         (units: au)
	y[7] = 7.92396305e-3; 		   //Y-Velocity (units: au/day)  !!!!!!!! e-6 original


	derivs(xmin, y, dydx);

	//file output streams
	ofstream ofAllData, ofPositionEarth, ofPositionJupiter;
	ofAllData.open("x-y_velocity.Earth.csv");
	ofPositionEarth.open("x-y_position.Earth.csv");
	ofPositionJupiter.open("x-y_position.Jupiter.csv");

	//variables for verifying Keppler's 2nd law
	double area;
	double oldx;
	double oldy;
	double radius;
	double velocity;
	double analytic;

	for(int k=0; k < kmax; k++)
	{
			x=xmin+k*h;

			rk4(y, dydx,  x, h, yout, derivs);

			/*
			//display output
			cout << "k = " << k
				 << "    Xe = " << yout[0] << "    X'e= " << yout[1]
				 << "    Ye = " << yout[2] << "    Y'e= " << yout[3] << endl

				 << "    Xj = " << yout[4] << "    X'j= " << yout[5]
				 << "    Yj = " << yout[6] << "    Y'j= " << yout[7] << endl << endl;
			 */
			//file output streams
			ofAllData << yout[1] << "," << yout[3] << endl; // << yout[4] << "," << yout[5] << "," << yout[6] << ","  << yout[7] << endl;
			ofPositionEarth << yout[0] << "," << yout[2] << endl;
			ofPositionJupiter << yout[4] << "," << yout[6] << endl;



			/*
			//Verifies Keppler's 2nd Law
			//Outputs areas swept out in 1st and 52nd week of Earth's orbit
			if(k==0 || k==338)
			{
			oldx = y[0];
			oldy = y[2];
			}

			if(k==13 || k==351)
			{
				cout << " a=" << oldx << "  d=" << y[2] << "  b=" << oldy << "  c=" << y[0] << endl;
				area =  fabs( oldx*y[2]-oldy*y[0] ) / 2;
				cout << "Calculated Area at " << k*h/7 << "-weeks =" << area << endl;

				radius = Re(y[0],y[2],0.5);
				cout << "radius=" << radius << endl;
				velocity = Re(y[1],y[3],0.5);
				cout << "velocity=" << velocity << endl;
				analytic = radius*velocity/2;
				cout << "Analytical Area = " << analytic << endl;
				cout <<"Error = " << (100*(area-analytic))/analytic << "%" << endl << endl;
			}
			 */


			y = yout;

			derivs(x,y,dydx);
	}

	//closes file output streams
	ofAllData.close();
	ofPositionEarth.close();
	ofPositionJupiter.close();

return 0;

}




void GRmercDerivs(const Doub x, VecDoub_I & y, VecDoub_O & dydx)

{
	//exponents for radius calculation in R(i,j,exp)
	double exp1 = 1.5;
	double exp2 = 2;

	//coefficient of relativity correction
	double a = 0.0001;//6.4e-7; //units: au^2 ?


///Mercury------------------------------------------------------------------------------------

	//Indices for Mercury's x and y positions in y[]
	int i = 0;
	int j = 2;

	//Radicand's of radius (x^2+y^2) between section and subscript bodies.
	 double Rsun = R(y[i],y[j]);

	//X-components
	dydx[i]   = y[i+1];
	dydx[i+1] = (1+a/pow(Rsun,exp2)) * (-GMsun * y[i]) / pow(Rsun,exp1);

	//Y-components
	dydx[j]   = y[j+1];
	dydx[j+1] = (1+a/pow(Rsun,exp2)) * (-GMsun * y[j]) / pow(Rsun,exp1);

}

int GRmercury()
{
		cout << "Started GR.Mercury -- Started GR.Mercury -- Started GR.Mercury" << "\n";


		VecDoub y(36);
		VecDoub dydx(36); //vector of positions & velocities for earth and
		VecDoub yout(36);


		double x;
		double xmin = 0;       //minimum starting position (units: au)
		double kmax = 20*87.969;  //max iterations (units: days)
		double h = 0.05;          //time step size (units: days)

		double constraint = (kmax-1518)/kmax;


		///Initial Conditions:

		//Mercury
		y[0] = Pmer;   //X-Position          (units: au)
		y[1] = 0;                      //X-Velocity (units: au/day)
		y[2] = 0;             //Y- Position         (units: au)
		y[3] = Vmer; 		   //Y-Velocity (units: au/day)

		GRmercDerivs(xmin, y, dydx);

			//file output stream
			ofstream ofPositionMercury;
			ofPositionMercury.open("40.csv");

		for(int k=0; k < kmax; k++)
			{
					x=xmin+k*h;

					rk4(y, dydx,  x, h, yout, GRmercDerivs);


					///debug output
					//cout << "k = " << k
					//	 << "    Xm = " << yout[0] << "    X'm= " << yout[1]
					//	 << "    Ym = " << yout[2] << "    Y'm= " << yout[3] << endl;


					//file output stream
					//if(k/kmax >= constraint)
					//{
					ofPositionMercury << yout[0] << "," << yout[2] << "\n";
					//}
					y = yout;

					GRmercDerivs(x,y,dydx);
			}

			//closes file output stream
			ofPositionMercury.close();


		cout << "\n" << "CoMpLeTe CoMpLeTe CoMpLeTe CoMpLeTe CoMpLeTe CoMpLeTe" << "\n";

	return 0;
}





int main()
{
	doublePend();
	//GRmercury();
	//sunEarthJupiter();
	return 0;
}














////////////////////THE WHOLE SOLAR SYSTEM////////////////////////

//////Relativistic orbit of Mercury in Solar System
/*
void GRmercDerivs(const Doub x, VecDoub_I & y, VecDoub_O & dydx)

{
	//exponents for radius calculation in R(i,j,exp)
	double exp1 = 1.5;
	double exp2 = 2;

	//coefficient of relativity correction
	double a = 0.01; //units: au^2 ?


///Mercury------------------------------------------------------------------------------------

	//Mercury x and y indices for y[]
	int i = 0;
	int j = 2;


		//Radicand's of radius (x^2+y^2) between section and subscript bodies.
		 double Rsun = R(y[i]       , y[j]);
		 double Rmer = R(y[i]-y[0]  , y[j]-y[2]);
		 double Rven = R(y[i]-y[4]  , y[j]-y[6]);
		 double Rear = R(y[i]-y[8]  , y[j]-y[10]);
		 double Rmar = R(y[i]-y[12] , y[j]-y[14]);
		 double Rjup = R(y[i]-y[16] , y[j]-y[18]);
		 double Rsat = R(y[i]-y[20] , y[j]-y[22]);
		 double Rura = R(y[i]-y[24] , y[j]-y[26]);
		 double Rplu = R(y[i]-y[28] , y[j]-y[30]);
		 double Rnep = R(y[i]-y[32] , y[j]-y[34]);

	//X-components
	dydx[i]   = y[i+1];
	dydx[i+1] = (1+a/pow(Rsun,exp2)) * (-GMsun * y[i]) / pow(Rsun,exp1)

			  + (1+a/pow(Rven,exp2)) * (-GMven * y[i]) / pow(Rven,exp1)
			  + (1+a/pow(Rear,exp2)) * (-GMear * y[i]) / pow(Rear,exp1)
			  + (1+a/pow(Rmar,exp2)) * (-GMmar * y[i]) / pow(Rmar,exp1)
			  + (1+a/pow(Rjup,exp2)) * (-GMjup * y[i]) / pow(Rjup,exp1)
			  + (1+a/pow(Rsat,exp2)) * (-GMsat * y[i]) / pow(Rsat,exp1)
			  + (1+a/pow(Rura,exp2)) * (-GMura * y[i]) / pow(Rura,exp1)
			  + (1+a/pow(Rplu,exp2)) * (-GMplu * y[i]) / pow(Rplu,exp1)
			  + (1+a/pow(Rnep,exp2)) * (-GMnep * y[i]) / pow(Rnep,exp1);

	//Y-components
	dydx[j]   = y[j+1];
	dydx[j+1] = (1+a/pow(Rsun,exp2)) * (-GMsun * y[j]) / pow(Rsun,exp1)

			  + (1+a/pow(Rven,exp2)) * (-GMven * y[j]) / pow(Rven,exp1)
			  + (1+a/pow(Rear,exp2)) * (-GMear * y[j]) / pow(Rear,exp1)
			  + (1+a/pow(Rmar,exp2)) * (-GMmar * y[j]) / pow(Rmar,exp1)
			  + (1+a/pow(Rjup,exp2)) * (-GMjup * y[j]) / pow(Rjup,exp1)
			  + (1+a/pow(Rsat,exp2)) * (-GMsat * y[j]) / pow(Rsat,exp1)
			  + (1+a/pow(Rura,exp2)) * (-GMura * y[j]) / pow(Rura,exp1)
			  + (1+a/pow(Rplu,exp2)) * (-GMplu * y[j]) / pow(Rplu,exp1)
			  + (1+a/pow(Rnep,exp2)) * (-GMnep * y[j]) / pow(Rnep,exp1);




///Venus--------------------------------------------------------------------------------------

	//position indices
	i = 4;
	j = 6;

	//Radicand's of radius (x^2+y^2) between section and subscript bodies.
	 Rsun = R(y[i]       , y[j]);
	 Rmer = R(y[i]-y[0]  , y[j]-y[2]);

	 Rear = R(y[i]-y[8]  , y[j]-y[10]);
	 Rmar = R(y[i]-y[12] , y[j]-y[14]);
	 Rjup = R(y[i]-y[16] , y[j]-y[18]);
	 Rsat = R(y[i]-y[20] , y[j]-y[22]);
	 Rura = R(y[i]-y[24] , y[j]-y[26]);
	 Rplu = R(y[i]-y[28] , y[j]-y[30]);
	 Rnep = R(y[i]-y[32] , y[j]-y[34]);

	//X-components
	dydx[i]   = y[i+1];
	dydx[i+1] = (1+a/pow(Rsun,exp2)) * (-GMsun * y[i]) / pow(Rsun,exp1)
		      + (1+a/pow(Rmer,exp2)) * (-GMmer * y[i]) / pow(Rmer,exp1)

		      + (1+a/pow(Rear,exp2)) * (-GMear * y[i]) / pow(Rear,exp1)
		      + (1+a/pow(Rmar,exp2)) * (-GMmar * y[i]) / pow(Rmar,exp1)
		      + (1+a/pow(Rjup,exp2)) * (-GMjup * y[i]) / pow(Rjup,exp1)
		      + (1+a/pow(Rsat,exp2)) * (-GMsat * y[i]) / pow(Rsat,exp1)
		      + (1+a/pow(Rura,exp2)) * (-GMura * y[i]) / pow(Rura,exp1)
		      + (1+a/pow(Rplu,exp2)) * (-GMplu * y[i]) / pow(Rplu,exp1)
		      + (1+a/pow(Rnep,exp2)) * (-GMnep * y[i]) / pow(Rnep,exp1);

	//Y-components
	dydx[j]   = y[j+1];
	dydx[j+1] = (1+a/pow(Rsun,exp2)) * (-GMsun * y[j]) / pow(Rsun,exp1)
			  + (1+a/pow(Rmer,exp2)) * (-GMmer * y[j]) / pow(Rmer,exp1)

			  + (1+a/pow(Rear,exp2)) * (-GMear * y[j]) / pow(Rear,exp1)
			  + (1+a/pow(Rmar,exp2)) * (-GMmar * y[j]) / pow(Rmar,exp1)
			  + (1+a/pow(Rjup,exp2)) * (-GMjup * y[j]) / pow(Rjup,exp1)
			  + (1+a/pow(Rsat,exp2)) * (-GMsat * y[j]) / pow(Rsat,exp1)
			  + (1+a/pow(Rura,exp2)) * (-GMura * y[j]) / pow(Rura,exp1)
			  + (1+a/pow(Rplu,exp2)) * (-GMplu * y[j]) / pow(Rplu,exp1)
			  + (1+a/pow(Rnep,exp2)) * (-GMnep * y[j]) / pow(Rnep,exp1);

///Earth---------------------------------------------------------------------------------------
i = 8;
j = 10;


//Radicand's of radius (x^2+y^2) between section and subscript bodies.
 Rsun = R(y[i]       , y[j]);
 Rmer = R(y[i]-y[0]  , y[j]-y[2]);
 Rven = R(y[i]-y[4]  , y[j]-y[6]);

 Rmar = R(y[i]-y[12] , y[j]-y[14]);
 Rjup = R(y[i]-y[16] , y[j]-y[18]);
 Rsat = R(y[i]-y[20] , y[j]-y[22]);
 Rura = R(y[i]-y[24] , y[j]-y[26]);
 Rplu = R(y[i]-y[28] , y[j]-y[30]);
 Rnep = R(y[i]-y[32] , y[j]-y[34]);

//X-components
dydx[i]   = y[i+1];
dydx[i+1] = (1+a/pow(Rsun,exp2)) * (-GMsun * y[i]) / pow(Rsun,exp1)
	      + (1+a/pow(Rmer,exp2)) * (-GMmer * y[i]) / pow(Rmer,exp1)
	      + (1+a/pow(Rven,exp2)) * (-GMven * y[i]) / pow(Rven,exp1)

	      + (1+a/pow(Rmar,exp2)) * (-GMmar * y[i]) / pow(Rmar,exp1)
	      + (1+a/pow(Rjup,exp2)) * (-GMjup * y[i]) / pow(Rjup,exp1)
	      + (1+a/pow(Rsat,exp2)) * (-GMsat * y[i]) / pow(Rsat,exp1)
	      + (1+a/pow(Rura,exp2)) * (-GMura * y[i]) / pow(Rura,exp1)
	      + (1+a/pow(Rplu,exp2)) * (-GMplu * y[i]) / pow(Rplu,exp1)
	      + (1+a/pow(Rnep,exp2)) * (-GMnep * y[i]) / pow(Rnep,exp1);

//Y-components
dydx[j]   = y[j+1];
dydx[j+1] = (1+a/pow(Rsun,exp2)) * (-GMsun * y[j]) / pow(Rsun,exp1)
		  + (1+a/pow(Rmer,exp2)) * (-GMmer * y[j]) / pow(Rmer,exp1)
		  + (1+a/pow(Rven,exp2)) * (-GMven * y[j]) / pow(Rven,exp1)

		  + (1+a/pow(Rmar,exp2)) * (-GMmar * y[j]) / pow(Rmar,exp1)
		  + (1+a/pow(Rjup,exp2)) * (-GMjup * y[j]) / pow(Rjup,exp1)
		  + (1+a/pow(Rsat,exp2)) * (-GMsat * y[j]) / pow(Rsat,exp1)
		  + (1+a/pow(Rura,exp2)) * (-GMura * y[j]) / pow(Rura,exp1)
		  + (1+a/pow(Rplu,exp2)) * (-GMplu * y[j]) / pow(Rplu,exp1)
		  + (1+a/pow(Rnep,exp2)) * (-GMnep * y[j]) / pow(Rnep,exp1);

///Mars-----------------------------------------------------------------------------------------
i = 12;
j = 14;



//Radicand's of radius (x^2+y^2) between section and subscript bodies.
 Rsun = R(y[i]       , y[j]);
 Rmer = R(y[i]-y[0]  , y[j]-y[2]);
 Rven = R(y[i]-y[4]  , y[j]-y[6]);
 Rear = R(y[i]-y[8]  , y[j]-y[10]);

 Rjup = R(y[i]-y[16] , y[j]-y[18]);
 Rsat = R(y[i]-y[20] , y[j]-y[22]);
 Rura = R(y[i]-y[24] , y[j]-y[26]);
 Rplu = R(y[i]-y[28] , y[j]-y[30]);
 Rnep = R(y[i]-y[32] , y[j]-y[34]);

//X-components
dydx[i]   = y[i+1];
dydx[i+1] = (1+a/pow(Rsun,exp2)) * (-GMsun * y[i]) / pow(Rsun,exp1)
	      + (1+a/pow(Rmer,exp2)) * (-GMmer * y[i]) / pow(Rmer,exp1)
	      + (1+a/pow(Rven,exp2)) * (-GMven * y[i]) / pow(Rven,exp1)
	      + (1+a/pow(Rear,exp2)) * (-GMear * y[i]) / pow(Rear,exp1)

	      + (1+a/pow(Rjup,exp2)) * (-GMjup * y[i]) / pow(Rjup,exp1)
	      + (1+a/pow(Rsat,exp2)) * (-GMsat * y[i]) / pow(Rsat,exp1)
	      + (1+a/pow(Rura,exp2)) * (-GMura * y[i]) / pow(Rura,exp1)
	      + (1+a/pow(Rplu,exp2)) * (-GMplu * y[i]) / pow(Rplu,exp1)
	      + (1+a/pow(Rnep,exp2)) * (-GMnep * y[i]) / pow(Rnep,exp1);

//Y-components
dydx[j]   = y[j+1];
dydx[j+1] = (1+a/pow(Rsun,exp2)) * (-GMsun * y[j]) / pow(Rsun,exp1)
		  + (1+a/pow(Rmer,exp2)) * (-GMmer * y[j]) / pow(Rmer,exp1)
		  + (1+a/pow(Rven,exp2)) * (-GMven * y[j]) / pow(Rven,exp1)
		  + (1+a/pow(Rear,exp2)) * (-GMear * y[j]) / pow(Rear,exp1)

		  + (1+a/pow(Rjup,exp2)) * (-GMjup * y[j]) / pow(Rjup,exp1)
		  + (1+a/pow(Rsat,exp2)) * (-GMsat * y[j]) / pow(Rsat,exp1)
		  + (1+a/pow(Rura,exp2)) * (-GMura * y[j]) / pow(Rura,exp1)
		  + (1+a/pow(Rplu,exp2)) * (-GMplu * y[j]) / pow(Rplu,exp1)
		  + (1+a/pow(Rnep,exp2)) * (-GMnep * y[j]) / pow(Rnep,exp1);

///Jupiter---------------------------------------------------------------------------------------
i = 16;
j = 18;



//Radicand's of radius (x^2+y^2) between section and subscript bodies.
 Rsun = R(y[i]       , y[j]);
 Rmer = R(y[i]-y[0]  , y[j]-y[2]);
 Rven = R(y[i]-y[4]  , y[j]-y[6]);
 Rear = R(y[i]-y[8]  , y[j]-y[10]);
 Rmar = R(y[i]-y[12] , y[j]-y[14]);

 Rsat = R(y[i]-y[20] , y[j]-y[22]);
 Rura = R(y[i]-y[24] , y[j]-y[26]);
 Rplu = R(y[i]-y[28] , y[j]-y[30]);
 Rnep = R(y[i]-y[32] , y[j]-y[34]);

//X-components
dydx[i]   = y[i+1];
dydx[i+1] = (1+a/pow(Rsun,exp2)) * (-GMsun * y[i]) / pow(Rsun,exp1)
	      + (1+a/pow(Rmer,exp2)) * (-GMmer * y[i]) / pow(Rmer,exp1)
	      + (1+a/pow(Rven,exp2)) * (-GMven * y[i]) / pow(Rven,exp1)
	      + (1+a/pow(Rear,exp2)) * (-GMear * y[i]) / pow(Rear,exp1)
	      + (1+a/pow(Rmar,exp2)) * (-GMmar * y[i]) / pow(Rmar,exp1)

	      + (1+a/pow(Rsat,exp2)) * (-GMsat * y[i]) / pow(Rsat,exp1)
	      + (1+a/pow(Rura,exp2)) * (-GMura * y[i]) / pow(Rura,exp1)
	      + (1+a/pow(Rplu,exp2)) * (-GMplu * y[i]) / pow(Rplu,exp1)
	      + (1+a/pow(Rnep,exp2)) * (-GMnep * y[i]) / pow(Rnep,exp1);

//Y-components
dydx[j]   = y[j+1];
dydx[j+1] = (1+a/pow(Rsun,exp2)) * (-GMsun * y[j]) / pow(Rsun,exp1)
		  + (1+a/pow(Rmer,exp2)) * (-GMmer * y[j]) / pow(Rmer,exp1)
		  + (1+a/pow(Rven,exp2)) * (-GMven * y[j]) / pow(Rven,exp1)
		  + (1+a/pow(Rear,exp2)) * (-GMear * y[j]) / pow(Rear,exp1)
		  + (1+a/pow(Rmar,exp2)) * (-GMmar * y[j]) / pow(Rmar,exp1)

		  + (1+a/pow(Rsat,exp2)) * (-GMsat * y[j]) / pow(Rsat,exp1)
		  + (1+a/pow(Rura,exp2)) * (-GMura * y[j]) / pow(Rura,exp1)
		  + (1+a/pow(Rplu,exp2)) * (-GMplu * y[j]) / pow(Rplu,exp1)
		  + (1+a/pow(Rnep,exp2)) * (-GMnep * y[j]) / pow(Rnep,exp1);

///Saturn---------------------------------------------------------------------------------------
i = 20;
j = 22;



//Radicand's of radius (x^2+y^2) between section and subscript bodies.
 Rsun = R(y[i]       , y[j]);
 Rmer = R(y[i]-y[0]  , y[j]-y[2]);
 Rven = R(y[i]-y[4]  , y[j]-y[6]);
 Rear = R(y[i]-y[8]  , y[j]-y[10]);
 Rmar = R(y[i]-y[12] , y[j]-y[14]);
 Rjup = R(y[i]-y[16] , y[j]-y[18]);

 Rura = R(y[i]-y[24] , y[j]-y[26]);
 Rplu = R(y[i]-y[28] , y[j]-y[30]);
 Rnep = R(y[i]-y[32] , y[j]-y[34]);

//X-components
dydx[i]   = y[i+1];
dydx[i+1] = (1+a/pow(Rsun,exp2)) * (-GMsun * y[i]) / pow(Rsun,exp1)
	      + (1+a/pow(Rmer,exp2)) * (-GMmer * y[i]) / pow(Rmer,exp1)
	      + (1+a/pow(Rven,exp2)) * (-GMven * y[i]) / pow(Rven,exp1)
	      + (1+a/pow(Rear,exp2)) * (-GMear * y[i]) / pow(Rear,exp1)
	      + (1+a/pow(Rmar,exp2)) * (-GMmar * y[i]) / pow(Rmar,exp1)
	      + (1+a/pow(Rjup,exp2)) * (-GMjup * y[i]) / pow(Rjup,exp1)

	      + (1+a/pow(Rura,exp2)) * (-GMura * y[i]) / pow(Rura,exp1)
	      + (1+a/pow(Rplu,exp2)) * (-GMplu * y[i]) / pow(Rplu,exp1)
	      + (1+a/pow(Rnep,exp2)) * (-GMnep * y[i]) / pow(Rnep,exp1);

//Y-components
dydx[j]   = y[j+1];
dydx[j+1] = (1+a/pow(Rsun,exp2)) * (-GMsun * y[j]) / pow(Rsun,exp1)
		  + (1+a/pow(Rmer,exp2)) * (-GMmer * y[j]) / pow(Rmer,exp1)
		  + (1+a/pow(Rven,exp2)) * (-GMven * y[j]) / pow(Rven,exp1)
		  + (1+a/pow(Rear,exp2)) * (-GMear * y[j]) / pow(Rear,exp1)
		  + (1+a/pow(Rmar,exp2)) * (-GMmar * y[j]) / pow(Rmar,exp1)
		  + (1+a/pow(Rjup,exp2)) * (-GMjup * y[j]) / pow(Rjup,exp1)

		  + (1+a/pow(Rura,exp2)) * (-GMura * y[j]) / pow(Rura,exp1)
		  + (1+a/pow(Rplu,exp2)) * (-GMplu * y[j]) / pow(Rplu,exp1)
		  + (1+a/pow(Rnep,exp2)) * (-GMnep * y[j]) / pow(Rnep,exp1);

///Uranus---------------------------------------------------------------------------------------
i = 24;
j = 26;


//Radicand's of radius (x^2+y^2) between section and subscript bodies.
 Rsun = R(y[i]       , y[j]);
 Rmer = R(y[i]-y[0]  , y[j]-y[2]);
 Rven = R(y[i]-y[4]  , y[j]-y[6]);
 Rear = R(y[i]-y[8]  , y[j]-y[10]);
 Rmar = R(y[i]-y[12] , y[j]-y[14]);
 Rjup = R(y[i]-y[16] , y[j]-y[18]);
 Rsat = R(y[i]-y[20] , y[j]-y[22]);

 Rplu = R(y[i]-y[28] , y[j]-y[30]);
 Rnep = R(y[i]-y[32] , y[j]-y[34]);

//X-components
dydx[i]   = y[i+1];
dydx[i+1] = (1+a/pow(Rsun,exp2)) * (-GMsun * y[i]) / pow(Rsun,exp1)
	      + (1+a/pow(Rmer,exp2)) * (-GMmer * y[i]) / pow(Rmer,exp1)
	      + (1+a/pow(Rven,exp2)) * (-GMven * y[i]) / pow(Rven,exp1)
	      + (1+a/pow(Rear,exp2)) * (-GMear * y[i]) / pow(Rear,exp1)
	      + (1+a/pow(Rmar,exp2)) * (-GMmar * y[i]) / pow(Rmar,exp1)
	      + (1+a/pow(Rjup,exp2)) * (-GMjup * y[i]) / pow(Rjup,exp1)
	      + (1+a/pow(Rsat,exp2)) * (-GMsat * y[i]) / pow(Rsat,exp1)

	      + (1+a/pow(Rplu,exp2)) * (-GMplu * y[i]) / pow(Rplu,exp1)
	      + (1+a/pow(Rnep,exp2)) * (-GMnep * y[i]) / pow(Rnep,exp1);

//Y-components
dydx[j]   = y[j+1];
dydx[j+1] = (1+a/pow(Rsun,exp2)) * (-GMsun * y[j]) / pow(Rsun,exp1)
		  + (1+a/pow(Rmer,exp2)) * (-GMmer * y[j]) / pow(Rmer,exp1)
		  + (1+a/pow(Rven,exp2)) * (-GMven * y[j]) / pow(Rven,exp1)
		  + (1+a/pow(Rear,exp2)) * (-GMear * y[j]) / pow(Rear,exp1)
		  + (1+a/pow(Rmar,exp2)) * (-GMmar * y[j]) / pow(Rmar,exp1)
		  + (1+a/pow(Rjup,exp2)) * (-GMjup * y[j]) / pow(Rjup,exp1)
		  + (1+a/pow(Rsat,exp2)) * (-GMsat * y[j]) / pow(Rsat,exp1)

		  + (1+a/pow(Rplu,exp2)) * (-GMplu * y[j]) / pow(Rplu,exp1)
		  + (1+a/pow(Rnep,exp2)) * (-GMnep * y[j]) / pow(Rnep,exp1);

///Pluto---------------------------------------------------------------------------------------
	i = 28;
	j = 30;



	//Radicand's of radius (x^2+y^2) between section and subscript bodies.
	 Rsun = R(y[i]       , y[j]);
	 Rmer = R(y[i]-y[0]  , y[j]-y[2]);
	 Rven = R(y[i]-y[4]  , y[j]-y[6]);
	 Rear = R(y[i]-y[8]  , y[j]-y[10]);
	 Rmar = R(y[i]-y[12] , y[j]-y[14]);
	 Rjup = R(y[i]-y[16] , y[j]-y[18]);
	 Rsat = R(y[i]-y[20] , y[j]-y[22]);
	 Rura = R(y[i]-y[24] , y[j]-y[26]);

	 Rnep = R(y[i]-y[32] , y[j]-y[34]);

	//X-components
	dydx[i]   = y[i+1];
	dydx[i+1] = (1+a/pow(Rsun,exp2)) * (-GMsun * y[i]) / pow(Rsun,exp1)
		      + (1+a/pow(Rmer,exp2)) * (-GMmer * y[i]) / pow(Rmer,exp1)
		      + (1+a/pow(Rven,exp2)) * (-GMven * y[i]) / pow(Rven,exp1)
		      + (1+a/pow(Rear,exp2)) * (-GMear * y[i]) / pow(Rear,exp1)
		      + (1+a/pow(Rmar,exp2)) * (-GMmar * y[i]) / pow(Rmar,exp1)
		      + (1+a/pow(Rjup,exp2)) * (-GMjup * y[i]) / pow(Rjup,exp1)
		      + (1+a/pow(Rsat,exp2)) * (-GMsat * y[i]) / pow(Rsat,exp1)
		      + (1+a/pow(Rura,exp2)) * (-GMura * y[i]) / pow(Rura,exp1)

		      + (1+a/pow(Rnep,exp2)) * (-GMnep * y[i]) / pow(Rnep,exp1);

	//Y-components
	dydx[j]   = y[j+1];
	dydx[j+1] = (1+a/pow(Rsun,exp2)) * (-GMsun * y[j]) / pow(Rsun,exp1)
			  + (1+a/pow(Rmer,exp2)) * (-GMmer * y[j]) / pow(Rmer,exp1)
			  + (1+a/pow(Rven,exp2)) * (-GMven * y[j]) / pow(Rven,exp1)
			  + (1+a/pow(Rear,exp2)) * (-GMear * y[j]) / pow(Rear,exp1)
			  + (1+a/pow(Rmar,exp2)) * (-GMmar * y[j]) / pow(Rmar,exp1)
			  + (1+a/pow(Rjup,exp2)) * (-GMjup * y[j]) / pow(Rjup,exp1)
			  + (1+a/pow(Rsat,exp2)) * (-GMsat * y[j]) / pow(Rsat,exp1)
			  + (1+a/pow(Rura,exp2)) * (-GMura * y[j]) / pow(Rura,exp1)

			  + (1+a/pow(Rnep,exp2)) * (-GMnep * y[j]) / pow(Rnep,exp1);

///Neptune---------------------------------------------------------------------------------------
	i = 32;
	j = 34;

	//Radicand's of radius (x^2+y^2) between section and subscript bodies.
	 Rsun = R(y[i]       , y[j]);
	 Rmer = R(y[i]-y[0]  , y[j]-y[2]);
	 Rven = R(y[i]-y[4]  , y[j]-y[6]);
	 Rear = R(y[i]-y[8]  , y[j]-y[10]);
	 Rmar = R(y[i]-y[12] , y[j]-y[14]);
	 Rjup = R(y[i]-y[16] , y[j]-y[18]);
	 Rsat = R(y[i]-y[20] , y[j]-y[22]);
	 Rura = R(y[i]-y[24] , y[j]-y[26]);
	 Rplu = R(y[i]-y[28] , y[j]-y[30]);
	 //

	//X-components
	dydx[i]   = y[i+1];
	dydx[i+1] = (1+a/pow(Rsun,exp2)) * (-GMsun * y[i]) / pow(Rsun,exp1)
		      + (1+a/pow(Rmer,exp2)) * (-GMmer * y[i]) / pow(Rmer,exp1)
		      + (1+a/pow(Rven,exp2)) * (-GMven * y[i]) / pow(Rven,exp1)
		      + (1+a/pow(Rear,exp2)) * (-GMear * y[i]) / pow(Rear,exp1)
		      + (1+a/pow(Rmar,exp2)) * (-GMmar * y[i]) / pow(Rmar,exp1)
		      + (1+a/pow(Rjup,exp2)) * (-GMjup * y[i]) / pow(Rjup,exp1)
		      + (1+a/pow(Rsat,exp2)) * (-GMsat * y[i]) / pow(Rsat,exp1)
		      + (1+a/pow(Rura,exp2)) * (-GMura * y[i]) / pow(Rura,exp1)
		      + (1+a/pow(Rplu,exp2)) * (-GMplu * y[i]) / pow(Rplu,exp1);
		      //

	//Y-components
	dydx[j]   = y[j+1];
	dydx[j+1] = (1+a/pow(Rsun,exp2)) * (-GMsun * y[j]) / pow(Rsun,exp1)
			  + (1+a/pow(Rmer,exp2)) * (-GMmer * y[j]) / pow(Rmer,exp1)
			  + (1+a/pow(Rven,exp2)) * (-GMven * y[j]) / pow(Rven,exp1)
			  + (1+a/pow(Rear,exp2)) * (-GMear * y[j]) / pow(Rear,exp1)
			  + (1+a/pow(Rmar,exp2)) * (-GMmar * y[j]) / pow(Rmar,exp1)
			  + (1+a/pow(Rjup,exp2)) * (-GMjup * y[j]) / pow(Rjup,exp1)
			  + (1+a/pow(Rsat,exp2)) * (-GMsat * y[j]) / pow(Rsat,exp1)
			  + (1+a/pow(Rura,exp2)) * (-GMura * y[j]) / pow(Rura,exp1)
			  + (1+a/pow(Rplu,exp2)) * (-GMplu * y[j]) / pow(Rplu,exp1);
			  //
}
*/

/*
int GRmercury()
{
		cout << "Started GR.Mercury -- GR.Started Mercury -- GR.Started Mercury" << "\n";


		VecDoub y(36);
		VecDoub dydx(36); //vector of positions & velocities for earth and
		VecDoub yout(36);


		double x;
		double xmin = 0;       //minimum starting position (units: au)
		double kmax = 226*87.969;  //max iterations (units: days)
		double h = .0005;          //time step size (units: days)


		///Initial Conditions:

		//Mercury
		y[0] = Pmer;   //X-Position          (units: au)
		y[1] = 0;                      //X-Velocity (units: au/day)
		y[2] = 0;             //Y- Position         (units: au)
		y[3] = Vmer; 		   //Y-Velocity (units: au/day)

		//Venus
		y[4] = Pven;   //X-Position          (units: au)
		y[5] = 0;                      //X-Velocity (units: au/day)
		y[6] = 0;             //Y- Position         (units: au)
		y[7] = Vven; 		   //Y-Velocity (units: au/day)

		//Earth
		y[8] = Pear;   //X-Position          (units: au)
		y[9] = 0;                      //X-Velocity (units: au/day)
		y[10] = 0;             //Y- Position         (units: au)
		y[11] = Vear; 		   //Y-Velocity (units: au/day)

		//Mars
		y[12] = Pmar;   //X-Position          (units: au)
		y[13] = 0;                      //X-Velocity (units: au/day)
		y[14] = 0;             //Y- Position         (units: au)
		y[15] = Vmar; 		   //Y-Velocity (units: au/day)

		//Jupiter
		y[16] = Pjup;   //X-Position          (units: au)
		y[17] = 0;                      //X-Velocity (units: au/day)
		y[18] = 0;             //Y- Position         (units: au)
		y[19] = Vjup; 		   //Y-Velocity (units: au/day)

		//Saturn
		y[20] = Psat;   //X-Position          (units: au)
		y[21] = 0;                      //X-Velocity (units: au/day)
		y[22] = 0;             //Y- Position         (units: au)
		y[23] = Vsat; 		   //Y-Velocity (units: au/day)

		//Uranus
		y[24] = Pura;   //X-Position          (units: au)
		y[25] = 0;                      //X-Velocity (units: au/day)
		y[26] = 0;             //Y- Position         (units: au)
		y[27] = Vura; 		   //Y-Velocity (units: au/day)

		//Pluto
		y[28] = Pplu;   //X-Position          (units: au)
		y[29] = 0;                      //X-Velocity (units: au/day)
		y[30] = 0;             //Y- Position         (units: au)
		y[31] = Vplu; 		   //Y-Velocity (units: au/day)

		//Neptune
		y[32] = Pnep;   //X-Position          (units: au)
		y[33] = 0;                      //X-Velocity (units: au/day)
		y[34] = 0;             //Y- Position         (units: au)
		y[35] = Vnep; 		   //Y-Velocity (units: au/day)


		GRmercDerivs(xmin, y, dydx);

			//file output stream
			ofstream ofPositionMercury;
			ofPositionMercury.open("position.GR_Mercury.csv");



		for(int k=0; k < kmax; k++)
			{
					x=xmin+k*h;

					rk4(y, dydx,  x, h, yout, GRmercDerivs);


					///debug output
					//cout << "k = " << k
					//	 << "    Xm = " << yout[0] << "    X'm= " << yout[1]
					//	 << "    Ym = " << yout[2] << "    Y'm= " << yout[3] << endl;


					//file output stream
					ofPositionMercury << yout[0] << "," << yout[2] << "\n";

					y = yout;

					GRmercDerivs(x,y,dydx);
			}

			//closes file output stream
			ofPositionMercury.close();


		cout << "\n" << "CoMpLeTe CoMpLeTe CoMpLeTe CoMpLeTe CoMpLeTe CoMpLeTe" << "\n";

	return 0;
}
*/



//--------/////Non-relativistic orbit in Solar System/////---------

/*
void mercDerivs(const Doub x, VecDoub_I & y, VecDoub_O & dydx)

{
	//exponents for radius calculation in R(i,j,exp)
	double exp1 = 1.5;


///Mercury------------------------------------------------------------------------------------
int i = 0;
int j = 2;

//X-components
dydx[0] = y[1];
dydx[1] = (-GMsun * y[i]) / Re(y[i]       , y[j]       , exp1)

		 + (-GMven * y[i]) / Re(y[i]-y[4]  , y[j]-y[6]  , exp1)
		 + (-GMear * y[i]) / Re(y[i]-y[8]  , y[j]-y[10] , exp1)
		 + (-GMmar * y[i]) / Re(y[i]-y[12] , y[j]-y[14] , exp1)
		 + (-GMjup * y[i]) / Re(y[i]-y[16] , y[j]-y[18] , exp1)
		 + (-GMsat * y[i]) / Re(y[i]-y[20] , y[j]-y[22] , exp1)
		 + (-GMura * y[i]) / Re(y[i]-y[24] , y[j]-y[26] , exp1)
		 + (-GMplu * y[i]) / Re(y[i]-y[28] , y[j]-y[30] , exp1)
		 + (-GMnep * y[i]) / Re(y[i]-y[32] , y[j]-y[34] , exp1);

//Y-components
dydx[2] = y[3];
dydx[3] = (-GMsun * y[j]) / Re(y[i]       , y[j]       , exp1)

		 + (-GMven * y[j]) / Re(y[i]-y[4]  , y[j]-y[6]  , exp1)
		 + (-GMear * y[j]) / Re(y[i]-y[8]  , y[j]-y[10] , exp1)
		 + (-GMmar * y[j]) / Re(y[i]-y[12] , y[j]-y[14] , exp1)
		 + (-GMjup * y[j]) / Re(y[i]-y[16] , y[j]-y[18] , exp1)
		 + (-GMsat * y[j]) / Re(y[i]-y[20] , y[j]-y[22] , exp1)
		 + (-GMura * y[j]) / Re(y[i]-y[24] , y[j]-y[26] , exp1)
		 + (-GMplu * y[j]) / Re(y[i]-y[28] , y[j]-y[30] , exp1)
		 + (-GMnep * y[j]) / Re(y[i]-y[32] , y[j]-y[34] , exp1);

///Venus--------------------------------------------------------------------------------------
i = 4;
j = 6;

//X-components
dydx[4] = y[5];
dydx[5] = (-GMsun * y[i]) / Re(y[i]       , y[j]       , exp1)
		 + (-GMmer * y[i]) / Re(y[i]-y[0]  , y[j]-y[2]  , exp1)

		 + (-GMear * y[i]) / Re(y[i]-y[8]  , y[j]-y[10] , exp1)
		 + (-GMmar * y[i]) / Re(y[i]-y[12] , y[j]-y[14] , exp1)
		 + (-GMjup * y[i]) / Re(y[i]-y[16] , y[j]-y[18] , exp1)
		 + (-GMsat * y[i]) / Re(y[i]-y[20] , y[j]-y[22] , exp1)
		 + (-GMura * y[i]) / Re(y[i]-y[24] , y[j]-y[26] , exp1)
		 + (-GMplu * y[i]) / Re(y[i]-y[28] , y[j]-y[30] , exp1)
		 + (-GMnep * y[i]) / Re(y[i]-y[32] , y[j]-y[34] , exp1);

//Y-components
dydx[6] = y[7];
dydx[7] = (-GMsun * y[j]) / Re(y[i]       , y[j]       , exp1)
		 + (-GMmer * y[j]) / Re(y[i]-y[0]  , y[j]-y[2]  , exp1)

		 + (-GMear * y[j]) / Re(y[i]-y[8]  , y[j]-y[10] , exp1)
		 + (-GMmar * y[j]) / Re(y[i]-y[12] , y[j]-y[14] , exp1)
		 + (-GMjup * y[j]) / Re(y[i]-y[16] , y[j]-y[18] , exp1)
		 + (-GMsat * y[j]) / Re(y[i]-y[20] , y[j]-y[22] , exp1)
		 + (-GMura * y[j]) / Re(y[i]-y[24] , y[j]-y[26] , exp1)
		 + (-GMplu * y[j]) / Re(y[i]-y[28] , y[j]-y[30] , exp1)
		 + (-GMnep * y[j]) / Re(y[i]-y[32] , y[j]-y[34] , exp1);

///Earth---------------------------------------------------------------------------------------
i = 8;
j = 10;

//X-components
dydx[8] = y[9];
dydx[9] = (-GMsun * y[i]) / Re(y[i]       , y[j]       , exp1)
		 + (-GMmer * y[i]) / Re(y[i]-y[0]  , y[j]-y[2]  , exp1)
		 + (-GMven * y[i]) / Re(y[i]-y[4]  , y[j]-y[6]  , exp1)

		 + (-GMmar * y[i]) / Re(y[i]-y[12] , y[j]-y[14] , exp1)
		 + (-GMjup * y[i]) / Re(y[i]-y[16] , y[j]-y[18] , exp1)
		 + (-GMsat * y[i]) / Re(y[i]-y[20] , y[j]-y[22] , exp1)
		 + (-GMura * y[i]) / Re(y[i]-y[24] , y[j]-y[26] , exp1)
		 + (-GMplu * y[i]) / Re(y[i]-y[28] , y[j]-y[30] , exp1)
		 + (-GMnep * y[i]) / Re(y[i]-y[32] , y[j]-y[34] , exp1);

//Y-components
dydx[10] = y[11];
dydx[11] = (-GMsun * y[j]) / Re(y[i]       , y[j]       , exp1)
		 + (-GMmer * y[j]) / Re(y[i]-y[0]  , y[j]-y[2]  , exp1)
		 + (-GMven * y[j]) / Re(y[i]-y[4]  , y[j]-y[6]  , exp1)

		 + (-GMmar * y[j]) / Re(y[i]-y[12] , y[j]-y[14] , exp1)
		 + (-GMjup * y[j]) / Re(y[i]-y[16] , y[j]-y[18] , exp1)
		 + (-GMsat * y[j]) / Re(y[i]-y[20] , y[j]-y[22] , exp1)
		 + (-GMura * y[j]) / Re(y[i]-y[24] , y[j]-y[26] , exp1)
		 + (-GMplu * y[j]) / Re(y[i]-y[28] , y[j]-y[30] , exp1)
		 + (-GMnep * y[j]) / Re(y[i]-y[32] , y[j]-y[34] , exp1);

///Mars-----------------------------------------------------------------------------------------
i = 12;
j = 14;

//X-components
dydx[12] = y[13];
dydx[13] = (-GMsun * y[i]) / Re(y[i]       , y[j]       , exp1)
		 + (-GMmer * y[i]) / Re(y[i]-y[0]  , y[j]-y[2]  , exp1)
		 + (-GMven * y[i]) / Re(y[i]-y[4]  , y[j]-y[6]  , exp1)
		 + (-GMear * y[i]) / Re(y[i]-y[8]  , y[j]-y[10] , exp1)

		 + (-GMjup * y[i]) / Re(y[i]-y[16] , y[j]-y[18] , exp1)
		 + (-GMsat * y[i]) / Re(y[i]-y[20] , y[j]-y[22] , exp1)
		 + (-GMura * y[i]) / Re(y[i]-y[24] , y[j]-y[26] , exp1)
		 + (-GMplu * y[i]) / Re(y[i]-y[28] , y[j]-y[30] , exp1)
		 + (-GMnep * y[i]) / Re(y[i]-y[32] , y[j]-y[34] , exp1);

//Y-components
dydx[14] = y[15];
dydx[15] = (-GMsun * y[j]) / Re(y[i]       , y[j]       , exp1)
		 + (-GMmer * y[j]) / Re(y[i]-y[0]  , y[j]-y[2]  , exp1)
		 + (-GMven * y[j]) / Re(y[i]-y[4]  , y[j]-y[6]  , exp1)
		 + (-GMear * y[j]) / Re(y[i]-y[8]  , y[j]-y[10] , exp1)

		 + (-GMjup * y[j]) / Re(y[i]-y[16] , y[j]-y[18] , exp1)
		 + (-GMsat * y[j]) / Re(y[i]-y[20] , y[j]-y[22] , exp1)
		 + (-GMura * y[j]) / Re(y[i]-y[24] , y[j]-y[26] , exp1)
		 + (-GMplu * y[j]) / Re(y[i]-y[28] , y[j]-y[30] , exp1)
		 + (-GMnep * y[j]) / Re(y[i]-y[32] , y[j]-y[34] , exp1);

///Jupiter---------------------------------------------------------------------------------------
i = 16;
j = 18;

//X-components
dydx[16] = y[17];
dydx[17] = (-GMsun * y[i]) / Re(y[i]       , y[j]       , exp1)
		 + (-GMmer * y[i]) / Re(y[i]-y[0]  , y[j]-y[2]  , exp1)
		 + (-GMven * y[i]) / Re(y[i]-y[4]  , y[j]-y[6]  , exp1)
		 + (-GMear * y[i]) / Re(y[i]-y[8]  , y[j]-y[10] , exp1)
		 + (-GMmar * y[i]) / Re(y[i]-y[12] , y[j]-y[14] , exp1)
		 //jupiter removed
		 + (-GMsat * y[i]) / Re(y[i]-y[20] , y[j]-y[22] , exp1)
		 + (-GMura * y[i]) / Re(y[i]-y[24] , y[j]-y[26] , exp1)
		 + (-GMplu * y[i]) / Re(y[i]-y[28] , y[j]-y[30] , exp1)
		 + (-GMnep * y[i]) / Re(y[i]-y[32] , y[j]-y[34] , exp1);

//Y-components
dydx[18] = y[19];
dydx[19] = (-GMsun * y[j]) / Re(y[i]       , y[j]       , exp1)
	   	 + (-GMmer * y[j]) / Re(y[i]-y[0]  , y[j]-y[2]  , exp1)
		 + (-GMven * y[j]) / Re(y[i]-y[4]  , y[j]-y[6]  , exp1)
		 + (-GMear * y[j]) / Re(y[i]-y[8]  , y[j]-y[10] , exp1)
		 + (-GMmar * y[j]) / Re(y[i]-y[12] , y[j]-y[14] , exp1)
		//jupiter removed
		 + (-GMsat * y[j]) / Re(y[i]-y[20] , y[j]-y[22] , exp1)
		 + (-GMura * y[j]) / Re(y[i]-y[24] , y[j]-y[26] , exp1)
		 + (-GMplu * y[j]) / Re(y[i]-y[28] , y[j]-y[30] , exp1)
		 + (-GMnep * y[j]) / Re(y[i]-y[32] , y[j]-y[34] , exp1);

///Saturn---------------------------------------------------------------------------------------
i = 20;
j = 22;

//X-components
dydx[20] = y[21];
dydx[21] = (-GMsun * y[i]) / Re(y[i]       , y[j]       , exp1)
		 + (-GMmer * y[i]) / Re(y[i]-y[0]  , y[j]-y[2]  , exp1)
		 + (-GMven * y[i]) / Re(y[i]-y[4]  , y[j]-y[6]  , exp1)
		 + (-GMear * y[i]) / Re(y[i]-y[8]  , y[j]-y[10] , exp1)
		 + (-GMmar * y[i]) / Re(y[i]-y[12] , y[j]-y[14] , exp1)
		 + (-GMjup * y[i]) / Re(y[i]-y[16] , y[j]-y[18] , exp1)
		 //saturn removed
		 + (-GMura * y[i]) / Re(y[i]-y[24] , y[j]-y[26] , exp1)
		 + (-GMplu * y[i]) / Re(y[i]-y[28] , y[j]-y[30] , exp1)
		 + (-GMnep * y[i]) / Re(y[i]-y[32] , y[j]-y[34] , exp1);

//Y-components
dydx[22] = y[23];
dydx[23] = (-GMsun * y[j]) / Re(y[i]       , y[j]       , exp1)
		 + (-GMmer * y[j]) / Re(y[i]-y[0]  , y[j]-y[2]  , exp1)
		 + (-GMven * y[j]) / Re(y[i]-y[4]  , y[j]-y[6]  , exp1)
		 + (-GMear * y[j]) / Re(y[i]-y[8]  , y[j]-y[10] , exp1)
		 + (-GMmar * y[j]) / Re(y[i]-y[12] , y[j]-y[14] , exp1)
		 + (-GMjup * y[j]) / Re(y[i]-y[16] , y[j]-y[18] , exp1)
		//saturn removed
		 + (-GMura * y[j]) / Re(y[i]-y[24] , y[j]-y[26] , exp1)
		 + (-GMplu * y[j]) / Re(y[i]-y[28] , y[j]-y[30] , exp1)
		 + (-GMnep * y[j]) / Re(y[i]-y[32] , y[j]-y[34] , exp1);

///Uranus---------------------------------------------------------------------------------------
i = 24;
j = 26;

//X-components
dydx[24] = y[25];
dydx[25] = (-GMsun * y[i]) / Re(y[i]       , y[j]       , exp1)
		 + (-GMmer * y[i]) / Re(y[i]-y[0]  , y[j]-y[2]  , exp1)
		 + (-GMven * y[i]) / Re(y[i]-y[4]  , y[j]-y[6]  , exp1)
		 + (-GMear * y[i]) / Re(y[i]-y[8]  , y[j]-y[10] , exp1)
		 + (-GMmar * y[i]) / Re(y[i]-y[12] , y[j]-y[14] , exp1)
		 + (-GMjup * y[i]) / Re(y[i]-y[16] , y[j]-y[18] , exp1)
		 + (-GMsat * y[i]) / Re(y[i]-y[20] , y[j]-y[22] , exp1)
		 //uranus removed
		 + (-GMplu * y[i]) / Re(y[i]-y[28] , y[j]-y[30] , exp1)
		 + (-GMnep * y[i]) / Re(y[i]-y[32] , y[j]-y[34] , exp1);

//Y-components
dydx[26] = y[27];
dydx[27] = (-GMsun * y[j]) / Re(y[i]       , y[j]       , exp1)
		 + (-GMmer * y[j]) / Re(y[i]-y[0]  , y[j]-y[2]  , exp1)
		 + (-GMven * y[j]) / Re(y[i]-y[4]  , y[j]-y[6]  , exp1)
		 + (-GMear * y[j]) / Re(y[i]-y[8]  , y[j]-y[10] , exp1)
		 + (-GMmar * y[j]) / Re(y[i]-y[12] , y[j]-y[14] , exp1)
		 + (-GMjup * y[j]) / Re(y[i]-y[16] , y[j]-y[18] , exp1)
		 + (-GMsat * y[j]) / Re(y[i]-y[20] , y[j]-y[22] , exp1)
		 //uranus removed
		 + (-GMplu * y[j]) / Re(y[i]-y[28] , y[j]-y[30] , exp1)
		 + (-GMnep * y[j]) / Re(y[i]-y[32] , y[j]-y[34] , exp1);

///Pluto---------------------------------------------------------------------------------------
i = 28;
j = 30;

//X-components
dydx[28] = y[29];
dydx[29] = (-GMsun * y[i]) / Re(y[i]       , y[j]       , exp1)
		 + (-GMmer * y[i]) / Re(y[i]-y[0]  , y[j]-y[2]  , exp1)
		 + (-GMven * y[i]) / Re(y[i]-y[4]  , y[j]-y[6]  , exp1)
		 + (-GMear * y[i]) / Re(y[i]-y[8]  , y[j]-y[10] , exp1)
		 + (-GMmar * y[i]) / Re(y[i]-y[12] , y[j]-y[14] , exp1)
		 + (-GMjup * y[i]) / Re(y[i]-y[16] , y[j]-y[18] , exp1)
		 + (-GMsat * y[i]) / Re(y[i]-y[20] , y[j]-y[22] , exp1)
		 + (-GMura * y[i]) / Re(y[i]-y[24] , y[j]-y[26] , exp1)
		 //pluto removed
		 + (-GMnep * y[i]) / Re(y[i]-y[32] , y[j]-y[34] , exp1);

//Y-components
dydx[30] = y[31];
dydx[31] = (-GMsun * y[j]) / Re(y[i]       , y[j]       , exp1)
		 + (-GMmer * y[j]) / Re(y[i]-y[0]  , y[j]-y[2]  , exp1)
		 + (-GMven * y[j]) / Re(y[i]-y[4]  , y[j]-y[6]  , exp1)
		 + (-GMear * y[j]) / Re(y[i]-y[8]  , y[j]-y[10] , exp1)
		 + (-GMmar * y[j]) / Re(y[i]-y[12] , y[j]-y[14] , exp1)
		 + (-GMjup * y[j]) / Re(y[i]-y[16] , y[j]-y[18] , exp1)
		 + (-GMsat * y[j]) / Re(y[i]-y[20] , y[j]-y[22] , exp1)
		 + (-GMura * y[j]) / Re(y[i]-y[24] , y[j]-y[26] , exp1)
		 //pluto removed
		 + (-GMnep * y[j]) / Re(y[i]-y[32] , y[j]-y[34] , exp1);

///Neptune---------------------------------------------------------------------------------------
i = 32;
j = 34;

//X-components
dydx[32] = y[33];
dydx[33] = (-GMsun * y[i]) / Re(y[i]       , y[j]       , exp1)
		 + (-GMmer * y[i]) / Re(y[i]-y[0]  , y[j]-y[2]  , exp1)
		 + (-GMven * y[i]) / Re(y[i]-y[4]  , y[j]-y[6]  , exp1)
		 + (-GMear * y[i]) / Re(y[i]-y[8]  , y[j]-y[10] , exp1)
		 + (-GMmar * y[i]) / Re(y[i]-y[12] , y[j]-y[14] , exp1)
		 + (-GMjup * y[i]) / Re(y[i]-y[16] , y[j]-y[18] , exp1)
		 + (-GMsat * y[i]) / Re(y[i]-y[20] , y[j]-y[22] , exp1)
		 + (-GMura * y[i]) / Re(y[i]-y[24] , y[j]-y[26] , exp1)
		 + (-GMplu * y[i]) / Re(y[i]-y[28] , y[j]-y[30] , exp1);
		//neptune removed

//Y-components
dydx[34] = y[35];
dydx[35] = (-GMsun * y[j]) / Re(y[i]       , y[j]       , exp1)
		 + (-GMmer * y[j]) / Re(y[i]-y[0]  , y[j]-y[2]  , exp1)
		 + (-GMven * y[j]) / Re(y[i]-y[4]  , y[j]-y[6]  , exp1)
		 + (-GMear * y[j]) / Re(y[i]-y[8]  , y[j]-y[10] , exp1)
		 + (-GMmar * y[j]) / Re(y[i]-y[12] , y[j]-y[14] , exp1)
		 + (-GMjup * y[j]) / Re(y[i]-y[16] , y[j]-y[18] , exp1)
		 + (-GMsat * y[j]) / Re(y[i]-y[20] , y[j]-y[22] , exp1)
		 + (-GMura * y[j]) / Re(y[i]-y[24] , y[j]-y[26] , exp1)
		 + (-GMplu * y[j]) / Re(y[i]-y[28] , y[j]-y[30] , exp1);
		//neptune removed

}
*/

/*
int mercury()
{
		cout << "Started Mercury -- Started Mercury -- Started Mercury" << "\n";



		VecDoub y(36);
		VecDoub dydx(36); //vector of positions & velocities for earth and
		VecDoub yout(36);


		double x;
		double xmin = 0;       //minimum starting position (units: au)
		double kmax = 10000*87.969;  //max iterations (units: days)
		double h = .5;          //time step size (units: days)


		///Initial Conditions:

		//Mercury
		y[0] = Pmer;   //X-Position          (units: au)
		y[1] = 0;                      //X-Velocity (units: au/day)
		y[2] = 0;             //Y- Position         (units: au)
		y[3] = Vmer; 		   //Y-Velocity (units: au/day)

		//Venus
		y[4] = Pven;   //X-Position          (units: au)
		y[5] = 0;                      //X-Velocity (units: au/day)
		y[6] = 0;             //Y- Position         (units: au)
		y[7] = Vven; 		   //Y-Velocity (units: au/day)

		//Earth
		y[8] = Pear;   //X-Position          (units: au)
		y[9] = 0;                      //X-Velocity (units: au/day)
		y[10] = 0;             //Y- Position         (units: au)
		y[11] = Vear; 		   //Y-Velocity (units: au/day)

		//Mars
		y[12] = Pmar;   //X-Position          (units: au)
		y[13] = 0;                      //X-Velocity (units: au/day)
		y[14] = 0;             //Y- Position         (units: au)
		y[15] = Vmar; 		   //Y-Velocity (units: au/day)

		//Jupiter
		y[16] = Pjup;   //X-Position          (units: au)
		y[17] = 0;                      //X-Velocity (units: au/day)
		y[18] = 0;             //Y- Position         (units: au)
		y[19] = Vjup; 		   //Y-Velocity (units: au/day)

		//Saturn
		y[20] = Psat;   //X-Position          (units: au)
		y[21] = 0;                      //X-Velocity (units: au/day)
		y[22] = 0;             //Y- Position         (units: au)
		y[23] = Vsat; 		   //Y-Velocity (units: au/day)

		//Uranus
		y[24] = Pura;   //X-Position          (units: au)
		y[25] = 0;                      //X-Velocity (units: au/day)
		y[26] = 0;             //Y- Position         (units: au)
		y[27] = Vura; 		   //Y-Velocity (units: au/day)

		//Pluto
		y[28] = Pplu;   //X-Position          (units: au)
		y[29] = 0;                      //X-Velocity (units: au/day)
		y[30] = 0;             //Y- Position         (units: au)
		y[31] = Vplu; 		   //Y-Velocity (units: au/day)

		//Neptune
		y[32] = Pnep;   //X-Position          (units: au)
		y[33] = 0;                      //X-Velocity (units: au/day)
		y[34] = 0;             //Y- Position         (units: au)
		y[35] = Vnep; 		   //Y-Velocity (units: au/day)


		mercDerivs(xmin, y, dydx);

			//file output stream
			ofstream ofPositionMercury, ofPerihelion;
			ofPositionMercury.open("x-y-positionMercury.csv");
			ofPerihelion.open("perihleionMercury.csv");


		for(int k=0; k < kmax; k++)
			{
					x=xmin+k*h;

					rk4(y, dydx,  x, h, yout, mercDerivs);


					//display output
					//cout << "k = " << k
					//	 << "    Xm = " << yout[0] << "    X'm= " << yout[1]
					//	 << "    Ym = " << yout[2] << "    Y'm= " << yout[3] << endl;


					//file output stream
					ofPositionMercury << yout[0] << "," << yout[2] << "\n";



					//Perihelion Tracker
					//P-mercury = 0.307491008
					double a = 0.307491000;
					double b = 0.30749101;
						if(y[0]>a && y[0]<b)
						{
							cout << "k=" <<k<< "  Perihelion @ (" << y[0] << ", " << y[2] << ")  y[0]=" << y[0] << endl;
							ofPerihelion << yout[0] << "," << yout[2] << "\n";
						}


						if(y[2]>a && y[2]<b)
						{
							cout << "k=" <<k<< "  Perihelion @ (" << y[0] << ", " << y[2] << ")  y[2]=" << y[2] << endl;
							ofPerihelion << yout[0] << "," << yout[2] << "\n";
						}


						if(pow((pow(y[0],2)+pow(y[2],2)),0.5)>a && pow((pow(y[0],2)+pow(y[2],2)),0.5)<b)
						{
							cout << "k=" <<k<< "  Perihelion @ (" << y[0] << ", " << y[2] << ")  pow(...,2)=" << pow((pow(y[0],2)+pow(y[2],2)),0.5) << endl;
							ofPerihelion << yout[0] << "," << yout[2] << "\n";
						}



					y = yout;

					mercDerivs(x,y,dydx);
			}

			//closes file output stream
			ofPositionMercury.close();
			ofPerihelion.close();

		cout << "\n" << "CoMpLeTe CoMpLeTe CoMpLeTe CoMpLeTe CoMpLeTe CoMpLeTe" << "\n";

	return 0;
}
*/

////////////////////THE WHOLE SOLAR SYSTEM////////////////////////














//working sun and mercury. no deviations in plot.
/*
 int sm()
{
		VecDoub y(4);
		VecDoub dydx(4); //vector of positions & velocities for earth and
		VecDoub yout(4);


		double x;
		double xmin = 0;       //minimum starting position (units: au)
		double kmax = 50*87.969;  //max iterations (units: days)
		double h = 1;          //time step size (units: days)


		///Initial Conditions:

		//Mercury
		y[0] = 0.307491008;   //X-Position          (units: au)
		y[1] = 0;                      //X-Velocity (units: au/day)
		y[2] = 0;             //Y- Position         (units: au)
		y[3] = 0.0340638003; 		   //Y-Velocity (units: au/day)

		mercDerivs(xmin, y, dydx);

			//file output stream
			ofstream ofPositionMercury;
			ofPositionMercury.open("x-y-positionMercury.csv");


		for(int k=0; k < kmax; k++)
			{
					x=xmin+k*h;

					rk4(y, dydx,  x, h, yout, mercDerivs);

					//display output
					cout << "k = " << k
						 << "    Xm = " << yout[0] << "    X'm= " << yout[1]
						 << "    Ym = " << yout[2] << "    Y'm= " << yout[3] << endl;

					//file output stream
					ofPositionMercury << yout[0] << "," << yout[2] << endl;

					y = yout;

					mercDerivs(x,y,dydx);
			}

			//closes file output stream
			ofPositionMercury.close();

	return 0;
}
 */

/*
//Last Working w/ N.M.G. example.

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


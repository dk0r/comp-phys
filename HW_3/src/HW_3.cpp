
#include <iostream>
#include <fstream>
#include <cmath>
#include <nr3.h>
#include <rk4.h>

using namespace std;


//////Physical Constants////

///Gravitational Constants
double Gsi  = 6.67398e-11; // (SI-units:  m^3 / kg*s^2)              [old# used: 6.67428e-11]
double Gastro = 1.48811382e-34; // (Astro-units: au^3 / kg*days^2)   [old# used: 1.48818071e-34]

//Celestial Mass  (SI-units: kg)
double Msun =    1.9891e30;
double Mmer =    0.3302e24;
double Mven =    4.8685e24;
double Mear =    5.9736e24;
double Mmar =    0.64185e24;
double Mjup = 1898.6e24;
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
double GMsun = 2.96014025e-4;
double GMmer = Gastro*Mmer;
double GMven = Gastro*Mven;
double GMear = 8.88979629e-10;
double GMmar = Gastro*Mmar;
double GMjup = 2.82545990e-7; //2.82545990e-7;
double GMsat = Gastro*Msat;
double GMura = Gastro*Mura;
double GMplu = Gastro*Mplu;
double GMnep = Gastro*Mnep;



//Calculates radius raised to some power (exp)
double R(double x, double y, double exp)
{
	double r;

	r = pow( pow(x,2)+pow(y,2) , exp );

	return r;
}






void mercDerivs(const Doub x, VecDoub_I & y, VecDoub_O & dydx)

{
	//exponents for radius calculation in R(i,j,exp)
	double exp1 = 1.5;
	double exp2 = 2.5;

///Mercury------------------------------------------------------------------------------------
int i = 0;
int j = 2;

//X-components
dydx[0] = y[1];
dydx[1] = (-GMsun * y[i]) / R(y[i]       , y[j]       , exp1)

		 + (-GMven * y[i]) / R(y[i]-y[4]  , y[j]-y[6]  , exp1)
		 + (-GMear * y[i]) / R(y[i]-y[8]  , y[j]-y[10] , exp1)
		 + (-GMmar * y[i]) / R(y[i]-y[12] , y[j]-y[14] , exp1)
		 + (-GMjup * y[i]) / R(y[i]-y[16] , y[j]-y[18] , exp1)
		 + (-GMsat * y[i]) / R(y[i]-y[20] , y[j]-y[22] , exp1)
		 + (-GMura * y[i]) / R(y[i]-y[24] , y[j]-y[26] , exp1)
		 + (-GMplu * y[i]) / R(y[i]-y[28] , y[j]-y[30] , exp1)
		 + (-GMnep * y[i]) / R(y[i]-y[32] , y[j]-y[34] , exp1);

//Y-components
dydx[2] = y[3];
dydx[3] = (-GMsun * y[j]) / R(y[i]       , y[j]       , exp1)

		 + (-GMven * y[j]) / R(y[i]-y[4]  , y[j]-y[6]  , exp1)
		 + (-GMear * y[j]) / R(y[i]-y[8]  , y[j]-y[10] , exp1)
		 + (-GMmar * y[j]) / R(y[i]-y[12] , y[j]-y[14] , exp1)
		 + (-GMjup * y[j]) / R(y[i]-y[16] , y[j]-y[18] , exp1)
		 + (-GMsat * y[j]) / R(y[i]-y[20] , y[j]-y[22] , exp1)
		 + (-GMura * y[j]) / R(y[i]-y[24] , y[j]-y[26] , exp1)
		 + (-GMplu * y[j]) / R(y[i]-y[28] , y[j]-y[30] , exp1)
		 + (-GMnep * y[j]) / R(y[i]-y[32] , y[j]-y[34] , exp1);

///Venus--------------------------------------------------------------------------------------
i = 4;
j = 6;

//X-components
dydx[4] = y[5];
dydx[5] = (-GMsun * y[i]) / R(y[i]       , y[j]       , exp1)
		 + (-GMmer * y[i]) / R(y[i]-y[0]  , y[j]-y[2]  , exp1)

		 + (-GMear * y[i]) / R(y[i]-y[8]  , y[j]-y[10] , exp1)
		 + (-GMmar * y[i]) / R(y[i]-y[12] , y[j]-y[14] , exp1)
		 + (-GMjup * y[i]) / R(y[i]-y[16] , y[j]-y[18] , exp1)
		 + (-GMsat * y[i]) / R(y[i]-y[20] , y[j]-y[22] , exp1)
		 + (-GMura * y[i]) / R(y[i]-y[24] , y[j]-y[26] , exp1)
		 + (-GMplu * y[i]) / R(y[i]-y[28] , y[j]-y[30] , exp1)
		 + (-GMnep * y[i]) / R(y[i]-y[32] , y[j]-y[34] , exp1);

//Y-components
dydx[6] = y[7];
dydx[7] = (-GMsun * y[j]) / R(y[i]       , y[j]       , exp1)
		 + (-GMmer * y[j]) / R(y[i]-y[0]  , y[j]-y[2]  , exp1)

		 + (-GMear * y[j]) / R(y[i]-y[8]  , y[j]-y[10] , exp1)
		 + (-GMmar * y[j]) / R(y[i]-y[12] , y[j]-y[14] , exp1)
		 + (-GMjup * y[j]) / R(y[i]-y[16] , y[j]-y[18] , exp1)
		 + (-GMsat * y[j]) / R(y[i]-y[20] , y[j]-y[22] , exp1)
		 + (-GMura * y[j]) / R(y[i]-y[24] , y[j]-y[26] , exp1)
		 + (-GMplu * y[j]) / R(y[i]-y[28] , y[j]-y[30] , exp1)
		 + (-GMnep * y[j]) / R(y[i]-y[32] , y[j]-y[34] , exp1);

///Earth---------------------------------------------------------------------------------------
i = 8;
j = 10;

//X-components
dydx[8] = y[9];
dydx[9] = (-GMsun * y[i]) / R(y[i]       , y[j]       , exp1)
		 + (-GMmer * y[i]) / R(y[i]-y[0]  , y[j]-y[2]  , exp1)
		 + (-GMven * y[i]) / R(y[i]-y[4]  , y[j]-y[6]  , exp1)

		 + (-GMmar * y[i]) / R(y[i]-y[12] , y[j]-y[14] , exp1)
		 + (-GMjup * y[i]) / R(y[i]-y[16] , y[j]-y[18] , exp1)
		 + (-GMsat * y[i]) / R(y[i]-y[20] , y[j]-y[22] , exp1)
		 + (-GMura * y[i]) / R(y[i]-y[24] , y[j]-y[26] , exp1)
		 + (-GMplu * y[i]) / R(y[i]-y[28] , y[j]-y[30] , exp1)
		 + (-GMnep * y[i]) / R(y[i]-y[32] , y[j]-y[34] , exp1);

//Y-components
dydx[10] = y[11];
dydx[11] = (-GMsun * y[j]) / R(y[i]       , y[j]       , exp1)
		 + (-GMmer * y[j]) / R(y[i]-y[0]  , y[j]-y[2]  , exp1)
		 + (-GMven * y[j]) / R(y[i]-y[4]  , y[j]-y[6]  , exp1)

		 + (-GMmar * y[j]) / R(y[i]-y[12] , y[j]-y[14] , exp1)
		 + (-GMjup * y[j]) / R(y[i]-y[16] , y[j]-y[18] , exp1)
		 + (-GMsat * y[j]) / R(y[i]-y[20] , y[j]-y[22] , exp1)
		 + (-GMura * y[j]) / R(y[i]-y[24] , y[j]-y[26] , exp1)
		 + (-GMplu * y[j]) / R(y[i]-y[28] , y[j]-y[30] , exp1)
		 + (-GMnep * y[j]) / R(y[i]-y[32] , y[j]-y[34] , exp1);

///Mars-----------------------------------------------------------------------------------------
i = 12;
j = 14;

//X-components
dydx[12] = y[13];
dydx[13] = (-GMsun * y[i]) / R(y[i]       , y[j]       , exp1)
		 + (-GMmer * y[i]) / R(y[i]-y[0]  , y[j]-y[2]  , exp1)
		 + (-GMven * y[i]) / R(y[i]-y[4]  , y[j]-y[6]  , exp1)
		 + (-GMear * y[i]) / R(y[i]-y[8]  , y[j]-y[10] , exp1)

		 + (-GMjup * y[i]) / R(y[i]-y[16] , y[j]-y[18] , exp1)
		 + (-GMsat * y[i]) / R(y[i]-y[20] , y[j]-y[22] , exp1)
		 + (-GMura * y[i]) / R(y[i]-y[24] , y[j]-y[26] , exp1)
		 + (-GMplu * y[i]) / R(y[i]-y[28] , y[j]-y[30] , exp1)
		 + (-GMnep * y[i]) / R(y[i]-y[32] , y[j]-y[34] , exp1);

//Y-components
dydx[14] = y[15];
dydx[15] = (-GMsun * y[j]) / R(y[i]       , y[j]       , exp1)
		 + (-GMmer * y[j]) / R(y[i]-y[0]  , y[j]-y[2]  , exp1)
		 + (-GMven * y[j]) / R(y[i]-y[4]  , y[j]-y[6]  , exp1)
		 + (-GMear * y[j]) / R(y[i]-y[8]  , y[j]-y[10] , exp1)

		 + (-GMjup * y[j]) / R(y[i]-y[16] , y[j]-y[18] , exp1)
		 + (-GMsat * y[j]) / R(y[i]-y[20] , y[j]-y[22] , exp1)
		 + (-GMura * y[j]) / R(y[i]-y[24] , y[j]-y[26] , exp1)
		 + (-GMplu * y[j]) / R(y[i]-y[28] , y[j]-y[30] , exp1)
		 + (-GMnep * y[j]) / R(y[i]-y[32] , y[j]-y[34] , exp1);

///Jupiter---------------------------------------------------------------------------------------
i = 16;
j = 18;

//X-components
dydx[16] = y[17];
dydx[17] = (-GMsun * y[i]) / R(y[i]       , y[j]       , exp1)
		 + (-GMmer * y[i]) / R(y[i]-y[0]  , y[j]-y[2]  , exp1)
		 + (-GMven * y[i]) / R(y[i]-y[4]  , y[j]-y[6]  , exp1)
		 + (-GMear * y[i]) / R(y[i]-y[8]  , y[j]-y[10] , exp1)
		 + (-GMmar * y[i]) / R(y[i]-y[12] , y[j]-y[14] , exp1)
		 //jupiter removed
		 + (-GMsat * y[i]) / R(y[i]-y[20] , y[j]-y[22] , exp1)
		 + (-GMura * y[i]) / R(y[i]-y[24] , y[j]-y[26] , exp1)
		 + (-GMplu * y[i]) / R(y[i]-y[28] , y[j]-y[30] , exp1)
		 + (-GMnep * y[i]) / R(y[i]-y[32] , y[j]-y[34] , exp1);

//Y-components
dydx[18] = y[19];
dydx[19] = (-GMsun * y[j]) / R(y[i]       , y[j]       , exp1)
	   	 + (-GMmer * y[j]) / R(y[i]-y[0]  , y[j]-y[2]  , exp1)
		 + (-GMven * y[j]) / R(y[i]-y[4]  , y[j]-y[6]  , exp1)
		 + (-GMear * y[j]) / R(y[i]-y[8]  , y[j]-y[10] , exp1)
		 + (-GMmar * y[j]) / R(y[i]-y[12] , y[j]-y[14] , exp1)
		//jupiter removed
		 + (-GMsat * y[j]) / R(y[i]-y[20] , y[j]-y[22] , exp1)
		 + (-GMura * y[j]) / R(y[i]-y[24] , y[j]-y[26] , exp1)
		 + (-GMplu * y[j]) / R(y[i]-y[28] , y[j]-y[30] , exp1)
		 + (-GMnep * y[j]) / R(y[i]-y[32] , y[j]-y[34] , exp1);

///Saturn---------------------------------------------------------------------------------------
i = 20;
j = 22;

//X-components
dydx[20] = y[21];
dydx[21] = (-GMsun * y[i]) / R(y[i]       , y[j]       , exp1)
		 + (-GMmer * y[i]) / R(y[i]-y[0]  , y[j]-y[2]  , exp1)
		 + (-GMven * y[i]) / R(y[i]-y[4]  , y[j]-y[6]  , exp1)
		 + (-GMear * y[i]) / R(y[i]-y[8]  , y[j]-y[10] , exp1)
		 + (-GMmar * y[i]) / R(y[i]-y[12] , y[j]-y[14] , exp1)
		 + (-GMjup * y[i]) / R(y[i]-y[16] , y[j]-y[18] , exp1)
		 //saturn removed
		 + (-GMura * y[i]) / R(y[i]-y[24] , y[j]-y[26] , exp1)
		 + (-GMplu * y[i]) / R(y[i]-y[28] , y[j]-y[30] , exp1)
		 + (-GMnep * y[i]) / R(y[i]-y[32] , y[j]-y[34] , exp1);

//Y-components
dydx[22] = y[23];
dydx[23] = (-GMsun * y[j]) / R(y[i]       , y[j]       , exp1)
		 + (-GMmer * y[j]) / R(y[i]-y[0]  , y[j]-y[2]  , exp1)
		 + (-GMven * y[j]) / R(y[i]-y[4]  , y[j]-y[6]  , exp1)
		 + (-GMear * y[j]) / R(y[i]-y[8]  , y[j]-y[10] , exp1)
		 + (-GMmar * y[j]) / R(y[i]-y[12] , y[j]-y[14] , exp1)
		 + (-GMjup * y[j]) / R(y[i]-y[16] , y[j]-y[18] , exp1)
		//saturn removed
		 + (-GMura * y[j]) / R(y[i]-y[24] , y[j]-y[26] , exp1)
		 + (-GMplu * y[j]) / R(y[i]-y[28] , y[j]-y[30] , exp1)
		 + (-GMnep * y[j]) / R(y[i]-y[32] , y[j]-y[34] , exp1);

///Uranus---------------------------------------------------------------------------------------
i = 24;
j = 26;

//X-components
dydx[24] = y[25];
dydx[25] = (-GMsun * y[i]) / R(y[i]       , y[j]       , exp1)
		 + (-GMmer * y[i]) / R(y[i]-y[0]  , y[j]-y[2]  , exp1)
		 + (-GMven * y[i]) / R(y[i]-y[4]  , y[j]-y[6]  , exp1)
		 + (-GMear * y[i]) / R(y[i]-y[8]  , y[j]-y[10] , exp1)
		 + (-GMmar * y[i]) / R(y[i]-y[12] , y[j]-y[14] , exp1)
		 + (-GMjup * y[i]) / R(y[i]-y[16] , y[j]-y[18] , exp1)
		 + (-GMsat * y[i]) / R(y[i]-y[20] , y[j]-y[22] , exp1)
		 //uranus removed
		 + (-GMplu * y[i]) / R(y[i]-y[28] , y[j]-y[30] , exp1)
		 + (-GMnep * y[i]) / R(y[i]-y[32] , y[j]-y[34] , exp1);

//Y-components
dydx[26] = y[27];
dydx[27] = (-GMsun * y[j]) / R(y[i]       , y[j]       , exp1)
		 + (-GMmer * y[j]) / R(y[i]-y[0]  , y[j]-y[2]  , exp1)
		 + (-GMven * y[j]) / R(y[i]-y[4]  , y[j]-y[6]  , exp1)
		 + (-GMear * y[j]) / R(y[i]-y[8]  , y[j]-y[10] , exp1)
		 + (-GMmar * y[j]) / R(y[i]-y[12] , y[j]-y[14] , exp1)
		 + (-GMjup * y[j]) / R(y[i]-y[16] , y[j]-y[18] , exp1)
		 + (-GMsat * y[j]) / R(y[i]-y[20] , y[j]-y[22] , exp1)
		 //uranus removed
		 + (-GMplu * y[j]) / R(y[i]-y[28] , y[j]-y[30] , exp1)
		 + (-GMnep * y[j]) / R(y[i]-y[32] , y[j]-y[34] , exp1);

///Pluto---------------------------------------------------------------------------------------
i = 28;
j = 30;

//X-components
dydx[28] = y[29];
dydx[29] = (-GMsun * y[i]) / R(y[i]       , y[j]       , exp1)
		 + (-GMmer * y[i]) / R(y[i]-y[0]  , y[j]-y[2]  , exp1)
		 + (-GMven * y[i]) / R(y[i]-y[4]  , y[j]-y[6]  , exp1)
		 + (-GMear * y[i]) / R(y[i]-y[8]  , y[j]-y[10] , exp1)
		 + (-GMmar * y[i]) / R(y[i]-y[12] , y[j]-y[14] , exp1)
		 + (-GMjup * y[i]) / R(y[i]-y[16] , y[j]-y[18] , exp1)
		 + (-GMsat * y[i]) / R(y[i]-y[20] , y[j]-y[22] , exp1)
		 + (-GMura * y[i]) / R(y[i]-y[24] , y[j]-y[26] , exp1)
		 //pluto removed
		 + (-GMnep * y[i]) / R(y[i]-y[32] , y[j]-y[34] , exp1);

//Y-components
dydx[30] = y[31];
dydx[31] = (-GMsun * y[j]) / R(y[i]       , y[j]       , exp1)
		 + (-GMmer * y[j]) / R(y[i]-y[0]  , y[j]-y[2]  , exp1)
		 + (-GMven * y[j]) / R(y[i]-y[4]  , y[j]-y[6]  , exp1)
		 + (-GMear * y[j]) / R(y[i]-y[8]  , y[j]-y[10] , exp1)
		 + (-GMmar * y[j]) / R(y[i]-y[12] , y[j]-y[14] , exp1)
		 + (-GMjup * y[j]) / R(y[i]-y[16] , y[j]-y[18] , exp1)
		 + (-GMsat * y[j]) / R(y[i]-y[20] , y[j]-y[22] , exp1)
		 + (-GMura * y[j]) / R(y[i]-y[24] , y[j]-y[26] , exp1)
		 //pluto removed
		 + (-GMnep * y[j]) / R(y[i]-y[32] , y[j]-y[34] , exp1);

///Neptune---------------------------------------------------------------------------------------
i = 32;
j = 34;

//X-components
dydx[32] = y[33];
dydx[33] = (-GMsun * y[i]) / R(y[i]       , y[j]       , exp1)
		 + (-GMmer * y[i]) / R(y[i]-y[0]  , y[j]-y[2]  , exp1)
		 + (-GMven * y[i]) / R(y[i]-y[4]  , y[j]-y[6]  , exp1)
		 + (-GMear * y[i]) / R(y[i]-y[8]  , y[j]-y[10] , exp1)
		 + (-GMmar * y[i]) / R(y[i]-y[12] , y[j]-y[14] , exp1)
		 + (-GMjup * y[i]) / R(y[i]-y[16] , y[j]-y[18] , exp1)
		 + (-GMsat * y[i]) / R(y[i]-y[20] , y[j]-y[22] , exp1)
		 + (-GMura * y[i]) / R(y[i]-y[24] , y[j]-y[26] , exp1)
		 + (-GMplu * y[i]) / R(y[i]-y[28] , y[j]-y[30] , exp1);
		//neptune removed

//Y-components
dydx[34] = y[35];
dydx[35] = (-GMsun * y[j]) / R(y[i]       , y[j]       , exp1)
		 + (-GMmer * y[j]) / R(y[i]-y[0]  , y[j]-y[2]  , exp1)
		 + (-GMven * y[j]) / R(y[i]-y[4]  , y[j]-y[6]  , exp1)
		 + (-GMear * y[j]) / R(y[i]-y[8]  , y[j]-y[10] , exp1)
		 + (-GMmar * y[j]) / R(y[i]-y[12] , y[j]-y[14] , exp1)
		 + (-GMjup * y[j]) / R(y[i]-y[16] , y[j]-y[18] , exp1)
		 + (-GMsat * y[j]) / R(y[i]-y[20] , y[j]-y[22] , exp1)
		 + (-GMura * y[j]) / R(y[i]-y[24] , y[j]-y[26] , exp1)
		 + (-GMplu * y[j]) / R(y[i]-y[28] , y[j]-y[30] , exp1);
		//neptune removed

}

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

/*

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

*/

					y = yout;

					mercDerivs(x,y,dydx);
			}

			//closes file output stream
			ofPositionMercury.close();
			ofPerihelion.close();

		cout << "\n" << "CoMpLeTe CoMpLeTe CoMpLeTe CoMpLeTe CoMpLeTe CoMpLeTe" << "\n";

	return 0;
}




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
	double kmax=365*12;  //max iterations (units: days)
	double h=1;       //time step size (units: days)

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
	ofAllData.open("all-data.csv");
	ofPositionEarth.open("x-y-positionEarth.csv");
	ofPositionJupiter.open("x-y-positionJupiter.csv");


	for(int k=0; k < kmax; k++)
	{
			x=xmin+k*h;

			rk4(y, dydx,  x, h, yout, derivs);

			//display output
			cout << "k = " << k
				 << "    Xe = " << yout[0] << "    X'e= " << yout[1]
				 << "    Ye = " << yout[2] << "    Y'e= " << yout[3] << endl

				 << "    Xj = " << yout[4] << "    X'j= " << yout[5]
				 << "    Yj = " << yout[6] << "    Y'j= " << yout[7] << endl << endl;

			//file output streams
			ofAllData << yout[0] << "," << yout[1] << "," << yout[2] << "," << yout[3] << yout[4] << "," << yout[5] << "," << yout[6] << ","  << yout[7] << endl;
			ofPositionEarth << yout[0] << "," << yout[2] << endl;
			ofPositionJupiter << yout[4] << "," << yout[6] << endl;

			y = yout;

			derivs(x,y,dydx);
	}

	//closes file output streams
	ofAllData.close();
	ofPositionEarth.close();
	ofPositionJupiter.close();

return 0;

}




int main()
{



	mercury();
	return 0;
}
























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


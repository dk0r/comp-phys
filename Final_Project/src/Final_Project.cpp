//============================================================================
// Name        : Final_Project.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>
#include <cmath>
#include <nr3.h>
#include <rk4.h>

using namespace std;

	//Pi
	double Pi = atan(1)*4;

	double g = 9.803; //gravitational acceleration near earth's surface

	//Mass1
	double M1 = 1;   //mass of Mass-1
	double R1 = 1; //pendulum length

	//Mass2
	double M2 = 1;	 //mass of Mass-2
	double R2 = 1; //pendulum length

	//Global Stepsize Counter (dt)
	double globStep = 0;



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
  //Mass-1
	//Angular Velocity
	dydx[0] = y[1];
	//Angular Acceleration
	dydx[1] = (  -g*(2*M1+M2)*sin(y[0])-M2*g*sin(y[0]-2*y[2])-2*sin(y[0]-y[2])*M2*(y[3]*y[3]*R2+y[1]*y[1]*R1*cos(y[0]-y[2]))  )
			/ (  R1*(2*M1+M2-M2*cos(2*y[0]-2*y[2]))  );

  //Mass-2
	//Angular Velocity
	dydx[2] = y[3];
	//Angular Acceleration
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
		int kmax=2000;  //max iterations (units: days)
		double h=0.02;       //time step size (units: days)

		///Initial Conditions:

	  //Mass-1
		//Angular-Position
		y[0] = (4*atan(1))+0.011;
		//Angular-Velocity
		y[1] = 0;

	  //Mass-2
		//Angular-Position
		y[2] = (4*atan(1))+0.011;
		//Angular-Velocity
		y[3] = .5;


		pendDerivs(xmin, y, dydx);

		//file output streams
		ofstream ofPositionM1, ofPositionM2, ofMidi;
		ofPositionM1.open("M1_position.csv");
		ofPositionM2.open("M2_position.csv");
		ofMidi.open("midi.csv");









		for(int k=0; k < kmax; k++)
		{

			//accuracy for determining when axes cross
			double xAcc = 0.009;
			double yAcc = 0.009;

			//track time-stepping for midi note play time

			//S!@#TG@$HG@$GH@$@ 0.002 SeConDs ?!?!?!?!?
			globStep += h; //S!@#TG@$HG@$GH@$@ 0.002 SeConDs ?!?!?!?!?//S!@#TG@$HG@$GH@$@ 0.002 SeConDs ?!?!?!?!?
			//S!@#TG@$HG@$GH@$@ 0.002 SeConDs ?!?!?!?!?

			//?????
			x=xmin+k*h;

			//performs runge-kutta4
			rk4(y, dydx,  x, h, yout, pendDerivs);

			/*
			//display output
			cout << "k = " << k
				 << "    Xe = " << yout[0] << "    X'e= " << yout[1]
				 << "    Ye = " << yout[2] << "    Y'e= " << yout[3] << endl

				 << "    Xj = " << yout[4] << "    X'j= " << yout[5]
				 << "    Yj = " << yout[6] << "    Y'j= " << yout[7] << endl << endl;
			 */


			//Stores x and y positions of Mass-1 & Mass-2
			double x1 =  R1*sin(y[0]);
			double y1 = -R1*cos(y[0]);
			double x2 = ( R1*sin(y[0]) + R2*sin(y[2]) );
			double y2 = ( -R1*cos(y[0]) - R2*cos(y[2]) );



			//midi Counter for midi[midiCount][] matrix;
			int midiCount = 0;


			//Midi matrix
			double midi[kmax][6];

			//Initializes Midi matrix
			for(int i=0; i<kmax; i++)
			{
				for(int j=0; j<6; j++)
				{
					midi[i][j] = 0;
				}
			}





			//Number of Octaves (limits range of scale)
			int octaves = 1;

			//Total number of notes
			int notes = octaves*12 + 1;


			//interval counter for knowing which interval to relate to the Notes[i] matrix
			int intCounter = 0;

			//Notes matrix filled w/ note#'s
			int Notes[notes];

			//Defines notes in matrix
			for(int i=0; i<notes; i++)
			{
				Notes[i] = 60 - i; //Middle C @ top to something Low
			}




			//holds previous iterations value of hint (the 1/2 interval size)
			double prev = 0;





			//Length of side of pendulum container
			double side = 2*(R1+R2);
			double halfSide = R1+R2;


			//Note interval size
			double interval = side/notes;
			//double halfInt = (R1+R2)/notes;




			//Determines whether X or Y axes crossed and assigns note# depending on
			//interval where crossing occurred.

			for(double i=halfSide; i >= -1*halfSide; i = i - interval )
			{

				//X-Axis Crossings (Both masses have ~~ same Y positions)
				if(fabs(y1-y2) <= xAcc)
				{

					cout << "X-Axis Crossing where y1 ~~ y2 @ " << y1 << "," << y2 << "\n";

					if(y1 <= i  &&  y1 >= i-interval)
					{
						//Defines midi note found in interval in above if() line.

						midi[midiCount][0] = 1; //track #
						midi[midiCount][1] = 1; //channel #
						midi[midiCount][2] = Notes[intCounter]; //midi note #
					  //midi[midiCount][3] volume velocity specified in matlab
						midi[midiCount][4] = globStep; //time which note will be played
						midi[midiCount][5] = 0.5; //duration of each note
						midiCount++;
					}

				}




				//Y-Axis Crossings (Both masses have ~~ same X positions)
				if(fabs(x1-x2) <= yAcc)
				{
					cout << "Y-Axis Crossing where x1 ~~ x2 @ " << x1 << "," << x2 << "\n";

					if(x1 < i  &&  x1 >= i-interval)
					{
						//Defines midi note found in interval in above if() line.

						midi[midiCount][0] = 1; //track #
						midi[midiCount][1] = 1; //channel #
						midi[midiCount][2] = Notes[intCounter]; //midi note #
					  //midi[midiCount][3] volume velocity specified in matlab
						midi[midiCount][4] = globStep; //time which note will be played
						midi[midiCount][5] = 0.5; //duration of each note
						midiCount++;
					}


				}

				//MoVe ThIs prev=i; To BoTtOm Of ThIs LoOOoOooOOOOOooOOOOOoOOOp
				prev = i;
				intCounter++;

			}


			//file output streams
			ofPositionM1 << R1*sin(y[0]) << "," << -R1*cos(y[0]) << endl;
			ofPositionM2 << R1*sin(y[0]) + R2*sin(y[2]) << "," << (-R1*cos(y[0]) - R2*cos(y[2])) << endl;
			for(int i=0; i<kmax; i++)
			{
				if(midi[i][0] != 1)
				{
					ofMidi << midi[i][0] << "," << midi[i][1] << "," << midi[i][2] << "," << midi[i][3] << "," << midi[i][4] << "," << midi[i][5] << endl;
				}
			}


			y = yout;

			pendDerivs(x,y,dydx);
		}

		//closes file output streams
		ofPositionM1.close();
		ofPositionM2.close();
		ofMidi.close();
}

int main()
{

	doublePend();
	return 0;

}

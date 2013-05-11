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


	double h=0.01;       //time step size      0.01
	int kmax=3000;  //max iterations        1000      100 in Mathematica


	//accuracy for determining when axes cross
	double crossAcc = 0.05;
	double repeatAcc = 0.05;
	double timeAcc = 0.2;

	//
	double prevX = 0;
	double prevY = 0;
	double prevTime = 0;

	//Number of Octaves (limits range of scale)
	int octaves = 1;
	//Total number of notes
	int notes = octaves*12 + 1;

	//Length of side of pendulum container
	double radius = R1+R2;

	//Note interval size
	double interval = 2*radius/notes;



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




double xNotes(double pos)
{
	int note = 91;
	double y = pos;

	cout << "xNotes @  y ~ " << pos << endl;

	if(y <= radius && y >= radius - interval*0 - interval*1)
		note = 73;

	if(y < radius - interval*1 && y >= radius - interval*2)
		note = 45;

	if(y < radius - interval*2 && y >= radius - interval*3)
		note = 44;

	if(y < radius - interval*3 && y >= radius - interval*4)
		note = 37;

	if(y < radius - interval*4 && y >= radius - interval*5)
		note = 64; //37

	if(y < radius - interval*5 && y >= radius - interval*6)
		note = 63;




	if(y < radius - interval*6 && y >= radius - interval*7)
		note = 61;




	if(y < radius - interval*7 && y >= radius - interval*8)
		note = 57;

	if(y < radius - interval*8 && y >= radius - interval*9)
		note = 56;

	if(y < radius - interval*9 && y >= radius - interval*10)
		note = 54;

	if(y < radius - interval*10 && y >= radius - interval*11)
		note = 52;

	if(y < radius - interval*11 && y >= radius - interval*12)
		note = 51;

	if(y < radius - interval*12 && y >= radius - interval*13)
		note = 68;

		return note;
}

double yNotes(double pos)
{
	double x = pos;

	int note = 91;


	if(x >= 0)
	{


		if(x >= 5.5*interval && x < radius )
		note = 69;


		if(x >= 4.5*interval && x < (radius - interval) )
		note = 69;


		if(x >= 3.5*interval && x < (radius - 2*interval) )
		note = 68;


		if(x >= 2.5*interval && x < (radius - 3*interval) )
		note = 37;


		if(x >= 1.5*interval && x < (radius - 4*interval) )
		note = 64;

		if(x >= 0.5*interval && x < (radius - 5*interval) )
		note = 63;


		if(x >= 0 && x < (radius - 6*interval) )
		note = 61;

		return note;
	}

	if(x < 0)
	{
		x = fabs(x);


		if(x >= 0 && x < (radius - 6*interval) )
		note = 61;

		if(x >= 0.5*interval && x < (radius - 5*interval) )
		note = 57;

		if(x >= 1.5*interval && x < (radius - 4*interval) )
		note = 56;

		if(x >= 2.5*interval && x < (radius - 3*interval) )
		note = 54;

		if(x >= 3.5*interval && x < (radius - 2*interval) )
		note = 52;

		if(x >= 4.5*interval && x < (radius - interval) )
		note = 51;

		if(x >= 5.5*interval && x < radius )
		note = 49;

		return note;
	}

}



void doublePend()
{
		VecDoub y(4);
		VecDoub dydx(4); //vector of positions & velocities for earth and
		VecDoub yout(4);

		double x;
		double xmin = 0;      //minimum starting position

	///Initial Conditions:

	  //Mass-1
		//Angular-Position
		y[0] = (3)*atan(1);//(4*atan(1))+0.011;
		//Angular-Velocity
		y[1] = 0;

	  //Mass-2
		//Angular-Position
		y[2] = (3)*atan(1);//(4*atan(1))+0.011;
		//Angular-Velocity
		y[3] = 0;


		pendDerivs(xmin, y, dydx);

		//file output streams
		ofstream ofPositionM1, ofPositionM2, ofMidi;
		ofPositionM1.open("M1_position.csv");
		ofPositionM2.open("M2_position.csv");
		ofMidi.open("midi.csv");

		//midi Counter for midi[midiCount][] matrix;
		int midiCount = 0;



		for(int k=0; k < kmax; k++)
		{

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

			//Midi matrix
			double midi[kmax][6];


			//Determines whether X or Y axes crossed and assigns note# depending on
			//interval where crossing occurred.

			//tracks time-stepping for midi note play time
			globStep += h;


		//X-Axis Crossings (Both masses have ~~ same Y positions)
			if( (fabs(y1-y2) <= crossAcc) && (fabs(globStep-prevTime) > timeAcc) )
				{
					prevTime = globStep;

					//	cout << "X-Cross @" << globStep << "steps. y1~y2 @ " << y1 << "," << y2 << "\n";
					//	cout <<"Entering X-Axis Midi.  midiCount = " << midiCount << " globStep = " << globStep
					//		 << "  midi-note: " << xNotes[y1] << endl << endl;
						//Defines midi note found in interval in above if() line.

						midi[midiCount][0] = 1; //track #
						midi[midiCount][1] = 1; //channel #
						midi[midiCount][2] = xNotes(x2); //midi note #
					    midi[midiCount][3] = 0; //volume velocity: to be specified in matlab
						midi[midiCount][4] = globStep; //time which note will be played
						midi[midiCount][5] = 0.5; //duration of each note

						ofMidi << midi[midiCount][0] << "," << midi[midiCount][1] << "," << midi[midiCount][2] << "," << midi[midiCount][3] << "," << midi[midiCount][4] << "," << midi[midiCount][5] << endl;

						midiCount++;

				}


			//Y-Axis Crossings (Both masses have ~~ same X positions)
				if( (fabs(x1-x2) <= crossAcc) && (fabs(globStep-prevTime) > timeAcc) )
				{
					prevTime = globStep;

						cout << "Y-Cross @" << globStep << "steps. x1~x2 @ " << x1 << "," << x2 << "\n";
						cout <<"Entering Y-Axis Midi.  midiCount = " << midiCount << " globStep = " << globStep
								<< "  midi-note: " << yNotes(y1) << endl << endl;
						//Defines midi note found in interval in above if() line.

						midi[midiCount][0] = 1; //track #
						midi[midiCount][1] = 1; //channel #
						midi[midiCount][2] = yNotes(y2); //midi note #
					  //midi[midiCount][3] volume velocity specified in matlab
						midi[midiCount][4] = globStep; //time which note will be played
						midi[midiCount][5] = 0.5; //duration of each note

						ofMidi << midi[midiCount][0] << "," << midi[midiCount][1] << "," << midi[midiCount][2] << "," << midi[midiCount][3] << "," << midi[midiCount][4] << "," << midi[midiCount][5] <<  endl;

						midiCount++;

				}


			//file output streams
			ofPositionM1 << R1*sin(y[0]) << "," << -R1*cos(y[0]) << endl;
			ofPositionM2 << R1*sin(y[0]) + R2*sin(y[2]) << "," << (-R1*cos(y[0]) - R2*cos(y[2])) << endl;


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

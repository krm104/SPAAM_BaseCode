/******************************************************************************
This program provides a very simple example of how the Singular Value
Decomposition algorithm can be used to solve for a 3x4 Projection matrix
using the SPAAM algorithm.

The 2 functions 'buildVt' and 'buildVf' are simply used to fill the list of
Correspondence_Pair objects. In your own implementation, these functions would 
most likely not be used, but you would instead create your own Correpsondence_Pair
objects based on the 2D location of screen pixels and the 3D location of a world point.

The main function simply calles the procedures to build the Correspondence_Pair
objects, load them into a full list (vector) and then perform the SVD operation
on the list.

The final 3x4 projection matrix is then printed to the console window.
******************************************************************************/

#include <vector>
#include <iostream>

#include "SPAAM_SVD.h"

using namespace std;

/************************************
Simple function to build a series of
2D point locations (dummy screen locations)
************************************/
void buildVt(vector<ublas::c_vector<double, 2>> &vt)
{
	ublas::c_vector<double, 2> v2;
	v2(0) = 192;
	v2(1) = 96;
	vt.push_back(v2);
	v2(0) = 576;
	v2(1) = 96;
	vt.push_back(v2);
	v2(0) = 960;
	v2(1) = 96;
	vt.push_back(v2);
	v2(0) = 1344;
	v2(1) = 96;
	vt.push_back(v2);
	v2(0) = 1728;
	v2(1) = 96;
	vt.push_back(v2);
	v2(0) = 192;
	v2(1) = 312;
	vt.push_back(v2);
	v2(0) = 576;
	v2(1) = 312;
	vt.push_back(v2);
	v2(0) = 960;
	v2(1) = 312;
	vt.push_back(v2);
	v2(0) = 1344;
	v2(1) = 312;
	vt.push_back(v2);
	v2(0) = 1728;
	v2(1) = 312;
	vt.push_back(v2);
	v2(0) = 192;
	v2(1) = 528;
	vt.push_back(v2);
	v2(0) = 576;
	v2(1) = 528;
	vt.push_back(v2);
	v2(0) = 960;
	v2(1) = 528;
	vt.push_back(v2);
	v2(0) = 1344;
	v2(1) = 528;
	vt.push_back(v2);
	v2(0) = 1728;
	v2(1) = 528;
	vt.push_back(v2);
	v2(0) = 192;
	v2(1) = 744;
	vt.push_back(v2);
	v2(0) = 576;
	v2(1) = 744;
	vt.push_back(v2);
	v2(0) = 960;
	v2(1) = 744;
	vt.push_back(v2);
	v2(0) = 1344;
	v2(1) = 744;
	vt.push_back(v2);
	v2(0) = 1728;
	v2(1) = 744;
	vt.push_back(v2);
	v2(0) = 192;
	v2(1) = 960;
	vt.push_back(v2);
	v2(0) = 576;
	v2(1) = 960;
	vt.push_back(v2);
	v2(0) = 960;
	v2(1) = 960;
	vt.push_back(v2);
	v2(0) = 1344;
	v2(1) = 960;
	vt.push_back(v2);
	v2(0) = 1728;
	v2(1) = 960;
	vt.push_back(v2);
	v2(0) = 192;
	v2(1) = 96;
	vt.push_back(v2);
	v2(0) = 576;
	v2(1) = 96;
	vt.push_back(v2);
	v2(0) = 960;
	v2(1) = 96;
	vt.push_back(v2);
	v2(0) = 1344;
	v2(1) = 96;
	vt.push_back(v2);
	v2(0) = 1728;
	v2(1) = 96;
	vt.push_back(v2);

	//cout << vt.size() << endl;
}

/************************************
Simple function to build a series of
3D point locations (dummy world locations)
************************************/
void buildVf(vector<ublas::c_vector<double, 3>> &vf)
{
	ublas::c_vector<double, 3> v1;
	v1(0) = -0.022621147155761717;
	v1(1) = 0.38774725341796873;
	v1(2) = 0.1269190673828125;
	vf.push_back(v1);
	v1(0) = -0.050548507690429685;
	v1(1) = 0.48674676513671877;
	v1(2) = 0.15059318542480468;
	vf.push_back(v1);
	v1(0) = -0.071583969116210941;
	v1(1) = 0.35663726806640628;
	v1(2) = 0.11513092041015625;
	vf.push_back(v1);
	v1(0) = -0.11799789428710937;
	v1(1) = 0.47049493408203125;
	v1(2) = 0.14779368591308595;
	vf.push_back(v1);
	v1(0) = -0.12474505615234376;
	v1(1) = 0.35620223999023437;
	v1(2) = 0.11392384338378907;
	vf.push_back(v1);
	v1(0) = -0.019622428894042969;
	v1(1) = 0.46816082763671873;
	v1(2) = 0.12061217498779297;
	vf.push_back(v1);
	v1(0) = -0.048498607635498049;
	v1(1) = 0.32230554199218753;
	v1(2) = 0.095608802795410161;
	vf.push_back(v1);
	v1(0) = -0.079887237548828122;
	v1(1) = 0.45576977539062502;
	v1(2) = 0.13440872192382813;
	vf.push_back(v1);
	v1(0) = -0.090922668457031244;
	v1(1) = 0.3146878967285156;
	v1(2) = 0.092620956420898432;
	vf.push_back(v1);
	v1(0) = -0.14409155273437499;
	v1(1) = 0.466261474609375;
	v1(2) = 0.12376845550537109;
	vf.push_back(v1);
	v1(0) = -0.017212881088256835;
	v1(1) = 0.30422473144531248;
	v1(2) = 0.082259178161621094;
	vf.push_back(v1);
	v1(0) = -0.045188476562500003;
	v1(1) = 0.39570046997070313;
	v1(2) = 0.10126533508300781;
	vf.push_back(v1);
	v1(0) = -0.070368377685546879;
	v1(1) = 0.33771701049804687;
	v1(2) = 0.086266128540039066;
	vf.push_back(v1);
	v1(0) = -0.10274946594238281;
	v1(1) = 0.39483236694335938;
	v1(2) = 0.098369087219238288;
	vf.push_back(v1);
	v1(0) = -0.11732191467285157;
	v1(1) = 0.34514050292968751;
	v1(2) = 0.08450955200195312;
	vf.push_back(v1);
	v1(0) = -0.020879354476928711;
	v1(1) = 0.43182464599609377;
	v1(2) = 0.088839363098144525;
	vf.push_back(v1);
	v1(0) = -0.046161254882812502;
	v1(1) = 0.3972061767578125;
	v1(2) = 0.071852272033691406;
	vf.push_back(v1);
	v1(0) = -0.067129570007324224;
	v1(1) = 0.34340292358398439;
	v1(2) = 0.075639427185058589;
	vf.push_back(v1);
	v1(0) = -0.10812483978271484;
	v1(1) = 0.44349353027343752;
	v1(2) = 0.083343391418457036;
	vf.push_back(v1);
	v1(0) = -0.10961792755126953;
	v1(1) = 0.31100454711914061;
	v1(2) = 0.061563499450683595;
	vf.push_back(v1);
	v1(0) = -0.021026494979858399;
	v1(1) = 0.45146182250976563;
	v1(2) = 0.076992515563964845;
	vf.push_back(v1);
	v1(0) = -0.042361728668212889;
	v1(1) = 0.31638238525390627;
	v1(2) = 0.054041759490966797;
	vf.push_back(v1);
	v1(0) = -0.072913192749023442;
	v1(1) = 0.4203233947753906;
	v1(2) = 0.058085742950439452;
	vf.push_back(v1);
	v1(0) = -0.091090324401855466;
	v1(1) = 0.31261776733398439;
	v1(2) = 0.04887359619140625;
	vf.push_back(v1);
	v1(0) = -0.13445088195800781;
	v1(1) = 0.42511129760742189;
	v1(2) = 0.060736068725585936;
	vf.push_back(v1);
	v1(0) = -0.016986076354980468;
	v1(1) = 0.27146170043945311;
	v1(2) = 0.10301142883300782;
	vf.push_back(v1);
	v1(0) = -0.038408023834228519;
	v1(1) = 0.35083627319335936;
	v1(2) = 0.11071772003173828;
	vf.push_back(v1);
	v1(0) = -0.061872913360595701;
	v1(1) = 0.28067010498046874;
	v1(2) = 0.091477111816406248;
	vf.push_back(v1);
	v1(0) = -0.10767545318603515;
	v1(1) = 0.4559621887207031;
	v1(2) = 0.14422778320312499;
	vf.push_back(v1);
	v1(0) = -0.11666415405273438;
	v1(1) = 0.34573306274414062;
	v1(2) = 0.11103826141357422;
	vf.push_back(v1);

	//cout << vf.size() << endl;
	  
}

/****************************************************************
Main Function and starting point of the program

Calls the functions to create the Correspondence_Pair objects
Iterates over the objects and loads them into an instance of 
a SPAAM_SVD object.
Calls the SVD operation on the SPAAM_SVD object and prints the
resulting 3x4 matric to the console window.
****************************************************************/
int main(int argc, char** argv)
{
	//create vectors to store Correspondence_Pairs//
	vector<ublas::c_vector<double, 3>> vf;
	vector<ublas::c_vector<double, 2>> vt;

	//Build sets of 2D and 3D points//
	buildVf(vf);
	buildVt(vt);
	
	//Instanciate SPAAM_SVD object//
	SPAAM_SVD<double, double> nsvd;

	//Put 2D and 3D points into the SPAAM_SVD object//
	for ( unsigned int i = 0; i < vf.size( ); i++ )
	{
		Correspondence_Pair<double, double> ncp(vf[i](0), vf[i](1), vf[i](2), vt[i](0), vt[i](1));
		nsvd.corr_points.push_back(ncp);
	}

	//Call the SVD operation and print the result//
	cout << nsvd.projectionDLTImpl( ) << endl;

	return 0;
}

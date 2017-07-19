#include <math.h>
#include <cstring>
#include <iostream>
#include "Utils.hpp"

class Derivatives{


    public:

	int Nx, Ny;
	double dx, dy;
	
	double alpha;
	double beta;
	double a, b, c;

	double *diagx, *offx;
	double *diagy, *offy;

    Derivatives(){
	Nx = 0; Ny = 0;
	dx = 0; dy = 0;

	alpha = 1.0/3.0;
	beta  = 0.0;
	a = 14.0/9.0;
	b = 1.0/9.0;
	c = 0.0;

	diagx = NULL; offx = NULL;
	diagy = NULL; offy = NULL;
    }

    Derivatives(int NX, int NY){
	Nx = NX; Ny = NY;

	alpha = 1.0/3.0;
	beta  = 0.0;
	a = 14.0/9.0;
	b = 1.0/9.0;
	c = 0.0;

	diagx = new double[NX]; 
	offx  = new double[NX];
	diagy = new double[NY]; 
	offy  = new double[NY];
    }

    void multRHSDeriv(double dh, double *phi, int N, double *RHSvec);
    void CompactDYPeriodic(double *phi, double *dphidy);
    void CompactDXPeriodic(double *phi, double *dphidx);

};

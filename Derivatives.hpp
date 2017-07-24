#include <math.h>
#include <cstring>
#include <iostream>
#include "Utils.hpp"

class Derivatives{


    public:

	int Nx, Ny;
	double dx, dy;

        //Interior Coefficients	
	double alpha;
	double beta;
	double a, b, c;

	//T6 Dirichlet Coefficients w/ T4 at edge
	double alpha1, alpha21, alpha22;
	double a1, b1, c1, d1, e1, f1;
	double a2, b2, c2, d2, e2; 

	double *diagx, *offx1, *offx2;
	double *diagy, *offy1, *offy2;

    Derivatives(){
	Nx = 0; Ny = 0;
	dx = 0; dy = 0;

	alpha = 1.0/3.0;
	beta  = 0.0;
	a = 14.0/9.0;
	b = 1.0/9.0;
	c = 0.0;

	alpha1 = 5.0;
	a1 = -197.0/60.0;
	b1 =   -5.0/12.0;
	c1 =    5.0;
	d1 =   -5.0/3.0;
	e1 =    5.0/12.0;
	f1 =   -1.0/20.0;

	alpha21 = 1.0/8.0;
	alpha22 = 3.0/4.0;
	a2 = -43.0/96.0;
	b2 = -5.0/6.0;
	c2 =  9.0/8.0;
	d2 =  1.0/6.0;
	e2  = -1.0/96.0;

	diagx = NULL; offx1 = NULL; offx2 = NULL;
	diagy = NULL; offy1 = NULL; offy2 = NULL;
    }

    Derivatives(int NX, int NY, double DX, double DY){
	Nx = NX; Ny = NY;
	dx = DX; dy = DY;

	alpha = 1.0/3.0;
	beta  = 0.0;
	a = 14.0/9.0;
	b = 1.0/9.0;
	c = 0.0;

	alpha1 = 5.0;
	a1 = -197.0/60.0;
	b1 =   -5.0/12.0;
	c1 =    5.0;
	d1 =   -5.0/3.0;
	e1 =    5.0/12.0;
	f1 =   -1.0/20.0;

	alpha21 = 1.0/8.0;
	alpha22 = 3.0/4.0;
	a2 = -43.0/96.0;
	b2 = -5.0/6.0;
	c2 =  9.0/8.0;
	d2 =  1.0/6.0;
	e2  = -1.0/96.0;

	diagx = new double[NX]; 
	offx1 = new double[NX];
	offx2 = new double[NX];
	diagy = new double[NY]; 
	offy1 = new double[NY];
	offy2 = new double[NY];

	for(int ip = 0; ip < Nx; ip++){
	    diagx[ip] = 1.0;
	    offx1[ip]  = alpha;
	    offx2[ip]  = alpha;
	}

	for(int ip = 0; ip < Ny; ip++){
	    diagy[ip] = 1.0;
	    offy1[ip]  = alpha;
	    offy2[ip]  = alpha;
	}

    }


    void multRHSDerivPeriodic(double dh, double *phi, int N, double *RHSvec);
    void CompactDYPeriodic(double *phi, double *dphidy);
    void CompactDXPeriodic(double *phi, double *dphidx);

    void multRHSDerivDirichlet(double dh, double *phi, int N, double *RHSvec);
    void CompactDYDirichlet(double *phi, double *dphidy);
    void CompactDXDirichlet(double *phi, double *dphidx);
};

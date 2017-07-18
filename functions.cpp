#include <omp.h>
#include <math.h>
#include <cstring>
#include <iostream>
#include "functions.h"

using namespace std;

void solvemu(double *T, double T_ref, double mu_ref, int Nx, int Ny, double *mu){
    #pragma omp parallel for 
    for(int ip = 0; ip < Ny*Nx; ip++){
	mu[ip] = mu_ref*pow(T[ip]/T_ref, 0.76);
    }
}

void solveU(double *rho, double *rhoU, int Nx, int Ny, double *U){
    #pragma omp parallel for 
    for(int ip = 0; ip < Ny*Nx; ip++){
	U[ip] = rhoU[ip]/rho[ip];
    }
}

void solvep(double *rho, double *rhoE, double *U, double *V, double gamma, int Nx, int Ny, double *p){
    #pragma omp parallel for 
    for(int ip = 0; ip < Ny*Nx; ip++){
	p[ip] = (gamma-1)*(rhoE[ip] - 0.5 * rho[ip]*(U[ip]*U[ip] + V[ip]*V[ip]));
    }
}

void solveT(double *rho, double *p, double R_gas, int Nx, int Ny, double *T){
    #pragma omp parallel for 
    for(int ip = 0; ip < Ny*Nx; ip++){
	T[ip] = p[ip]/(rho[ip]*R_gas);
    }
}

void solveSOS(double *rho, double *p, double gamma, int Nx, int Ny, double *SOS){
    #pragma omp parallel for 
    for(int ip = 0; ip < Ny*Nx; ip++){
	SOS[ip] = sqrt(gamma*p[ip]/rho[ip]);
    }
}

void solveTri(double a[], double b[], double c[], double d[], double x[], double *work, int size)
{
	memcpy(work, b, size*sizeof(double));
	memcpy(x, d, size*sizeof(double));

	for(int ip = 1; ip < size; ip++){
	    double m = a[ip]/work[ip-1];
	    work[ip] = work[ip] - m*c[ip-1];
	    x[ip] = x[ip] - m*x[ip-1];
	}

	x[size-1] /= work[size-1];
	for(int ip = size-2; ip >= 0; ip--){
	    x[ip] = (x[ip] - c[ip]*x[ip+1])/work[ip];
	}


}

void cyclic(double *a, double *b, double *c, double alpha, double beta,
	double *r, int n, double *x)
{
	unsigned long i;
	double fact,gamma,*bb,*u,*z;

	if (n <= 2) cout << "n too small in cyclic" << endl;

	bb = new double[n];
	u  = new double[n];
	z  = new double[n];

	gamma = -b[0];
	for (i=0;i<n;i++) bb[i]=b[i];
	bb[0]=b[0]-gamma;
	bb[n-1]=b[n-1]-alpha*beta/gamma;

	double *work = new double[n];
  	solveTri(a,bb,c,r,x,work,n);

	for (i=0;i<n;i++) u[i]=0.0;
	u[0]=gamma;
	u[n-1]=alpha;

	solveTri(a,bb,c,u,z,work,n);

	fact=(x[0]+beta*x[n-1]/gamma)/(1.0+z[0]+beta*z[n-1]/gamma);

	for (i=0;i<n;i++) x[i] -= fact*z[i];

	delete[] z;
	delete[] u;
	delete[] bb;
	delete[] work;
}

void multRHSDeriv(double a, double b, double c, double dh, double *phi, int N, double *RHSvec){

    

    double c1 = -c/6.0;
    double c2 = -b/4.0;
    double c3 = -a/2.0;
    double c4 =  a/2.0;
    double c5 =  b/4.0;
    double c6 =  c/6.0;

    RHSvec[0] = c1*phi[N-3] + c2*phi[N-2] + c3*phi[N-1] + \
			c4*phi[1] + c5*phi[2] + c6*phi[3];
    RHSvec[1] = c1*phi[N-2] + c2*phi[N-1] + c3*phi[0] + \
			c4*phi[2] + c5*phi[3] + c6*phi[4];
    RHSvec[2] = c1*phi[N-1] + c2*phi[0] + c3*phi[1] + \
			c4*phi[3] + c5*phi[4] + c6*phi[5];

    for(int ip = 0; ip < N; ip++){
        RHSvec[ip] = c1*phi[ip-3] + c2*phi[ip-2] + c3*phi[ip-1] + \
			c4*phi[ip+1] + c5*phi[ip+2] + c6*phi[ip+3]; 

    } 

    RHSvec[N-3] = c1*phi[N-6] + c2*phi[N-5] + c3*phi[N-4] + \
			c4*phi[N-2] + c5*phi[N-1] + c6*phi[0]; 
    RHSvec[N-2] = c1*phi[N-5] + c2*phi[N-4] + c3*phi[N-3] + \
			c4*phi[N-1] + c5*phi[0] + c6*phi[1]; 
    RHSvec[N-1] = c1*phi[N-4] + c2*phi[N-3] + c3*phi[N-2] + \
			c4*phi[0] + c5*phi[1] + c6*phi[2]; 

    for(int ip = 0; ip < N; ip++){
        RHSvec[ip] /= dh;
    }

}

void CompactDYPeriodic(double *phi, double *diagy, double *offy, double alpha, double a, double b, double c, double dy, int Nx, int Ny, double *dphidy){


    double RHSvec[Ny];
    #pragma omp parallel for private(RHSvec)
    for(int ip = 0 ; ip < Nx; ip++){
	double *phiPointer = &phi[ip*Ny];
	multRHSDeriv(a,b,c,dy, phiPointer, Ny, RHSvec);
	cyclic(offy, diagy, offy, alpha, alpha, RHSvec, Ny, &dphidy[ip*Ny]);
    }

}

void transposeMatrix(double *in, int Nx, int Ny, double *out){

    #pragma omp parallel for collapse(2)  
    for(int i = 0; i < Ny; ++i)
        for(int j = 0; j < Nx; ++j)
            out[i*Nx + j] =  in[j*Ny + i];

}

void transposeMatrix_Fast1(const double *in, int n, int p, double *out, int block){
    #pragma omp parallel for 
    for (int i = 0; i < n; i += block) {
        for(int j = 0; j < n; ++j) {
            for(int b = 0; b < block && i + b < n; ++b) {
                out[j*n + i + b] = in[(i + b)*n + j];
            }
        }
    }
}

void transposeMatrix_Fast2(const double *in, int n, int p, double *out, int blocksize){
    int i, j, row, col;
    #pragma omp parallel for private(i, j, row, col) collapse(2) // schedule(static, 2)
    for ( i = 0; i < n; i += blocksize) {
        for ( j = 0; j < p; j += blocksize) {
            for (row = i; row < i + blocksize && row < n; row++) {
                for (col = j; col < j + blocksize && col < p; col++) {
                    out[row*p + col] = in[col*n + row];
                }
            }
        }
    }
}

void CompactDXPeriodic(double *phi, double *diagx, double *offx, double alpha, double a, double b, double c, double dx, int Nx, int Ny, double *dphidx){

    double *phiTrans = new double[Nx*Ny];
    double *dphidxTrans = new double[Nx*Ny];
    //transposeMatrix(phi, Nx, Ny, phiTrans);
    transposeMatrix_Fast2(phi, Nx, Ny, phiTrans, 60);

    double RHSvec[Nx];
    #pragma omp parallel for private(RHSvec)
    for(int ip = 0 ; ip < Ny; ip++){
	double *phiPointer = &phiTrans[ip*Nx];
	multRHSDeriv(a,b,c,dx, phiPointer, Nx, RHSvec);
	cyclic(offx, diagx, offx, alpha, alpha, RHSvec, Nx, &dphidxTrans[ip*Nx]);
    }
    delete[] phiTrans;

    //transposeMatrix(dphidxTrans, Ny, Nx, dphidx);
    transposeMatrix_Fast2(dphidxTrans, Ny, Nx, dphidx, 60);
    delete[] dphidxTrans;
}

void multRHSFilter(double a0, double a1, double a2, double a3, double a4, double a5, double *phi, int N, double *RHSvec){


    double c0 = a0;
    double c1 = a1/2.0;
    double c2 = a2/2.0;
    double c3 = a3/2.0;
    double c4 = a4/2.0;
    double c5 = a5/2.0;

    for(int ip = 0; ip < N; ip++){
        if(ip == 0){
            RHSvec[ip] = c5*phi[N-5]  + c4*phi[N-4]  + c3*phi[N-3]  +
			 c2*phi[N-2]  + c1*phi[N-1]  + c0*phi[ip]   +
			 c1*phi[ip+1] + c2*phi[ip+2] + c3*phi[ip+3] +
			 c4*phi[ip+4] + c5*phi[ip+5];
        }else if(ip == 1){
            RHSvec[ip] = c5*phi[N-4]  + c4*phi[N-3]  + c3*phi[N-2]  +
			 c2*phi[N-1]  + c1*phi[ip-1] + c0*phi[ip]   +
			 c1*phi[ip+1] + c2*phi[ip+2] + c3*phi[ip+3] +
			 c4*phi[ip+4] + c5*phi[ip+5];
        }else if(ip == 2){
            RHSvec[ip] = c5*phi[N-3]  + c4*phi[N-2]  + c3*phi[N-1]  +
			 c2*phi[ip-2] + c1*phi[ip-1] + c0*phi[ip]   +
			 c1*phi[ip+1] + c2*phi[ip+2] + c3*phi[ip+3] +
			 c4*phi[ip+4] + c5*phi[ip+5];
        }else if(ip == 3){
            RHSvec[ip] = c5*phi[N-2]  + c4*phi[N-1]  + c3*phi[ip-3] +
			 c2*phi[ip-2] + c1*phi[ip-1] + c0*phi[ip]   +
			 c1*phi[ip+1] + c2*phi[ip+2] + c3*phi[ip+3] +
			 c4*phi[ip+4] + c5*phi[ip+5];
        }else if(ip == 4){
            RHSvec[ip] = c5*phi[N-1]  + c4*phi[ip-4] + c3*phi[ip-3] +
			 c2*phi[ip-2] + c1*phi[ip-1] + c0*phi[ip]   +
			 c1*phi[ip+1] + c2*phi[ip+2] + c3*phi[ip+3] +
			 c4*phi[ip+4] + c5*phi[ip+5];
        }else if(ip == N-5){
            RHSvec[ip] = c5*phi[ip-5] + c4*phi[ip-4] + c3*phi[ip-3] +
			 c2*phi[ip-2] + c1*phi[ip-1] + c0*phi[ip]   +
			 c1*phi[ip+1] + c2*phi[ip+2] + c3*phi[ip+3] +
			 c4*phi[ip+4] + c5*phi[0];
        }else if(ip == N-4){
            RHSvec[ip] = c5*phi[ip-5] + c4*phi[ip-4] + c3*phi[ip-3] +
			 c2*phi[ip-2] + c1*phi[ip-1] + c0*phi[ip]   +
			 c1*phi[ip+1] + c2*phi[ip+2] + c3*phi[ip+3] +
			 c4*phi[0]    + c5*phi[1];
        }else if(ip == N-3){
            RHSvec[ip] = c5*phi[ip-5] + c4*phi[ip-4] + c3*phi[ip-3] +
			 c2*phi[ip-2] + c1*phi[ip-1] + c0*phi[ip]   +
			 c1*phi[ip+1] + c2*phi[ip+2] + c3*phi[0]    +
			 c4*phi[1]    + c5*phi[2];
        }else if(ip == N-2){
            RHSvec[ip] = c5*phi[ip-5] + c4*phi[ip-4] + c3*phi[ip-3] +
			 c2*phi[ip-2] + c1*phi[ip-1] + c0*phi[ip]   +
			 c1*phi[ip+1] + c2*phi[0]    + c3*phi[1]    +
			 c4*phi[2]    + c5*phi[3];
        }else if(ip == N-1){
            RHSvec[ip] = c5*phi[ip-5] + c4*phi[ip-4] + c3*phi[ip-3] +
			 c2*phi[ip-2] + c1*phi[ip-1] + c0*phi[ip]   +
			 c1*phi[0]    + c2*phi[1]    + c3*phi[2]    +
			 c4*phi[3] + c5*phi[4];
        }else{
            RHSvec[ip] = c5*phi[ip-5] + c4*phi[ip-4] + c3*phi[ip-3] +
			 c2*phi[ip-2] + c1*phi[ip-1] + c0*phi[ip]   +
			 c1*phi[ip+1] + c2*phi[ip+2] + c3*phi[ip+3] +
			 c4*phi[ip+4] + c5*phi[ip+5];
        }
    }

}

void FilterPeriodicY(double *phi, double *diagFy, double *offFy, double alphaF, double a0, double a1, double a2, double a3, double a4, double a5, int Nx, int Ny, double *phiF){

    double RHSvec[Ny];
    #pragma omp parallel for private(RHSvec)
    for(int ip = 0 ; ip < Nx; ip++){
	double *phiPointer = &phi[ip*Ny];
	multRHSFilter(a0, a1, a2, a3, a4, a5, phiPointer, Ny, RHSvec);
	cyclic(offFy, diagFy, offFy, alphaF, alphaF, RHSvec, Ny, &phiF[ip*Ny]);
    }
}

void FilterPeriodicX(double *phi, double *diagFx, double *offFx, double alphaF, double a0, double a1, double a2, double a3, double a4, double a5, int Nx, int Ny, double *phiF){

    double *phiTrans = new double[Nx*Ny];
    double *phiFTrans = new double[Nx*Ny];
    //transposeMatrix(phi, Nx, Ny, phiTrans);
    transposeMatrix_Fast2(phi, Nx, Ny, phiTrans, 60);

    double RHSvec[Nx];
    #pragma omp parallel for private(RHSvec) 
    for(int ip = 0 ; ip < Ny; ip++){
        double *phiPointer = &phiTrans[ip*Nx];
        multRHSFilter(a0, a1, a2, a3, a4, a5, phiPointer, Nx, RHSvec);
        cyclic(offFx, diagFx, offFx, alphaF, alphaF, RHSvec, Nx, &phiFTrans[ip*Nx]);
    }
    delete[] phiTrans;

    //transposeMatrix(phiFTrans, Ny, Nx, phiF);
    transposeMatrix_Fast2(phiFTrans, Ny, Nx, phiF, 60);
    delete[] phiFTrans;
}

void getRange(double *phi, std::string dataName, int Nx, int Ny){
    double dataMin = 1000000;
    double dataMax = -1000000;
    for(int ip = 0; ip < Nx*Ny; ip++){
	if(phi[ip] > dataMax){
	    dataMax = phi[ip];
	}

	if(phi[ip] < dataMin){
	    dataMin = phi[ip];
	}
    }

    cout << "Range of " << dataName << ": " << dataMin << ":" << dataMax << endl;

}

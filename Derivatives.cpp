#include <math.h>
#include <cstring>
#include <iostream>
#include "Utils.hpp"
#include "Derivatives.hpp"

using namespace std;

void Derivatives::multRHSDeriv(double dh, double *phi, int N, double *RHSvec){



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

void Derivatives::CompactDYPeriodic(double *phi, double *dphidy){

    double RHSvec[Ny];
//    #pragma omp parallel for private(RHSvec)
    for(int ip = 0 ; ip < Nx; ip++){
        double *phiPointer = &phi[ip*Ny];
        multRHSDeriv(dy, phiPointer, Ny, RHSvec);
        cyclic(offy, diagy, offy, alpha, alpha, RHSvec, Ny, &dphidy[ip*Ny]);
    }

}

void Derivatives::CompactDXPeriodic(double *phi, double *dphidx){

    double *phiTrans = new double[Nx*Ny];
    double *dphidxTrans = new double[Nx*Ny];
    //transposeMatrix(phi, Nx, Ny, phiTrans);
    transposeMatrix_Fast2(phi, Nx, Ny, phiTrans, 60);

    double RHSvec[Nx];
    #pragma omp parallel for private(RHSvec)
    for(int ip = 0 ; ip < Ny; ip++){
        double *phiPointer = &phiTrans[ip*Nx];
        multRHSDeriv(dx, phiPointer, Nx, RHSvec);
        cyclic(offx, diagx, offx, alpha, alpha, RHSvec, Nx, &dphidxTrans[ip*Nx]);
    }
    delete[] phiTrans;

    //transposeMatrix(dphidxTrans, Ny, Nx, dphidx);
    transposeMatrix_Fast2(dphidxTrans, Ny, Nx, dphidx, 60);
    delete[] dphidxTrans;
}

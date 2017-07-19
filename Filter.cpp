#include "Filter.hpp"

using namespace std;

void Filter::multRHSFilter(double *phi, int N, double *RHSvec){

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

void Filter::FilterPeriodicY(double *phi, double *phiF){

    double RHSvec[Ny];
//    #pragma omp parallel for private(RHSvec)
    for(int ip = 0 ; ip < Nx; ip++){
        double *phiPointer = &phi[ip*Ny];
        multRHSFilter(phiPointer, Ny, RHSvec);
        cyclic(offFy, diagFy, offFy, alphaF, alphaF, RHSvec, Ny, &phiF[ip*Ny]);
    }
}

void Filter::FilterPeriodicX(double *phi, double *phiF){

    double *phiTrans = new double[Nx*Ny];
    double *phiFTrans = new double[Nx*Ny];
    //transposeMatrix(phi, Nx, Ny, phiTrans);
    transposeMatrix_Fast2(phi, Nx, Ny, phiTrans, 60);

    double RHSvec[Nx];
    #pragma omp parallel for private(RHSvec)
    for(int ip = 0 ; ip < Ny; ip++){
        double *phiPointer = &phiTrans[ip*Nx];
        multRHSFilter(phiPointer, Nx, RHSvec);
        cyclic(offFx, diagFx, offFx, alphaF, alphaF, RHSvec, Nx, &phiFTrans[ip*Nx]);
    }
    delete[] phiTrans;

    //transposeMatrix(phiFTrans, Ny, Nx, phiF);
    transposeMatrix_Fast2(phiFTrans, Ny, Nx, phiF, 60);
    delete[] phiFTrans;
}

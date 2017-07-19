#include <math.h>
#include "IdealGas.hpp"

void IdealGas::solveMu(double *T, double *mu){
    for(int ip = 0; ip < Ny*Nx; ip++){
        mu[ip] = mu_ref*pow(T[ip]/T_ref, 0.76);
    }
}

void IdealGas::solveU(double *rho, double *rhoU, double *U){
    for(int ip = 0; ip < Ny*Nx; ip++){
        U[ip] = rhoU[ip]/rho[ip];
    }
}

void IdealGas::solvep(double *rho, double *rhoE, double *U, double *V, double *p){
    for(int ip = 0; ip < Ny*Nx; ip++){
        p[ip] = (gamma-1)*(rhoE[ip] - 0.5 * rho[ip]*(U[ip]*U[ip] + V[ip]*V[ip]));
    }
}

void IdealGas::solveT(double *rho, double *p, double *T){
    for(int ip = 0; ip < Ny*Nx; ip++){
        T[ip] = p[ip]/(rho[ip]*R_gas);
    }
}

void IdealGas::solveSOS(double *rho, double *p, double *SOS){
    for(int ip = 0; ip < Ny*Nx; ip++){
        SOS[ip] = sqrt(gamma*p[ip]/rho[ip]);
    }
}


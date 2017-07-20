#include "Solver.hpp"

void Solver::hellotest(){
	cout << "Testing Testing" << endl;
}

void Solver::applyInitialCondition(){
    for(int ip = 0; ip < Nx*Ny; ip++){
	U[ip]     = U0[ip];
	V[ip]     = V0[ip];
	rho1[ip]  = rho0[ip];
	p[ip]     = p0[ip];
	rhoU1[ip] = rho0[ip]*U0[ip];
	rhoV1[ip] = rho0[ip]*V0[ip];
	T[ip]     = p0[ip]/(rho0[ip]*idealGas->R_gas);
    }

    idealGas->solverhoE(rho0, p0, U0, V0, rhoE1);
    idealGas->solveMu(T, mu);
    idealGas->solveSOS(rho0, p0, SOS);

}

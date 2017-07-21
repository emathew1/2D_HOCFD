#include "Solver.hpp"

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

void Solver::computeDtFromCFL_advanceTime(){
    double dtTemp = 100000000;
    double UChar, VChar;

    for(int ip = 0; ip < Nx*Ny; ip++){

	UChar = fabs(U[ip]) + SOS[ip];
	VChar = fabs(V[ip]) + SOS[ip];

	double dtTemp2 = CFL/(UChar/dx + VChar/dy);

	if(dtTemp2 < dtTemp){
	    dtTemp = dtTemp2;
	}
    }

    dt = dtTemp;

    time += dt;

}

void Solver::computeVelocityTemperatureGradients(){

    derivatives->CompactDXPeriodic(U, Ux);
    derivatives->CompactDXPeriodic(V, Vx);
    derivatives->CompactDYPeriodic(U, Uy);
    derivatives->CompactDYPeriodic(V, Vy);

    derivatives->CompactDXPeriodic(T, Tx);
    derivatives->CompactDYPeriodic(T, Ty);

}

void Solver::computeContinuity(double *rhoU, double *rhoV){
   
    derivatives->CompactDXPeriodic(rhoU, rhsDxOut);
    derivatives->CompactDYPeriodic(rhoV, rhsDyOut);

    for(int ip = 0; ip < Nx*Ny; ip++){
	rhok[ip] = -dt*(rhsDxOut[ip] + rhsDyOut[ip]);
    }

}

void Solver::computeXMomentum(double *rhoU, double *rhoV){
    
    for(int jp = 0; jp < Nx*Ny; jp++){
        rhsDxIn[jp] = -(rhoU[jp]*U[jp] + p[jp] - 2.0*mu[jp]*Ux[jp] + (2.0/3.0)*mu[jp]*(Ux[jp] + Vy[jp]));
        rhsDyIn[jp] = -(rhoV[jp]*U[jp] - mu[jp]*(Uy[jp] + Vx[jp]));
    }

    derivatives->CompactDXPeriodic(rhsDxIn, rhsDxOut);
    derivatives->CompactDYPeriodic(rhsDyIn, rhsDyOut);
    
    for(int jp = 0; jp < Nx*Ny; jp++){
	rhoUk[jp] = dt*(rhsDxOut[jp]+rhsDyOut[jp]);
    }
}

void Solver::computeYMomentum(double *rhoU, double *rhoV){

    for(int jp = 0; jp < Nx*Ny; jp++){
        rhsDxIn[jp] = -(rhoU[jp]*V[jp] - mu[jp]*(Uy[jp] + Vx[jp]));
        rhsDyIn[jp] = -(rhoV[jp]*V[jp] + p[jp] - 2.0*mu[jp]*Vy[jp] + (2.0/3.0)*mu[jp]*(Ux[jp] + Vy[jp]));
    }

    derivatives->CompactDXPeriodic(rhsDxIn, rhsDxOut);
    derivatives->CompactDYPeriodic(rhsDyIn, rhsDyOut);

    for(int jp = 0; jp < Nx*Ny; jp++){
        rhoVk[jp] = dt*(rhsDxOut[jp]+rhsDyOut[jp]);
    }
}

void Solver::computeEnergy(double *rhoE){

    for(int jp = 0; jp < Nx*Ny; jp++){
        rhsDxIn[jp] = -(rhoE[jp]*U[jp] + U[jp]*p[jp] - (mu[jp]/idealGas->Pr/(idealGas->gamma-1.0))*Tx[jp] +
                            -U[jp]*(2.0*mu[jp]*Ux[jp] - (2.0/3.0)*mu[jp]*(Ux[jp]+Vy[jp])) +
                            -V[jp]*(mu[jp]*(Vx[jp] + Uy[jp])));
        rhsDyIn[jp] = -(rhoE[jp]*V[jp] + V[jp]*p[jp] - (mu[jp]/idealGas->Pr/(idealGas->gamma-1.0))*Ty[jp] +
                            -V[jp]*(2.0*mu[jp]*Vy[jp] - (2.0/3.0)*mu[jp]*(Ux[jp]+Vy[jp])) +
                            -U[jp]*(mu[jp]*(Vx[jp] + Uy[jp])));
    }

    derivatives->CompactDXPeriodic(rhsDxIn, rhsDxOut);
    derivatives->CompactDYPeriodic(rhsDyIn, rhsDyOut);

    for(int jp = 0; jp < Nx*Ny; jp++){
        rhoEk[jp] = dt*(rhsDxOut[jp]+rhsDyOut[jp]);
    }


}

void Solver::updateSolutionRKStep1(){

    //Update the final solution
    for(int jp = 0; jp < Nx*Ny; jp++){
        rho2[jp]  = rho1[jp]  + rhok[jp]/6.0;
        rhoU2[jp] = rhoU1[jp] + rhoUk[jp]/6.0;
        rhoV2[jp] = rhoV1[jp] + rhoVk[jp]/6.0;
        rhoE2[jp] = rhoE1[jp] + rhoEk[jp]/6.0;
    }

    //Get the intermediate solution
    for(int jp = 0; jp < Nx*Ny; jp++){
	rho1k[jp]  = rho1[jp] + rhok[jp]/2.0;
        rhoU1k[jp] = rhoU1[jp] + rhoUk[jp]/2.0;
        rhoV1k[jp] = rhoV1[jp] + rhoVk[jp]/2.0;
        rhoE1k[jp] = rhoE1[jp] + rhoEk[jp]/2.0;
    }

    //update primative and fluid properties
    idealGas->solveU(rho1k, rhoU1k, U);
    idealGas->solveU(rho1k, rhoV1k, V);
    idealGas->solvep(rho1k, rhoE1k, U, V, p);
    idealGas->solveT(rho1k, p, T);
    idealGas->solveMu(T, mu);

}

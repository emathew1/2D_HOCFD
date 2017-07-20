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

void Solver::computeContinuity(){
   
    derivatives->CompactDXPeriodic(rhoU1, rhsDxOut);
    derivatives->CompactDYPeriodic(rhoV1, rhsDyOut);

    for(int ip = 0; ip < Nx*Ny; ip++){
	rhok[ip] = -dt*(rhsDxOut[ip] + rhsDyOut[ip]);
    }

}

void Solver::computeXMomentum(){
    
    for(int jp = 0; jp < Nx*Ny; jp++){
        rhsDxIn[jp] = -(rhoU1[jp]*U[jp] + p[jp] - 2.0*mu[jp]*Ux[jp] + (2.0/3.0)*mu[jp]*(Ux[jp] + Vy[jp]));
        rhsDyIn[jp] = -(rhoV1[jp]*U[jp] - mu[jp]*(Uy[jp] + Vx[jp]));
    }

    derivatives->CompactDXPeriodic(rhsDxIn, rhsDxOut);
    derivatives->CompactDYPeriodic(rhsDyIn, rhsDyOut);
    
    for(int jp = 0; jp < Nx*Ny; jp++){
	rhoUk[jp] = dt*(rhsDxOut[jp]+rhsDyOut[jp]);
    }
}

void Solver::computeYMomentum(){

    for(int jp = 0; jp < Nx*Ny; jp++){
        rhsDxIn[jp] = -(rhoU1[jp]*V[jp] - mu[jp]*(Uy[jp] + Vx[jp]));
        rhsDyIn[jp] = -(rhoV1[jp]*V[jp] + p[jp] - 2.0*mu[jp]*Vy[jp] + (2.0/3.0)*mu[jp]*(Ux[jp] + Vy[jp]));
    }

    derivatives->CompactDXPeriodic(rhsDxIn, rhsDxOut);
    derivatives->CompactDYPeriodic(rhsDyIn, rhsDyOut);

    for(int jp = 0; jp < Nx*Ny; jp++){
        rhoVk[jp] = dt*(rhsDxOut[jp]+rhsDyOut[jp]);
    }
}

void Solver::computeEnergy(){

    for(int jp = 0; jp < Nx*Ny; jp++){
        rhsDxIn[jp] = -(rhoE1[jp]*U[jp] + U[jp]*p[jp] - (mu[jp]/idealGas->Pr/(idealGas->gamma-1.0))*Tx[jp] +
                            -U[jp]*(2.0*mu[jp]*Ux[jp] - (2.0/3.0)*mu[jp]*(Ux[jp]+Vy[jp])) +
                            -V[jp]*(mu[jp]*(Vx[jp] + Uy[jp])));
        rhsDyIn[jp] = -(rhoE1[jp]*V[jp] + V[jp]*p[jp] - (mu[jp]/idealGas->Pr/(idealGas->gamma-1.0))*Ty[jp] +
                            -V[jp]*(2.0*mu[jp]*Vy[jp] - (2.0/3.0)*mu[jp]*(Ux[jp]+Vy[jp])) +
                            -U[jp]*(mu[jp]*(Vx[jp] + Uy[jp])));
    }

    derivatives->CompactDXPeriodic(rhsDxIn, rhsDxOut);
    derivatives->CompactDYPeriodic(rhsDyIn, rhsDyOut);

    for(int jp = 0; jp < Nx*Ny; jp++){
        rhoEk[jp] = dt*(rhsDxOut[jp]+rhsDyOut[jp]);
    }


}

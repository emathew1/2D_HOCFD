#include "Solver.hpp"

void Solver::setBCForDerivatives(){

    if(bcX0 == SPONGE || bcX0 ==  WALL || bcX0 ==  MOVING_WALL){
	derivatives->offx2[0] = derivatives->alpha1;
	derivatives->offx2[1] = derivatives->alpha22;
	derivatives->offx1[1] = derivatives->alpha21;
	filter->alphaFx[0] = 0.0;
	filter->offFx2[0] = 0.0;
    }

    if(bcX1 == SPONGE || bcX1 ==  WALL || bcX1 ==  MOVING_WALL){
	derivatives->offx2[Nx-2] = derivatives->alpha21;
	derivatives->offx1[Nx-2] = derivatives->alpha22;
	derivatives->offx1[Nx-1] = derivatives->alpha1;
	filter->alphaFx[Nx-1] = 0.0;
	filter->offFx1[Nx-1] = 0.0;
    }

    if(bcY0 == SPONGE || bcY0 ==  WALL || bcY0 ==  MOVING_WALL){
	derivatives->offy2[0] = derivatives->alpha1;
	derivatives->offy2[1] = derivatives->alpha22;
	derivatives->offy1[1] = derivatives->alpha21;
	filter->alphaFy[0] = 0.0;
	filter->offFy2[0] = 0.0;
    }

    if(bcY1 == SPONGE || bcY1 ==  WALL || bcY1 ==  MOVING_WALL){
	derivatives->offy2[Ny-2] = derivatives->alpha21;
	derivatives->offy1[Ny-2] = derivatives->alpha22;
	derivatives->offy1[Ny-1] = derivatives->alpha1;
	filter->alphaFy[Ny-1] = 0.0;
	filter->offFy1[Ny-1] = 0.0;
    }

    if(bcY0 == SPONGE || bcY1 == SPONGE || bcX0 == SPONGE || bcX1 == SPONGE){ 
        //Sponge Stuff
        spongeAvgT = 1.0;
        spongeEpsP = 0.005;
        spongeAvgRho     = new double[Nx*Ny];
	spongeRhoSource  = new double[Nx*Ny];
        spongeAvgRhoU    = new double[Nx*Ny];
	spongeRhoUSource = new double[Nx*Ny];
        spongeAvgRhoV    = new double[Nx*Ny];
	spongeRhoVSource = new double[Nx*Ny];
        spongeAvgRhoE    = new double[Nx*Ny];
	spongeRhoESource = new double[Nx*Ny];
        spongeSigma 	 = new double[Nx*Ny];
        spongeStrength = 12.0;
        spongeLengthX0 = 0;   
        spongeLengthX1 = 0;   
        spongeLengthY0 = 0;   
        spongeLengthY1 = 0;   
    }

}

void Solver::initSpongeStuff(){


    //Initialize the sponge distribution at the edges of the domains with sponges
    for(int ip = 0; ip < Nx; ip++){
	for(int jp = 0; jp < Ny; jp++){

	    //Initialize sponge strength as zero everywhere
	    spongeSigma[ip*Ny + jp] = 0.0;	


	    if(bcX0 == SPONGE){
		if(x[ip] < spongeLengthX0 && spongeLengthX0 != 0){
		    double spongeX = (spongeLengthX0-x[ip])/spongeLengthX0;
		    spongeSigma[ip*Ny + jp] = spongeStrength*(0.068*pow(spongeX, 2.0) + 0.845*pow(spongeX, 8.0));
		}
	    }

	    if(bcX1 == SPONGE){
		if(x[ip] > (Lx - spongeLengthX1) && spongeLengthX1 != 0){
		    double spongeX = (x[ip] - (Lx - spongeLengthX1))/spongeLengthX1;
		    spongeSigma[ip*Ny + jp] = spongeStrength*(0.068*pow(spongeX, 2.0) + 0.845*pow(spongeX, 8.0));
		}
	    }

	    if(bcY0 == SPONGE){
		if(y[jp] < spongeLengthY0 && spongeLengthY0 != 0){
		    double spongeY = (spongeLengthY0-y[jp])/spongeLengthY0;
		    spongeSigma[ip*Ny + jp] = spongeStrength*(0.068*pow(spongeY, 2.0) + 0.845*pow(spongeY, 8.0));
		}
	    }

	    if(bcY1 == SPONGE){
		if(y[jp] > (Ly - spongeLengthY1) && spongeLengthY1 != 0){
		    double spongeY = (y[jp] - (Ly - spongeLengthY1))/spongeLengthY1;
		    spongeSigma[ip*Ny + jp] = spongeStrength*(0.068*pow(spongeY, 2.0) + 0.845*pow(spongeY, 8.0));
		}
	    } 
	}
    }   

   //Initialize the sponge average as the initial condition
   for(int ip = 0; ip < Nx*Ny; ip++){
       spongeAvgRho[ip]  = rho1[ip];   
       spongeAvgRhoU[ip] = rhoU1[ip];   
       spongeAvgRhoV[ip] = rhoV1[ip];   
       spongeAvgRhoE[ip] = rhoE1[ip];   
   }

}

void Solver::computeSpongeSource(double *rhoIn, double *rhoUIn, double *rhoVIn, double *rhoEIn){


   if(bcX0 == SPONGE || bcX1 == SPONGE || bcY0 == SPONGE || bcY1 == SPONGE){
       cout << "TEST" << endl;
       for(int ip = 0; ip < Nx*Ny; ip++){
           spongeRhoSource[ip]  = spongeSigma[ip]*(spongeAvgRho[ip]  - rhoIn[ip]);
           spongeRhoUSource[ip] = spongeSigma[ip]*(spongeAvgRhoU[ip] - rhoUIn[ip]);
           spongeRhoVSource[ip] = spongeSigma[ip]*(spongeAvgRhoV[ip] - rhoVIn[ip]);
           spongeRhoESource[ip] = spongeSigma[ip]*(spongeAvgRhoE[ip] - rhoEIn[ip]);
       }
   }
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
//    dt = 0.008;

    time += dt;
    step++;
}

void Solver::computeCompactDY(double *phi, double *dphidy){

    if(bcY0 == PERIODIC || bcY1 == PERIODIC){
        derivatives->CompactDYPeriodic(phi, dphidy);
    }else{
	derivatives->CompactDYDirichlet(phi, dphidy);
    }

}

void Solver::computeCompactDX(double *phi, double *dphidx){

    if(bcX0 == PERIODIC || bcX1 == PERIODIC){
        derivatives->CompactDXPeriodic(phi, dphidx);
    }else{
	derivatives->CompactDXDirichlet(phi, dphidx);
    }
}

void Solver::computeVelocityTemperatureGradients(){


        computeCompactDY(U, Uy);
        computeCompactDY(V, Vy);
        computeCompactDY(T, Ty);

        computeCompactDX(U, Ux);
        computeCompactDX(V, Vx);
        computeCompactDX(T, Tx);

}

void Solver::computeContinuity(double *rhoU, double *rhoV){

    computeCompactDX(rhoU, rhsDxOut);
    computeCompactDY(rhoV, rhsDyOut);

    for(int ip = 0; ip < Nx*Ny; ip++){
	rhok[ip] = dt*(-rhsDxOut[ip] -rhsDyOut[ip]);
    }

}

void Solver::computeXMomentum(double *rhoU, double *rhoV){
    
    for(int jp = 0; jp < Nx*Ny; jp++){
        rhsDxIn[jp] = -(rhoU[jp]*U[jp] + p[jp] - 2.0*mu[jp]*Ux[jp] + (2.0/3.0)*mu[jp]*(Ux[jp] + Vy[jp]));
        rhsDyIn[jp] = -(rhoV[jp]*U[jp] - mu[jp]*(Uy[jp] + Vx[jp]));
    }

    computeCompactDX(rhsDxIn, rhsDxOut);
    computeCompactDY(rhsDyIn, rhsDyOut);
    
    for(int jp = 0; jp < Nx*Ny; jp++){
	rhoUk[jp] = dt*(rhsDxOut[jp]+rhsDyOut[jp]);
    }
}

void Solver::computeYMomentum(double *rhoU, double *rhoV){

    for(int jp = 0; jp < Nx*Ny; jp++){
        rhsDxIn[jp] = -(rhoU[jp]*V[jp] - mu[jp]*(Uy[jp] + Vx[jp]));
        rhsDyIn[jp] = -(rhoV[jp]*V[jp] + p[jp] - 2.0*mu[jp]*Vy[jp] + (2.0/3.0)*mu[jp]*(Ux[jp] + Vy[jp]));
    }

    computeCompactDX(rhsDxIn, rhsDxOut);
    computeCompactDY(rhsDyIn, rhsDyOut);

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

    computeCompactDX(rhsDxIn, rhsDxOut);
    computeCompactDY(rhsDyIn, rhsDyOut);

    for(int jp = 0; jp < Nx*Ny; jp++){
        rhoEk[jp] = dt*(rhsDxOut[jp]+rhsDyOut[jp]);
    }


}

void Solver::computeRhs(){

    if(bcX0 == SPONGE || bcX1 == SPONGE || bcY0 == SPONGE || bcY1 == SPONGE){
        for(int jp = 0; jp < Nx*Ny; jp++){
	    rhok[jp]  += dt*spongeRhoSource[jp];	
	    rhoUk[jp] += dt*spongeRhoUSource[jp];	
	    rhoVk[jp] += dt*spongeRhoVSource[jp];	
	    rhoEk[jp] += dt*spongeRhoESource[jp];	
        }
    }

}

void Solver::enforceBCs(){

    if(bcX0 != PERIODIC){
        int ip = 0;
        for(int jp = 0; jp < Ny; jp++){
	    rhok[ip*Ny + jp] = 0.0;
	    rhoUk[ip*Ny + jp] = 0.0;
	    rhoVk[ip*Ny + jp] = 0.0;
	    rhoEk[ip*Ny + jp] = 0.0;
	}
    }

    if(bcX1 != PERIODIC){
        int ip = Nx-1;
        for(int jp = 0; jp < Ny; jp++){
	    rhok[ip*Ny + jp] = 0.0;
	    rhoUk[ip*Ny + jp] = 0.0;
	    rhoVk[ip*Ny + jp] = 0.0;
	    rhoEk[ip*Ny + jp] = 0.0;
	}
    }

    if(bcY0 != PERIODIC){
        for(int ip = 0; ip < Nx; ip++){
            int jp = 0; 
	    rhok[ip*Ny + jp] = 0.0;
	    rhoUk[ip*Ny + jp] = 0.0;
	    rhoVk[ip*Ny + jp] = 0.0;
	    rhoEk[ip*Ny + jp] = 0.0;
	}
    }


    if(bcY1 != PERIODIC){
        for(int ip = 0; ip < Nx; ip++){
            int jp = Ny-1; 
	    rhok[ip*Ny + jp] = 0.0;
	    rhoUk[ip*Ny + jp] = 0.0;
	    rhoVk[ip*Ny + jp] = 0.0;
	    rhoEk[ip*Ny + jp] = 0.0;
	}
    }

}

void Solver::updateSolutionRKStep1(){

    
    //Don't want to disturb the BC's if we're non periodic
    enforceBCs();

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
    idealGas->solveSOS(rho1k, p, SOS);

}

void Solver::updateSolutionRKStep2(){

    //Don't want to disturb the BC's if we're non periodic
    enforceBCs();

    //Update the final solution
    for(int jp = 0; jp < Nx*Ny; jp++){
        rho2[jp]  += rhok[jp]/3.0;
        rhoU2[jp] += rhoUk[jp]/3.0;
        rhoV2[jp] += rhoVk[jp]/3.0;
        rhoE2[jp] += rhoEk[jp]/3.0;
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
    idealGas->solveSOS(rho1k, p, SOS);
}

void Solver::updateSolutionRKStep3(){

    //Don't want to disturb the BC's if we're non periodic
    enforceBCs();

    //Update the final solution
    for(int jp = 0; jp < Nx*Ny; jp++){
        rho2[jp]  += rhok[jp]/3.0;
        rhoU2[jp] += rhoUk[jp]/3.0;
        rhoV2[jp] += rhoVk[jp]/3.0;
        rhoE2[jp] += rhoEk[jp]/3.0;
    }

    //Get the intermediate solution
    for(int jp = 0; jp < Nx*Ny; jp++){
	rho1k[jp]  = rho1[jp] + rhok[jp];
        rhoU1k[jp] = rhoU1[jp] + rhoUk[jp];
        rhoV1k[jp] = rhoV1[jp] + rhoVk[jp];
        rhoE1k[jp] = rhoE1[jp] + rhoEk[jp];
    }

    //update primative and fluid properties
    idealGas->solveU(rho1k, rhoU1k, U);
    idealGas->solveU(rho1k, rhoV1k, V);
    idealGas->solvep(rho1k, rhoE1k, U, V, p);
    idealGas->solveT(rho1k, p, T);
    idealGas->solveMu(T, mu);
    idealGas->solveSOS(rho1k, p, SOS);

}

void Solver::updateSolutionRKStep4(){

    //Don't want to disturb the BC's if we're non periodic
    enforceBCs();

    //Update the final solution
    for(int jp = 0; jp < Nx*Ny; jp++){
        rho2[jp]  += rhok[jp]/6.0;
        rhoU2[jp] += rhoUk[jp]/6.0;
        rhoV2[jp] += rhoVk[jp]/6.0;
        rhoE2[jp] += rhoEk[jp]/6.0;
    }

}

void Solver::filterCompactY(double *phi, double *phiF){

    if(bcY0 == PERIODIC || bcY1 == PERIODIC){
        filter->FilterPeriodicY(phi, phiF);
    }else{
	filter->FilterFiniteDomainY(phi, phiF);
    }

}

void Solver::filterCompactX(double *phi, double *phiF){

    if(bcX0 == PERIODIC || bcX1 == PERIODIC){
        filter->FilterPeriodicX(phi, phiF);
    }else{
	filter->FilterFiniteDomainX(phi, phiF);
    }
}

void Solver::filterAndUpdateSolution(){

    if(step%filterStep == 0){
	if(filterCount%2 == 0){
	    double *rhoTemp  = new double[Nx*Ny];
	    double *rhoUTemp = new double[Nx*Ny];
	    double *rhoVTemp = new double[Nx*Ny];
	    double *rhoETemp = new double[Nx*Ny];

            filterCompactX(rho2,  rhoTemp);
            filterCompactX(rhoU2, rhoUTemp);
            filterCompactX(rhoV2, rhoVTemp);
            filterCompactX(rhoE2, rhoETemp);

            filterCompactY(rhoTemp,  rho1);
            filterCompactY(rhoUTemp, rhoU1);
            filterCompactY(rhoVTemp, rhoV1);
            filterCompactY(rhoETemp, rhoE1);

	    delete[] rhoTemp;
	    delete[] rhoUTemp;
	    delete[] rhoVTemp;
	    delete[] rhoETemp;
	}else{

            double *rhoTemp = new double[Nx*Ny];
            double *rhoUTemp = new double[Nx*Ny];
            double *rhoVTemp = new double[Nx*Ny];
            double *rhoETemp = new double[Nx*Ny];
         
	    filterCompactY(rho2,  rhoTemp);
            filterCompactY(rhoU2, rhoUTemp);
            filterCompactY(rhoV2, rhoVTemp);
            filterCompactY(rhoE2, rhoETemp);

            filterCompactX(rhoTemp,  rho1);
            filterCompactX(rhoUTemp, rhoU1);
            filterCompactX(rhoVTemp, rhoV1);
            filterCompactX(rhoETemp, rhoE1);

            delete[] rhoTemp;
            delete[] rhoUTemp;
            delete[] rhoVTemp;
            delete[] rhoETemp;
	}
	filterCount++;	
    }else{

        //Just Update the solutions
        memcpy(rho1,  rho2,  Nx*Ny*sizeof(double));   
        memcpy(rhoU1, rhoU2, Nx*Ny*sizeof(double)); 
        memcpy(rhoV1, rhoV2, Nx*Ny*sizeof(double)); 
        memcpy(rhoE1, rhoE2, Nx*Ny*sizeof(double)); 

    }

}

void Solver::updateEndOfStepPrimAndTemp(){

    idealGas->solveU(rho1, rhoU1, U);
    idealGas->solveU(rho1, rhoV1, V);
    idealGas->solvep(rho1, rhoE1, U, V, p);
    idealGas->solveT(rho1, p, T);
    idealGas->solveMu(T, mu);
    idealGas->solveSOS(rho1, p, SOS);

}

void Solver::updateSpongeBCs(){

    if(bcX0 == SPONGE || bcX1 == SPONGE || bcY0 == SPONGE || bcY1 == SPONGE){

        double eps = 1.0/(spongeAvgT/dt + 1.0);
   
    	for(int ip = 0; ip < Nx*Ny; ip++){
            spongeAvgRho[ip]  += eps*(rho1[ip]  - spongeAvgRho[ip]);
            spongeAvgRhoU[ip] += eps*(rhoU1[ip] - spongeAvgRhoU[ip]);
            spongeAvgRhoV[ip] += eps*(rhoV1[ip] - spongeAvgRhoV[ip]);
            spongeAvgRhoE[ip] += eps*(rhoE1[ip] - spongeAvgRhoE[ip]);

	    spongeAvgRhoE[ip] = spongeEpsP*spongeAvgRhoE[ip] + 
	        (1.0-spongeEpsP)*(spongeP/(1.0-idealGas->gamma) + 0.5*(pow(spongeAvgRhoU[ip],2.0) + 
	        pow(spongeAvgRhoV[ip],2.0))/spongeAvgRho[ip]);
        }
        //Need to set boundary conditions on the sponge for the next step

        if(bcX0 == SPONGE){
	    for(int ip = 0; ip < Ny; ip++){
	    	int ii = ip*Nx;
	    	rho1[ii]  = spongeAvgRho[ii];
	    	rhoU1[ii] = spongeAvgRhoU[ii];
	    	rhoV1[ii] = spongeAvgRhoV[ii];
	    	rhoE1[ii] = spongeAvgRhoE[ii];
	    }
    	}

    	if(bcX1 == SPONGE){
	    for(int ip = 0; ip < Ny; ip++){
	    	int ii = ip*Nx + (Ny-1);
	    	rho1[ii]  = spongeAvgRho[ii];
	    	rhoU1[ii] = spongeAvgRhoU[ii];
	    	rhoV1[ii] = spongeAvgRhoV[ii];
	    	rhoE1[ii] = spongeAvgRhoE[ii];
	    }
    	}

    	if(bcY0 == SPONGE){
	    for(int ip = 0; ip < Nx; ip++){
	    	int ii = ip;
	    	rho1[ii]  = spongeAvgRho[ii];
	    	rhoU1[ii] = spongeAvgRhoU[ii];
	    	rhoV1[ii] = spongeAvgRhoV[ii];
	    	rhoE1[ii] = spongeAvgRhoE[ii];
	    }
    	}

        if(bcY1 == SPONGE){
	    for(int ip = 0; ip < Nx; ip++){
	    	int ii = (Ny-1)*Nx + ip;
	    	rho1[ii]  = spongeAvgRho[ii];
	    	rhoU1[ii] = spongeAvgRhoU[ii];
	    	rhoV1[ii] = spongeAvgRhoV[ii];
	        rhoE1[ii] = spongeAvgRhoE[ii];
	    }
        }

    }
}

void Solver::checkSolution(){

    if(step%checkStep == 0){
        t2 = std::chrono::system_clock::now();
        cout << "-------------------------------------------------" << endl;
        cout << " Step = "<< step << ", time = " << time << ", dt = " << dt << endl; 
        cout << "-------------------------------------------------" << endl;
        cout << "  Time since last timestep = " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()/(double)1000000000 << endl;
        getRange(rho1, "RHO", Nx, Ny);
        getRange(U, "U", Nx, Ny);
        getRange(V, "V", Nx, Ny);
        getRange(p, "P", Nx, Ny);
        getRange(T, "T", Nx, Ny);
        getRange(mu, "mu", Nx, Ny);
        getRange(rhoE1, "RHOE", Nx, Ny);
        getRange(rhoU1, "RHOU", Nx, Ny);
        getRange(SOS, "SOS", Nx, Ny);
        cout << endl;

        t1 = std::chrono::system_clock::now();
    }
}

void Solver::dumpSolution(){

        cout << endl;
        cout << "===============" << endl;
        cout << " DUMPING FIELD " << endl;
        cout << "===============" << endl;


        ofstream outfile;
        outfile.precision(17);
        string outputFileName;
        outputFileName = "rho.out.";
        outputFileName.append(to_string(step)); 
        outfile.open(outputFileName);
        for(int jp = 0; jp < Nx; jp++){
            for(int kp = 0; kp < Ny; kp++){
                outfile << rho1[jp*Ny+kp] << " ";
            }
            outfile << endl;
        }
        outfile.close();

        outputFileName = "rhoU.out."; 
        outputFileName.append(to_string(step)); 
        outfile.open(outputFileName);
        outfile.precision(17);
        for(int jp = 0; jp < Nx; jp++){
            for(int kp = 0; kp < Ny; kp++){
                outfile << rhoU1[jp*Ny+kp] << " ";
            }
            outfile << endl;
        }
        outfile.close();

        outputFileName = "rhoV.out."; 
        outputFileName.append(to_string(step)); 
        outfile.open(outputFileName);
        outfile.precision(17);
        for(int jp = 0; jp < Nx; jp++){
            for(int kp = 0; kp < Ny; kp++){
                outfile << rhoV1[jp*Ny+kp] << " ";
            }
            outfile << endl;
        }
        outfile.close();

        outputFileName = "rhoE.out."; 
        outputFileName.append(to_string(step)); 
        outfile.open(outputFileName);
        outfile.precision(17);
        for(int jp = 0; jp < Nx; jp++){
            for(int kp = 0; kp < Ny; kp++){
                outfile << rhoE1[jp*Ny+kp]*SOS[jp*Ny+kp] << " ";
            }
            outfile << endl;
        }
        outfile.close();

}

void Solver::checkEnd(){

    if(time >= timeEnd){
	cout << "=================" << endl;
	cout << " HIT END OF TIME " << endl;
	cout << "=================" << endl;

	endFlag = 1;
    }

    if(step >= timeStep){
	cout << "=================" << endl;
	cout << " HIT END OF TIME " << endl;
	cout << "=================" << endl;

	endFlag = 1;

    } 

    if(endFlag){
	dumpSolution();
    }

}


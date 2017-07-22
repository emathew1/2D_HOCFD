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
    dt = 0.008;

    time += dt;
    step++;
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
    idealGas->solveSOS(rho1k, p, SOS);

}

void Solver::updateSolutionRKStep2(){

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

    //Update the final solution
    for(int jp = 0; jp < Nx*Ny; jp++){
        rho2[jp]  += rhok[jp]/6.0;
        rhoU2[jp] += rhoUk[jp]/6.0;
        rhoV2[jp] += rhoVk[jp]/6.0;
        rhoE2[jp] += rhoEk[jp]/6.0;
    }

}

void Solver::filterAndUpdateSolution(){

    if(step%filterStep == 0){
	if(filterCount%2 == 0){
	    double *rhoTemp  = new double[Nx*Ny];
	    double *rhoUTemp = new double[Nx*Ny];
	    double *rhoVTemp = new double[Nx*Ny];
	    double *rhoETemp = new double[Nx*Ny];

            filter->FilterPeriodicX(rho2,  rhoTemp);
            filter->FilterPeriodicX(rhoU2, rhoUTemp);
            filter->FilterPeriodicX(rhoV2, rhoVTemp);
            filter->FilterPeriodicX(rhoE2, rhoETemp);

            filter->FilterPeriodicY(rhoTemp,  rho1);
            filter->FilterPeriodicY(rhoUTemp, rhoU1);
            filter->FilterPeriodicY(rhoVTemp, rhoV1);
            filter->FilterPeriodicY(rhoETemp, rhoE1);

	    delete[] rhoTemp;
	    delete[] rhoUTemp;
	    delete[] rhoVTemp;
	    delete[] rhoETemp;
	}else{

            double *rhoTemp = new double[Nx*Ny];
            double *rhoUTemp = new double[Nx*Ny];
            double *rhoVTemp = new double[Nx*Ny];
            double *rhoETemp = new double[Nx*Ny];
         
	    filter->FilterPeriodicY(rho2,  rhoTemp);
            filter->FilterPeriodicY(rhoU2, rhoUTemp);
            filter->FilterPeriodicY(rhoV2, rhoVTemp);
            filter->FilterPeriodicY(rhoE2, rhoETemp);

            filter->FilterPeriodicX(rhoTemp,  rho1);
            filter->FilterPeriodicX(rhoUTemp, rhoU1);
            filter->FilterPeriodicX(rhoVTemp, rhoV1);
            filter->FilterPeriodicX(rhoETemp, rhoE1);

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

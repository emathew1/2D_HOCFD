#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <chrono>

#include "Utils.hpp"
#include "Solver.hpp"


using namespace std;
using namespace std::chrono;


int main(int argc, char *argv[]){


    cout << endl;
    cout << "--------------------------------------" << endl;
    cout << " 2-D High Order Solver " << endl;
    cout << "--------------------------------------" << endl;
    cout << endl;

    int Nx = 200; double Lx = 4.0;//2*M_PI - 2*M_PI/((double)Nx);
    int Ny = 200; double Ly = 4.0;//2*M_PI - 2*M_PI/((double)Ny);

    Solver *solver = new Solver(Nx, Ny, Lx, Ly);

    //Set some of the solver parameters
    solver->timeStep    = 100000;
    solver->timeEnd     = 0.25;
    solver->filterStep  = 1;
    solver->checkStep   = 1;
    solver->outputStep  = 100;
    solver->CFL		= 0.25;

    //Change fluid properties
    solver->idealGas->mu_ref = 0.0000001;

    //Set the BC's for the solver
    solver->bcX0 = Solver::PERIODIC;
    solver->bcX1 = Solver::PERIODIC;
    solver->bcY0 = Solver::PERIODIC;
    solver->bcY1 = Solver::PERIODIC;
    solver->setBCForDerivatives();

    //Need to set alpha = 0 at sponge ends

    //Allocate an initial condition
    for(int ip = 0; ip < Nx; ip++){
	for(int jp = 0; jp < Ny; jp++){
	    solver->U0[ip*Ny + jp] = 0.0;
	    solver->V0[ip*Ny + jp] = 0.0;
	    solver->p0[ip*Ny + jp] = 1/solver->idealGas->gamma + 
			10.0*exp(-(pow(solver->x[ip]-2.0,2.0) + pow(solver->y[jp] - 2.0,2.0))/0.25);
	    solver->rho0[ip*Ny + jp] = 1.0;

/*
	    if(solver->x[ip] < 2.0){
	        solver->p0[ip*Ny + jp]   = 1.0/solver->idealGas->gamma;
	        solver->U0[ip*Ny + jp]   = 0.0;
	        solver->V0[ip*Ny + jp]   = 0.0;
	        solver->rho0[ip*Ny + jp] = 1.0;
	    }else{
	        solver->p0[ip*Ny + jp]   = 0.1/solver->idealGas->gamma;
	        solver->U0[ip*Ny + jp]   = 0.0;
	        solver->V0[ip*Ny + jp]   = 0.0;
	        solver->rho0[ip*Ny + jp] = 0.125;
	    }
*/
	}
    }

    solver->applyInitialCondition();

/*
    //Provide non-standard attributes for the sponge layers
    solver->spongeAvgT = 10.0;
    solver->spongeLengthX0 = 0.75;
    solver->spongeLengthY0 = 0.75;
    solver->spongeLengthX1 = 0.75;
    solver->spongeLengthY1 = 0.75;
    solver->spongeP = 1.0/solver->idealGas->gamma;
    solver->initSpongeStuff();
*/

/*
        cout.precision(17);


    for(int ip = 0; ip < Ny; ip++){
        cout <<  solver->filter->alphaFy[ip] << endl;;
    }
    cout << endl;

    solver->filter->calcFilterCoefficients(solver->filter->alphaFx);
    double RHSvec[solver->Nx];
    double *work = new double[solver->Nx];

    double *test = new double[solver->Nx];
    for(int ip = 0; ip < Nx; ip++){
	if(ip < (double)Nx/2.0-0.5){
	    test[ip] = (double)ip;
	}else{
	    test[ip] = 0.0;
	}
    }    

    for(int ip = 0; ip < Nx; ip++){
	cout << test[ip] << endl;
    } 
    cout << endl;

    //solver->filter->multRHSFilterFiniteDomain(test, Nx, RHSvec);
    solver->filter->multRHSFilter(test, Nx, RHSvec);

    for(int ip = 0; ip < Nx; ip++){
	cout <<  RHSvec[ip] << endl;
    }
    cout << endl;


    double *filterTest = new double[solver->Nx];
    //solveTri(solver->filter->offFx, solver->filter->diagFx, solver->filter->offFx, RHSvec, filterTest, work, Nx);
    cyclic(solver->filter->offFy, solver->filter->diagFy, solver->filter->offFy, 0.49,0.49,  RHSvec, Ny, filterTest);



    for(int ip = 0; ip < Nx; ip++){
	cout << filterTest[ip] << endl;
    }
    double *rhoTemp  = new double[Nx*Ny];

    solver->filterCompactX(solver->rho1, rhoTemp);
    memcpy(solver->rho1, rhoTemp, Nx*Ny*(sizeof(double)));
    solver->computeVelocityTemperatureGradients();
    solver->computeXMomentum(solver->rhoU1, solver->rhoV1);

    memcpy(solver->rhoU1, solver->rhoUk, Nx*Ny*(sizeof(double))); 
*/ 
    solver->dumpSolution();

    while(!solver->endFlag){

        //At the beginning of every step
        solver->computeDtFromCFL_advanceTime();


	//=======================
        // RK Step 1
	//=======================

        solver->computeVelocityTemperatureGradients();

	solver->computeSpongeSource(solver->rho1, solver->rhoU1, solver->rhoV1, solver->rhoE1);	


        solver->computeContinuity(solver->rhoU1, solver->rhoV1);

        solver->computeXMomentum(solver->rhoU1, solver->rhoV1);

        solver->computeYMomentum(solver->rhoU1, solver->rhoV1);

        solver->computeEnergy(solver->rhoE1);

	solver->computeRhs();
    
        solver->updateSolutionRKStep1();


	//=======================
        // RK Step 2
	//=======================

        solver->computeVelocityTemperatureGradients();
	
	solver->computeSpongeSource(solver->rho1k, solver->rhoU1k, solver->rhoV1k, solver->rhoE1k);	
        solver->computeContinuity(solver->rhoU1k, solver->rhoV1k);

        solver->computeXMomentum(solver->rhoU1k, solver->rhoV1k);

        solver->computeYMomentum(solver->rhoU1k, solver->rhoV1k);

        solver->computeEnergy(solver->rhoE1k);
 
	solver->computeRhs();

        solver->updateSolutionRKStep2();



	//=======================
        // RK Step 3
	//=======================

        solver->computeVelocityTemperatureGradients();
	
	solver->computeSpongeSource(solver->rho1k, solver->rhoU1k, solver->rhoV1k, solver->rhoE1k);	
        solver->computeContinuity(solver->rhoU1k, solver->rhoV1k);

        solver->computeXMomentum(solver->rhoU1k, solver->rhoV1k);

        solver->computeYMomentum(solver->rhoU1k, solver->rhoV1k);

        solver->computeEnergy(solver->rhoE1k);

	solver->computeRhs();

        solver->updateSolutionRKStep3();


	//=======================
        // RK Step 4
	//=======================

        solver->computeVelocityTemperatureGradients();
	
	solver->computeSpongeSource(solver->rho1k, solver->rhoU1k, solver->rhoV1k, solver->rhoE1k);	
        solver->computeContinuity(solver->rhoU1k, solver->rhoV1k);

        solver->computeXMomentum(solver->rhoU1k, solver->rhoV1k);

        solver->computeYMomentum(solver->rhoU1k, solver->rhoV1k);

        solver->computeEnergy(solver->rhoE1k);

	solver->computeRhs();

        solver->updateSolutionRKStep4();
	

	//=======================
        // End of timestep stuff
        //======================= 

        // Filter Solution
        solver->filterAndUpdateSolution();

        // Update primative and temperature
        solver->updateEndOfStepPrimAndTemp();

	// Update sponge BC averages and pressure
	solver->updateSpongeBCs();

        solver->checkSolution();

        if(solver->step%solver->outputStep == 0){
	    solver->dumpSolution();
        }

        solver->checkEnd();

    }
    return 0;
}

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
    int Ny = 100; double Ly = 4.0;//2*M_PI - 2*M_PI/((double)Ny);

    Solver *solver = new Solver(Nx, Ny, Lx, Ly);

    //Set some of the solver parameters
    solver->timeStep    = 1500;
    solver->timeEnd     = 10;
    solver->filterStep  = 1;
    solver->checkStep   = 1;
    solver->outputStep  = 100;
    solver->CFL		= 0.25;

    //Change fluid properties
    solver->idealGas->mu_ref = 0.000000001;

    //Set the BC's for the solver
    solver->bcX0 = Solver::SPONGE; //Solver::ADIABATIC_WALL;
    solver->bcX1 = Solver::SPONGE; // SPONGE;
    solver->bcY0 = Solver::PERIODIC;
    solver->bcY1 = Solver::PERIODIC;
    solver->setBCForDerivatives();

    //Allocate an initial condition
    for(int ip = 0; ip < Nx; ip++){
	for(int jp = 0; jp < Ny; jp++){
	    solver->U0[ip*Ny + jp] = 0.0;
	    solver->V0[ip*Ny + jp] = 0.0;
	    solver->p0[ip*Ny + jp] = 1/solver->idealGas->gamma + 
			0.035*exp(-(pow(solver->x[ip]-2.0,2.0))/0.01);

	    solver->rho0[ip*Ny + jp] = 1.0;

	}
    }

    //Apply the initial conditions
    solver->applyInitialCondition();

    //Provide non-standard attributes for the sponge layers
    solver->spongeStrength = 12.0;
    solver->spongeAvgT = 10.0;
    solver->spongeLengthX0 = 1.0;
    solver->spongeLengthY0 = 1.0;
    solver->spongeLengthX1 = 1.0;
    solver->spongeLengthY1 = 1.0;
    solver->spongeP = 1.0/solver->idealGas->gamma;


    //Need to run this with sponge BC's to initialize 
    solver->initSpongeStuff();

    //Dump up the initial condition
    solver->dumpSolution();


    double *test = new double[100];
    for(int ip = 0; ip < 100; ip++){
	test[ip] = 10.0;
    }


    while(!solver->endFlag){

        //At the beginning of every step
        solver->computeDtFromCFL_advanceTime();

	//=======================
        // RK Step 1
	//=======================

	solver->rkStep = 1;

	solver->handleBCs();

        solver->computeVelocityTemperatureGradients();

	solver->computeSpongeSource(solver->rho1, solver->rhoU1, solver->rhoV1, solver->rhoE1);	

        solver->computeContinuity(solver->rhoU1, solver->rhoV1);

        solver->computeXMomentum(solver->rhoU1, solver->rhoV1);

        solver->computeYMomentum(solver->rhoU1, solver->rhoV1);

        solver->computeEnergy(solver->rhoE1);

	solver->enforceBCs();

	solver->computeRhs();
    
        solver->updateSolutionRKStep1();


	//=======================
        // RK Step 2
	//=======================

	solver->rkStep = 2;

	solver->handleBCs();

        solver->computeVelocityTemperatureGradients();
	
	solver->computeSpongeSource(solver->rho1k, solver->rhoU1k, solver->rhoV1k, solver->rhoE1k);	

        solver->computeContinuity(solver->rhoU1k, solver->rhoV1k);

        solver->computeXMomentum(solver->rhoU1k, solver->rhoV1k);

        solver->computeYMomentum(solver->rhoU1k, solver->rhoV1k);

        solver->computeEnergy(solver->rhoE1k);

	solver->enforceBCs();
 
	solver->computeRhs();

        solver->updateSolutionRKStep2();



	//=======================
        // RK Step 3
	//=======================

	solver->rkStep = 3;

	solver->handleBCs();

        solver->computeVelocityTemperatureGradients();
	
	solver->computeSpongeSource(solver->rho1k, solver->rhoU1k, solver->rhoV1k, solver->rhoE1k);	

        solver->computeContinuity(solver->rhoU1k, solver->rhoV1k);

        solver->computeXMomentum(solver->rhoU1k, solver->rhoV1k);

        solver->computeYMomentum(solver->rhoU1k, solver->rhoV1k);

        solver->computeEnergy(solver->rhoE1k);

	solver->enforceBCs();

	solver->computeRhs();

        solver->updateSolutionRKStep3();


	//=======================
        // RK Step 4
	//=======================

	solver->rkStep = 4;

	solver->handleBCs();

        solver->computeVelocityTemperatureGradients();
	
	solver->computeSpongeSource(solver->rho1k, solver->rhoU1k, solver->rhoV1k, solver->rhoE1k);	

        solver->computeContinuity(solver->rhoU1k, solver->rhoV1k);

        solver->computeXMomentum(solver->rhoU1k, solver->rhoV1k);

        solver->computeYMomentum(solver->rhoU1k, solver->rhoV1k);

        solver->computeEnergy(solver->rhoE1k);

	solver->enforceBCs();

	solver->computeRhs();

        solver->updateSolutionRKStep4();
	

	//=======================
        // End of timestep stuff
        //======================= 

        // Filter Solution
        solver->filterAndUpdateSolution();

	// Update sponge BC averages and pressure
	solver->updateSpongeBCs();

        // Update primative and temperature
        solver->updateEndOfStepPrimAndTemp();

        solver->checkSolution();

        if(solver->step%solver->outputStep == 0){
	    solver->dumpSolution();
        }

        solver->checkEnd();

    }

    return 0;
}

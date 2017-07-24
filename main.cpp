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

    int Nx = 10; double Lx = 1; //2*M_PI - 2*M_PI/((double)Nx);
    int Ny = 10; double Ly = 1; //2*M_PI - 2*M_PI/((double)Ny);

    Solver *solver = new Solver(Nx, Ny, Lx, Ly);

    //Set some of the solver parameters
    solver->timeStep    = 100;
    solver->timeEnd     = 10.0;
    solver->filterStep  = 5;
    solver->checkStep   = 1;
    solver->outputStep  = 100;
    solver->CFL		= 0.3;

    //Change fluid properties
    solver->idealGas->mu_ref = 0.00001;

    //Set the BC's for the solver
    solver->bcX0 = Solver::SPONGE;
    solver->bcX1 = Solver::SPONGE;
    solver->bcY0 = Solver::PERIODIC;
    solver->bcY1 = Solver::PERIODIC;
    solver->setBCForDerivatives();

    //Allocate an initial condition
    for(int ip = 0; ip < Nx; ip++){
	for(int jp = 0; jp < Ny; jp++){
	    solver->p0[ip*Ny + jp]   = 1.0/1.4;
	    solver->U0[ip*Ny + jp]   = solver->x[ip];
	    solver->V0[ip*Ny + jp]   = solver->y[jp];
	    solver->rho0[ip*Ny + jp] = 1.0;
	}
    }

    solver->applyInitialCondition();

    solver->derivatives->CompactDYPeriodic(solver->U, solver->Uy);

    for(int ip = 0; ip < Nx; ip++){
	for(int jp = 0; jp < Ny; jp++){
	    cout << solver->Uy[ip*Ny + jp] << " ";
	}
	cout << endl;
    }
    cout << endl;

/*
    while(!solver->endFlag){

        //At the beginning of every step
        solver->computeDtFromCFL_advanceTime();

	//=======================
        // RK Step 1
	//=======================

        solver->computeVelocityTemperatureGradients();

        solver->computeContinuity(solver->rhoU1, solver->rhoV1);

        solver->computeXMomentum(solver->rhoU1, solver->rhoV1);

        solver->computeYMomentum(solver->rhoU1, solver->rhoV1);

        solver->computeEnergy(solver->rhoE1);
    
        solver->updateSolutionRKStep1();



	//=======================
        // RK Step 2
	//=======================

        solver->computeVelocityTemperatureGradients();
 
        solver->computeContinuity(solver->rhoU1k, solver->rhoV1k);

        solver->computeXMomentum(solver->rhoU1k, solver->rhoV1k);

        solver->computeYMomentum(solver->rhoU1k, solver->rhoV1k);

        solver->computeEnergy(solver->rhoE1k);
 
        solver->updateSolutionRKStep2();



	//=======================
        // RK Step 3
	//=======================

        solver->computeVelocityTemperatureGradients();

        solver->computeContinuity(solver->rhoU1k, solver->rhoV1k);

        solver->computeXMomentum(solver->rhoU1k, solver->rhoV1k);

        solver->computeYMomentum(solver->rhoU1k, solver->rhoV1k);

        solver->computeEnergy(solver->rhoE1k);

        solver->updateSolutionRKStep3();


	//=======================
        // RK Step 4
	//=======================

        solver->computeVelocityTemperatureGradients();

        solver->computeContinuity(solver->rhoU1k, solver->rhoV1k);

        solver->computeXMomentum(solver->rhoU1k, solver->rhoV1k);

        solver->computeYMomentum(solver->rhoU1k, solver->rhoV1k);

        solver->computeEnergy(solver->rhoE1k);

        solver->updateSolutionRKStep4();


	//=======================
        // End of timestep stuff
        //======================= 

        // Filter Solution
        solver->filterAndUpdateSolution();

        // Update primative and temperature
        solver->updateEndOfStepPrimAndTemp();

        solver->checkSolution();

        if(solver->step%solver->outputStep == 0){
	    solver->dumpSolution();
        }

        solver->checkEnd();

    }
*/

    return 0;
}

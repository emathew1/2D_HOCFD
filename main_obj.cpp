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
    cout << " Testing 2-D High Order Solver w/ OOP " << endl;
    cout << "--------------------------------------" << endl;

    int Nx = 10; double Lx = 2*M_PI - 2*M_PI/((double)Nx);
    int Ny = 10; double Ly = 2*M_PI - 2*M_PI/((double)Ny);

    Solver *solver = new Solver(Nx, Ny, Lx, Ly);

    //Set some of the solver parameters
    solver->timeStep    = 100000;
    solver->timeEnd     = 10.0;
    solver->filterStep  = 5;
    solver->checkStep   = 1;
    solver->outputStep  = 1;
    solver->CFL		= 0.3;

    //Change fluid properties
    solver->idealGas->mu_ref = 0.00001;

    //Allocate an initial condition
    for(int ip = 0; ip < Nx; ip++){
	for(int jp = 0; jp < Ny; jp++){
	    solver->p0[ip*Ny + jp]   = 1.0/1.4;
	    solver->U0[ip*Ny + jp]   = cos(solver->x[ip]);
	    solver->V0[ip*Ny + jp]   = sin(solver->x[ip]);;
	    solver->rho0[ip*Ny + jp] = 1.0;
	}
    }


    solver->applyInitialCondition();
    solver->computeDtFromCFL_advanceTime();

    solver->computeVelocityTemperatureGradients();

    std::cout.precision(18);

    for(int ip = 0; ip < Nx; ip++){
	for(int jp = 0; jp < Ny; jp++){
	    cout << solver->Vx[ip*Ny + jp] << " ";
	}
	cout << endl;
    }

/*
    //Time and step stuff
    double dt   = 0;
    double time = 0;
    int    step = 0;
    int    checkStep = 1;

    //Mesh Parameters
    const int Nx   = 3072;
    const int Ny   = 3072;
    const double Lx = 1.0;
    const double Ly = 1.0;
    cout << " -> Creating " << Nx << "x" << Ny << " mesh " << endl;
    cout << " -> Dimensions: " << Lx << "x" << Ly << endl;

    //Fluid Parameters
    const double gamma   = 1.4;
    const double Pr      = 0.7;
    const double p_ref   = 1.0/gamma;
    const double rho_ref = 1.0;
    const double T_ref   = 1.0;
    const double mu_ref  = 0.00000001;
    const double R_gas   = p_ref/rho_ref/T_ref;
    const double cp      = R_gas*gamma/(gamma-1.0);
    cout << " -> Using gamma = " << gamma << endl; 
    cout << " ->       Pr = " << Pr << endl;

    //Solver Parameters
    const int timeStep = 100000;
    const double timeEnd = 10.0;
    const int filterStep = 5;
    const int outputStep = 1;
    const double CFL = 0.3;
    cout << " -> Running " << timeStep << " time steps" << endl; 
    cout << " -> Filtering every " << filterStep << " time steps" << endl; 
    cout << " -> constant CFL of " << CFL << endl; 

    //Start with 6th order tridiagonal solver
    const double alpha = 1.0/3.0;
    const double beta  = 0.0;
    const double a     = 14.0/9.0; //(1.0/6.0) *(     alpha + 9.0);
    const double b     = 1.0/9.0;  //(1.0/15.0)*(32.0*alpha - 9.0);
    const double c     = 0.0; // (1.0/10.0)*(-3.0*alpha + 1.0); 
    cout << " -> Using 6th order tridiagonal solver " << endl;

    //8th-order tridiagonal filter
    const double alphaF = 0.49;
    const double betaF  = 0.0;
    const double a0     = (93.0 + 70.0*alphaF)/128.0;
    const double a1     = ( 7.0 + 18.0*alphaF)/16.0;
    const double a2     = (-7.0 + 14.0*alphaF)/32.0;
    const double a3     = ( 1.0 -  2.0*alphaF)/16.0;
    const double a4     = (-1.0 +  2.0*alphaF)/128.0;
    const double a5     = 0.0;
    cout << " -> Using 8th order tridiagonal filter " << endl;

    //Allocate the data
 
    //Just the linear x and y grids
    double *x  = new double[Nx];
    double *y  = new double[Ny];
   
    for(int ip = 0; ip < Nx; ip++){
	x[ip] = (((double)ip)/((double)Nx - 1.0))*Lx;
    }
    

    for(int ip = 0; ip < Ny; ip++){
	y[ip] = (((double)ip)/((double)Ny - 1.0))*Ly;
    }

    double dx = x[1]-x[0]; 
    double dy = y[1]-y[0]; 

    //Keep track of number of arrays we've allocated
    int nArray = 0;
    cout << " -> Allocating the data..." << endl; 

    //Initial Condition Arrays
    double *rho0 = new double[Nx*Ny]; nArray++;
    double *U0   = new double[Nx*Ny]; nArray++;
    double *V0   = new double[Nx*Ny]; nArray++;
    double *p0   = new double[Nx*Ny]; nArray++;

    //Solver Data Arrays 
    //Derived Data
    double *mu   = new double[Nx*Ny]; nArray++;
    double *U    = new double[Nx*Ny]; nArray++;
    double *V    = new double[Nx*Ny]; nArray++;
    double *SOS  = new double[Nx*Ny]; nArray++;
    double *p    = new double[Nx*Ny]; nArray++;
    double *T    = new double[Nx*Ny]; nArray++;
    double *Ux   = new double[Nx*Ny]; nArray++;
    double *Uy   = new double[Nx*Ny]; nArray++;
    double *Vx   = new double[Nx*Ny]; nArray++;
    double *Vy   = new double[Nx*Ny]; nArray++;
    double *Tx   = new double[Nx*Ny]; nArray++;
    double *Ty   = new double[Nx*Ny]; nArray++; 


    //SolverStepData
    double *rho1   = new double[Nx*Ny]; nArray++;
    double *rhok   = new double[Nx*Ny]; nArray++;
    double *rho1k  = new double[Nx*Ny]; nArray++;
    double *rho2   = new double[Nx*Ny]; nArray++;

    double *rhoU1  = new double[Nx*Ny]; nArray++;
    double *rhoUk  = new double[Nx*Ny]; nArray++;
    double *rhoU1k = new double[Nx*Ny]; nArray++;
    double *rhoU2  = new double[Nx*Ny]; nArray++;
 
    double *rhoV1  = new double[Nx*Ny]; nArray++;
    double *rhoVk  = new double[Nx*Ny]; nArray++;
    double *rhoV1k = new double[Nx*Ny]; nArray++;
    double *rhoV2  = new double[Nx*Ny]; nArray++;

    double *rhoE1  = new double[Nx*Ny]; nArray++;
    double *rhoEk  = new double[Nx*Ny]; nArray++;
    double *rhoE1k = new double[Nx*Ny]; nArray++;
    double *rhoE2  = new double[Nx*Ny]; nArray++;

    double *rhsDxIn   = new double[Nx*Ny]; nArray++;
    double *rhsDxOut  = new double[Nx*Ny]; nArray++;
    double *rhsDyIn   = new double[Nx*Ny]; nArray++;
    double *rhsDyOut  = new double[Nx*Ny]; nArray++;
 
    //As of now we've allocation nArrays Nx*Ny arrays, how much data is that?
    long long unsigned int memAllocation = 0;
    memAllocation = nArray*Nx*Ny; 
    cout << " -> Allocated " << nArray << " arrays for the solver ";
    cout << ", requiring " << (double)memAllocation/1024.0/1024.0/1024.0 << " Gb of memory" << endl;

    cout << " -> Setting RUP initial condition..." << endl;
    //Set initial condition
    for(int ip = 0; ip < Nx; ip++){
	for(int jp = 0; jp < Ny; jp++){


            p0[ip*Ny + jp]   = 2.5/gamma;

            if(y[jp] > 0.25*Ly && y[jp] < 0.75*Ly){
	        U0[ip*Ny + jp] = 0.5;
                rho0[ip*Ny + jp] = 2.0;
            }else{
	        U0[ip*Ny + jp] = -0.5;
                rho0[ip*Ny + jp] = 1.0;
            }

	    double sigma = 0.05/sqrt(2.0);
	    V0[ip*Ny + jp] = 0.1*sin(4*M_PI*x[ip])*(exp(-(pow(y[jp]-0.25,2)/(2*sigma*sigma))) + 
				exp(-(pow(y[jp]-0.75,2)/(2*sigma*sigma))));
	}
    }

    //Set up Tridiagonal Solvers...

    cout << " -> Setting up compact diff and filter matrices..." << endl;
    //Set initial condition
    //LHS Differencing in X
    double *diagx = new double[Nx];
    double *offx  = new double[Nx];

    //LHS Differencing in Y
    double *diagy = new double[Ny];
    double *offy  = new double[Ny];

    //LHS Filtering in X
    double *diagFx = new double[Nx];
    double *offFx  = new double[Nx];

    //LHS Filtering in Y
    double *diagFy = new double[Ny];
    double *offFy  = new double[Ny];

    for (int i = 0; i < Nx; i++)
    {
	diagx[i]  = 1.0;
	diagFx[i] = 1.0;
	offx[i]   = alpha;
	offFx[i]  = alphaF;
    }

    for (int i = 0; i < Ny; i++)
    {
	diagy[i]  = 1.0;
	diagFy[i] = 1.0;
	offy[i]   = alpha;
	offFy[i]  = alphaF;
    }
    cout << " -> Starting time stepping! " << endl << endl;

*/

/*
    auto t1 = high_resolution_clock::now();
    auto t2 = high_resolution_clock::now();
    int filterCount = 0;
    for(int ip = 0; ip < timeStep; ip++){

	//Apply the initial condition
	if(ip == 0){
	    for(int jp = 0; jp < Nx*Ny; jp++){
		U[jp]     = U0[jp];
		V[jp]     = V0[jp];
		rho1[jp]  = rho0[jp];
		p[jp]     = p0[jp];
		rhoU1[jp] = rho0[jp]*U0[jp];
		rhoV1[jp] = rho0[jp]*V0[jp];
		rhoE1[jp] = p0[jp]/(gamma-1) + 0.5*rho0[jp]*(U0[jp]*U0[jp]+V0[jp]*V0[jp]);
		T[jp]     = p0[jp]/(rho0[jp]*R_gas);
		mu[jp]    = mu_ref*pow(T[jp]/T_ref, 0.76);
		SOS[jp]   = sqrt(gamma*p0[jp]/rho0[jp]);

	    }
   	    delete[] rho0;
	    delete[] p0;
	    delete[] U0;
	    delete[] V0;
	}	

	//Compute the dt from CFL condition
	double UChar = -100000; double tempU;
	double VChar = -100000; double tempV;
	for(int jp = 0; jp < Nx*Ny; jp++){
	    tempU = abs(U[jp]+SOS[jp]);
	    if(tempU >= UChar){
		UChar = tempU;
	    }
	    tempV = abs(V[jp]+SOS[jp]);
	    if(tempV >= VChar){
		VChar = tempV;
	    }
	}
	dt = fabs(min(CFL*dx/UChar, CFL*dy/VChar));
	time += dt;


	////////////RK STEP 1///////////////

	//Compute the Velocity Gradients
        CompactDXPeriodic(U, diagx, offx, alpha, a, b, c, dx, Nx, Ny, Ux);
        CompactDXPeriodic(V, diagx, offx, alpha, a, b, c, dx, Nx, Ny, Vx);
        CompactDYPeriodic(U, diagy, offy, alpha, a, b, c, dy, Nx, Ny, Uy);
        CompactDYPeriodic(V, diagy, offy, alpha, a, b, c, dy, Nx, Ny, Vy);
  
        //Compute the Temperature Gradients	
        CompactDXPeriodic(T, diagx, offx, alpha, a, b, c, dx, Nx, Ny, Tx);
        CompactDYPeriodic(T, diagy, offy, alpha, a, b, c, dy, Nx, Ny, Ty);

	//Continuity
        CompactDXPeriodic(rhoU1, diagx, offx, alpha, a, b, c, dx, Nx, Ny, rhsDxOut);
        CompactDYPeriodic(rhoV1, diagy, offy, alpha, a, b, c, dy, Nx, Ny, rhsDyOut);
	for(int jp = 0; jp < Nx*Ny; jp++){
	    rhok[jp] = dt*(-rhsDxOut[jp]-rhsDyOut[jp]);
	}	


	//X-Momentum
	for(int jp = 0; jp < Nx*Ny; jp++){
	    rhsDxIn[jp] = -(rhoU1[jp]*U[jp] + p[jp] - 2.0*mu[jp]*Ux[jp] + (2.0/3.0)*mu[jp]*(Ux[jp] + Vy[jp])); 
	    rhsDyIn[jp] = -(rhoV1[jp]*U[jp] - mu[jp]*(Uy[jp] + Vx[jp])); 
	}
        CompactDXPeriodic(rhsDxIn, diagx, offx, alpha, a, b, c, dx, Nx, Ny, rhsDxOut);
        CompactDYPeriodic(rhsDyIn, diagy, offy, alpha, a, b, c, dy, Nx, Ny, rhsDyOut);
	for(int jp = 0; jp < Nx*Ny; jp++){
	    rhoUk[jp] = dt*(rhsDxOut[jp]+rhsDyOut[jp]);
	}	


	//Y-Momentum
	for(int jp = 0; jp < Nx*Ny; jp++){
	    rhsDxIn[jp] = -(rhoU1[jp]*V[jp] - mu[jp]*(Uy[jp] + Vx[jp])); 
	    rhsDyIn[jp] = -(rhoV1[jp]*V[jp] + p[jp] - 2.0*mu[jp]*Vy[jp] + (2.0/3.0)*mu[jp]*(Ux[jp] + Vy[jp])); 
	}
        CompactDXPeriodic(rhsDxIn, diagx, offx, alpha, a, b, c, dx, Nx, Ny, rhsDxOut);
        CompactDYPeriodic(rhsDyIn, diagy, offy, alpha, a, b, c, dy, Nx, Ny, rhsDyOut);
	for(int jp = 0; jp < Nx*Ny; jp++){
	    rhoVk[jp] = dt*(rhsDxOut[jp]+rhsDyOut[jp]);
	}	

	//Energy
	for(int jp = 0; jp < Nx*Ny; jp++){
	    rhsDxIn[jp] = -(rhoE1[jp]*U[jp] + U[jp]*p[jp] - (mu[jp]/Pr/(gamma-1.0))*Tx[jp] +
				-U[jp]*(2.0*mu[jp]*Ux[jp] - (2.0/3.0)*mu[jp]*(Ux[jp]+Vy[jp])) +
				-V[jp]*(mu[jp]*(Vx[jp] + Uy[jp]))); 
	    rhsDyIn[jp] = -(rhoE1[jp]*V[jp] + V[jp]*p[jp] - (mu[jp]/Pr/(gamma-1.0))*Ty[jp] +
				-V[jp]*(2.0*mu[jp]*Vy[jp] - (2.0/3.0)*mu[jp]*(Ux[jp]+Vy[jp])) +
				-U[jp]*(mu[jp]*(Vx[jp] + Uy[jp]))); 
	}
        CompactDXPeriodic(rhsDxIn, diagx, offx, alpha, a, b, c, dx, Nx, Ny, rhsDxOut);
        CompactDYPeriodic(rhsDyIn, diagy, offy, alpha, a, b, c, dy, Nx, Ny, rhsDyOut);
	for(int jp = 0; jp < Nx*Ny; jp++){
	    rhoEk[jp] = dt*(rhsDxOut[jp]+rhsDyOut[jp]);
	}	

	//Update solution at the end of step 1
	for(int jp = 0; jp < Nx*Ny; jp++){
	    rho2[jp]  = rho1[jp]  + rhok[jp]/6.0;
	    rhoU2[jp] = rhoU1[jp] + rhoUk[jp]/6.0;
	    rhoV2[jp] = rhoV1[jp] + rhoVk[jp]/6.0;
	    rhoE2[jp] = rhoE1[jp] + rhoEk[jp]/6.0;
	}

	//Solve for derived data

	//intermediate solution at the end of step 1
	for(int jp = 0; jp < Nx*Ny; jp++){
	    rho1k[jp]  = rho1[jp] + rhok[jp]/2.0;
	    rhoU1k[jp] = rhoU1[jp] + rhoUk[jp]/2.0;
	    rhoV1k[jp] = rhoV1[jp] + rhoVk[jp]/2.0;
	    rhoE1k[jp] = rhoE1[jp] + rhoEk[jp]/2.0;
	}
	
	//Solve for U and V velocities
	solveU(rho1k, rhoU1k, Nx, Ny, U);
	solveU(rho1k, rhoV1k, Nx, Ny, V);

	//Solve for P
	solvep(rho1k, rhoE1k, U, V, gamma, Nx, Ny, p);

	//Solve for T
	solveT(rho1k, p, R_gas, Nx, Ny, T);

	//Solve for mu
	solvemu(T, T_ref, mu_ref, Nx, Ny, mu);


	///////////RK STEP 2////////////////////
        //Compute the Velocity Gradients
        CompactDXPeriodic(U, diagx, offx, alpha, a, b, c, dx, Nx, Ny, Ux);
        CompactDXPeriodic(V, diagx, offx, alpha, a, b, c, dx, Nx, Ny, Vx);
        CompactDYPeriodic(U, diagy, offy, alpha, a, b, c, dy, Nx, Ny, Uy);
        CompactDYPeriodic(V, diagy, offy, alpha, a, b, c, dy, Nx, Ny, Vy);

        //Compute the Temperature Gradients     
        CompactDXPeriodic(T, diagx, offx, alpha, a, b, c, dx, Nx, Ny, Tx);
        CompactDYPeriodic(T, diagy, offy, alpha, a, b, c, dy, Nx, Ny, Ty);

        //Continuity
        CompactDXPeriodic(rhoU1k, diagx, offx, alpha, a, b, c, dx, Nx, Ny, rhsDxOut);
        CompactDYPeriodic(rhoV1k, diagy, offy, alpha, a, b, c, dy, Nx, Ny, rhsDyOut);
        for(int jp = 0; jp < Nx*Ny; jp++){
            rhok[jp] = dt*(-rhsDxOut[jp]-rhsDyOut[jp]);
        }

        //X-Momentum
        for(int jp = 0; jp < Nx*Ny; jp++){
            rhsDxIn[jp] = -(rhoU1k[jp]*U[jp] + p[jp] - 2.0*mu[jp]*Ux[jp] + (2.0/3.0)*mu[jp]*(Ux[jp] + Vy[jp]));
            rhsDyIn[jp] = -(rhoV1k[jp]*U[jp] - mu[jp]*(Uy[jp] + Vx[jp]));
        }
        CompactDXPeriodic(rhsDxIn, diagx, offx, alpha, a, b, c, dx, Nx, Ny, rhsDxOut);
        CompactDYPeriodic(rhsDyIn, diagy, offy, alpha, a, b, c, dy, Nx, Ny, rhsDyOut);
        for(int jp = 0; jp < Nx*Ny; jp++){
            rhoUk[jp] = dt*(rhsDxOut[jp]+rhsDyOut[jp]);
        }      


        //Y-Momentum
        for(int jp = 0; jp < Nx*Ny; jp++){
            rhsDxIn[jp] = -(rhoU1k[jp]*V[jp] - mu[jp]*(Uy[jp] + Vx[jp]));
            rhsDyIn[jp] = -(rhoV1k[jp]*V[jp] + p[jp] - 2.0*mu[jp]*Vy[jp] + (2.0/3.0)*mu[jp]*(Ux[jp] + Vy[jp]));
        }
        CompactDXPeriodic(rhsDxIn, diagx, offx, alpha, a, b, c, dx, Nx, Ny, rhsDxOut);
        CompactDYPeriodic(rhsDyIn, diagy, offy, alpha, a, b, c, dy, Nx, Ny, rhsDyOut);
        for(int jp = 0; jp < Nx*Ny; jp++){
            rhoVk[jp] = dt*(rhsDxOut[jp]+rhsDyOut[jp]);
        }

        //Energy
        for(int jp = 0; jp < Nx*Ny; jp++){
            rhsDxIn[jp] = -(rhoE1k[jp]*U[jp] + U[jp]*p[jp] - (mu[jp]/Pr/(gamma-1.0))*Tx[jp] +
                                -U[jp]*(2.0*mu[jp]*Ux[jp] - (2.0/3.0)*mu[jp]*(Ux[jp]+Vy[jp])) +
                                -V[jp]*(mu[jp]*(Vx[jp] + Uy[jp])));
            rhsDyIn[jp] = -(rhoE1k[jp]*V[jp] + V[jp]*p[jp] - (mu[jp]/Pr/(gamma-1.0))*Ty[jp] +
                                -V[jp]*(2.0*mu[jp]*Vy[jp] - (2.0/3.0)*mu[jp]*(Ux[jp]+Vy[jp])) +
                                -U[jp]*(mu[jp]*(Vx[jp] + Uy[jp])));
        }
        CompactDXPeriodic(rhsDxIn, diagx, offx, alpha, a, b, c, dx, Nx, Ny, rhsDxOut);
        CompactDYPeriodic(rhsDyIn, diagy, offy, alpha, a, b, c, dy, Nx, Ny, rhsDyOut);
        for(int jp = 0; jp < Nx*Ny; jp++){
            rhoEk[jp] = dt*(rhsDxOut[jp]+rhsDyOut[jp]);
        }

	//Update final solution at the end of step 2
	for(int jp = 0; jp < Nx*Ny; jp++){
	    rho2[jp]  += rhok[jp]/3.0;
	    rhoU2[jp] += rhoUk[jp]/3.0;
	    rhoV2[jp] += rhoVk[jp]/3.0;
	    rhoE2[jp] += rhoEk[jp]/3.0;
	}

	//Update intermediate solution at the end of step 2
	for(int jp = 0; jp < Nx*Ny; jp++){
	    rho1k[jp]  = rho1[jp] + rhok[jp]/2.0;
	    rhoU1k[jp] = rhoU1[jp] + rhoUk[jp]/2.0;
	    rhoV1k[jp] = rhoV1[jp] + rhoVk[jp]/2.0;
	    rhoE1k[jp] = rhoE1[jp] + rhoEk[jp]/2.0;
	}
	
	//Solve for U and V velocities
	solveU(rho1k, rhoU1k, Nx, Ny, U);
	solveU(rho1k, rhoV1k, Nx, Ny, V);

	//Solve for P
	solvep(rho1k, rhoE1k, U, V, gamma, Nx, Ny, p);

	//Solve for T
	solveT(rho1k, p, R_gas, Nx, Ny, T);

	//Solve for mu
	solvemu(T, T_ref, mu_ref, Nx, Ny, mu);


	///////////RK STEP 3////////////////////
        //Compute the Velocity Gradients
        CompactDXPeriodic(U, diagx, offx, alpha, a, b, c, dx, Nx, Ny, Ux);
        CompactDXPeriodic(V, diagx, offx, alpha, a, b, c, dx, Nx, Ny, Vx);
        CompactDYPeriodic(U, diagy, offy, alpha, a, b, c, dy, Nx, Ny, Uy);
        CompactDYPeriodic(V, diagy, offy, alpha, a, b, c, dy, Nx, Ny, Vy);

        //Compute the Temperature Gradients     
        CompactDXPeriodic(T, diagx, offx, alpha, a, b, c, dx, Nx, Ny, Tx);
        CompactDYPeriodic(T, diagy, offy, alpha, a, b, c, dy, Nx, Ny, Ty);

        //Continuity
        CompactDXPeriodic(rhoU1k, diagx, offx, alpha, a, b, c, dx, Nx, Ny, rhsDxOut);
        CompactDYPeriodic(rhoV1k, diagy, offy, alpha, a, b, c, dy, Nx, Ny, rhsDyOut);
        for(int jp = 0; jp < Nx*Ny; jp++){
            rhok[jp] = dt*(-rhsDxOut[jp]-rhsDyOut[jp]);
        }

        //X-Momentum
        for(int jp = 0; jp < Nx*Ny; jp++){
            rhsDxIn[jp] = -(rhoU1k[jp]*U[jp] + p[jp] - 2.0*mu[jp]*Ux[jp] + (2.0/3.0)*mu[jp]*(Ux[jp] + Vy[jp]));
            rhsDyIn[jp] = -(rhoV1k[jp]*U[jp] - mu[jp]*(Uy[jp] + Vx[jp]));
        }
        CompactDXPeriodic(rhsDxIn, diagx, offx, alpha, a, b, c, dx, Nx, Ny, rhsDxOut);
        CompactDYPeriodic(rhsDyIn, diagy, offy, alpha, a, b, c, dy, Nx, Ny, rhsDyOut);
        for(int jp = 0; jp < Nx*Ny; jp++){
            rhoUk[jp] = dt*(rhsDxOut[jp]+rhsDyOut[jp]);
        }      


        //Y-Momentum
        for(int jp = 0; jp < Nx*Ny; jp++){
            rhsDxIn[jp] = -(rhoU1k[jp]*V[jp] - mu[jp]*(Uy[jp] + Vx[jp]));
            rhsDyIn[jp] = -(rhoV1k[jp]*V[jp] + p[jp] - 2.0*mu[jp]*Vy[jp] + (2.0/3.0)*mu[jp]*(Ux[jp] + Vy[jp]));
        }
        CompactDXPeriodic(rhsDxIn, diagx, offx, alpha, a, b, c, dx, Nx, Ny, rhsDxOut);
        CompactDYPeriodic(rhsDyIn, diagy, offy, alpha, a, b, c, dy, Nx, Ny, rhsDyOut);
        for(int jp = 0; jp < Nx*Ny; jp++){
            rhoVk[jp] = dt*(rhsDxOut[jp]+rhsDyOut[jp]);
        }

        //Energy
        for(int jp = 0; jp < Nx*Ny; jp++){
            rhsDxIn[jp] = -(rhoE1k[jp]*U[jp] + U[jp]*p[jp] - (mu[jp]/Pr/(gamma-1.0))*Tx[jp] +
                                -U[jp]*(2.0*mu[jp]*Ux[jp] - (2.0/3.0)*mu[jp]*(Ux[jp]+Vy[jp])) +
                                -V[jp]*(mu[jp]*(Vx[jp] + Uy[jp])));
            rhsDyIn[jp] = -(rhoE1k[jp]*V[jp] + V[jp]*p[jp] - (mu[jp]/Pr/(gamma-1.0))*Ty[jp] +
                                -V[jp]*(2.0*mu[jp]*Vy[jp] - (2.0/3.0)*mu[jp]*(Ux[jp]+Vy[jp])) +
                                -U[jp]*(mu[jp]*(Vx[jp] + Uy[jp])));
        }
        CompactDXPeriodic(rhsDxIn, diagx, offx, alpha, a, b, c, dx, Nx, Ny, rhsDxOut);
        CompactDYPeriodic(rhsDyIn, diagy, offy, alpha, a, b, c, dy, Nx, Ny, rhsDyOut);
        for(int jp = 0; jp < Nx*Ny; jp++){
            rhoEk[jp] = dt*(rhsDxOut[jp]+rhsDyOut[jp]);
        }

	//Update final solution at the end of step 3
	for(int jp = 0; jp < Nx*Ny; jp++){
	    rho2[jp]  += rhok[jp]/3.0;
	    rhoU2[jp] += rhoUk[jp]/3.0;
	    rhoV2[jp] += rhoVk[jp]/3.0;
	    rhoE2[jp] += rhoEk[jp]/3.0;
	}

	//Update intermediate solution at the end of step 3
	for(int jp = 0; jp < Nx*Ny; jp++){
	    rho1k[jp]  = rho1[jp] + rhok[jp];
	    rhoU1k[jp] = rhoU1[jp] + rhoUk[jp];
	    rhoV1k[jp] = rhoV1[jp] + rhoVk[jp];
	    rhoE1k[jp] = rhoE1[jp] + rhoEk[jp];
	}
	
	//Solve for U and V velocities
	solveU(rho1k, rhoU1k, Nx, Ny, U);
	solveU(rho1k, rhoV1k, Nx, Ny, V);

	//Solve for P
	solvep(rho1k, rhoE1k, U, V, gamma, Nx, Ny, p);

	//Solve for T
	solveT(rho1k, p, R_gas, Nx, Ny, T);

	//Solve for mu
	solvemu(T, T_ref, mu_ref, Nx, Ny, mu);


	///////////RK STEP 2////////////////////
        //Compute the Velocity Gradients
        CompactDXPeriodic(U, diagx, offx, alpha, a, b, c, dx, Nx, Ny, Ux);
        CompactDXPeriodic(V, diagx, offx, alpha, a, b, c, dx, Nx, Ny, Vx);
        CompactDYPeriodic(U, diagy, offy, alpha, a, b, c, dy, Nx, Ny, Uy);
        CompactDYPeriodic(V, diagy, offy, alpha, a, b, c, dy, Nx, Ny, Vy);

        //Compute the Temperature Gradients     
        CompactDXPeriodic(T, diagx, offx, alpha, a, b, c, dx, Nx, Ny, Tx);
        CompactDYPeriodic(T, diagy, offy, alpha, a, b, c, dy, Nx, Ny, Ty);

        //Continuity
        CompactDXPeriodic(rhoU1k, diagx, offx, alpha, a, b, c, dx, Nx, Ny, rhsDxOut);
        CompactDYPeriodic(rhoV1k, diagy, offy, alpha, a, b, c, dy, Nx, Ny, rhsDyOut);
        for(int jp = 0; jp < Nx*Ny; jp++){
            rhok[jp] = dt*(-rhsDxOut[jp]-rhsDyOut[jp]);
        }

        //X-Momentum
        for(int jp = 0; jp < Nx*Ny; jp++){
            rhsDxIn[jp] = -(rhoU1k[jp]*U[jp] + p[jp] - 2.0*mu[jp]*Ux[jp] + (2.0/3.0)*mu[jp]*(Ux[jp] + Vy[jp]));
            rhsDyIn[jp] = -(rhoV1k[jp]*U[jp] - mu[jp]*(Uy[jp] + Vx[jp]));
        }
        CompactDXPeriodic(rhsDxIn, diagx, offx, alpha, a, b, c, dx, Nx, Ny, rhsDxOut);
        CompactDYPeriodic(rhsDyIn, diagy, offy, alpha, a, b, c, dy, Nx, Ny, rhsDyOut);
        for(int jp = 0; jp < Nx*Ny; jp++){
            rhoUk[jp] = dt*(rhsDxOut[jp]+rhsDyOut[jp]);
        }      


        //Y-Momentum
        for(int jp = 0; jp < Nx*Ny; jp++){
            rhsDxIn[jp] = -(rhoU1k[jp]*V[jp] - mu[jp]*(Uy[jp] + Vx[jp]));
            rhsDyIn[jp] = -(rhoV1k[jp]*V[jp] + p[jp] - 2.0*mu[jp]*Vy[jp] + (2.0/3.0)*mu[jp]*(Ux[jp] + Vy[jp]));
        }
        CompactDXPeriodic(rhsDxIn, diagx, offx, alpha, a, b, c, dx, Nx, Ny, rhsDxOut);
        CompactDYPeriodic(rhsDyIn, diagy, offy, alpha, a, b, c, dy, Nx, Ny, rhsDyOut);
        for(int jp = 0; jp < Nx*Ny; jp++){
            rhoVk[jp] = dt*(rhsDxOut[jp]+rhsDyOut[jp]);
        }

        //Energy
        for(int jp = 0; jp < Nx*Ny; jp++){
            rhsDxIn[jp] = -(rhoE1k[jp]*U[jp] + U[jp]*p[jp] - (mu[jp]/Pr/(gamma-1.0))*Tx[jp] +
                                -U[jp]*(2.0*mu[jp]*Ux[jp] - (2.0/3.0)*mu[jp]*(Ux[jp]+Vy[jp])) +
                                -V[jp]*(mu[jp]*(Vx[jp] + Uy[jp])));
            rhsDyIn[jp] = -(rhoE1k[jp]*V[jp] + V[jp]*p[jp] - (mu[jp]/Pr/(gamma-1.0))*Ty[jp] +
                                -V[jp]*(2.0*mu[jp]*Vy[jp] - (2.0/3.0)*mu[jp]*(Ux[jp]+Vy[jp])) +
                                -U[jp]*(mu[jp]*(Vx[jp] + Uy[jp])));
        }
        CompactDXPeriodic(rhsDxIn, diagx, offx, alpha, a, b, c, dx, Nx, Ny, rhsDxOut);
        CompactDYPeriodic(rhsDyIn, diagy, offy, alpha, a, b, c, dy, Nx, Ny, rhsDyOut);
        for(int jp = 0; jp < Nx*Ny; jp++){
            rhoEk[jp] = dt*(rhsDxOut[jp]+rhsDyOut[jp]);
        }

	//Get final solution at the end of step 4
	for(int jp = 0; jp < Nx*Ny; jp++){
	    rho2[jp]  += rhok[jp]/6.0;
	    rhoU2[jp] += rhoUk[jp]/6.0;
	    rhoV2[jp] += rhoVk[jp]/6.0;
	    rhoE2[jp] += rhoEk[jp]/6.0;
	}

	

	if(ip%filterStep == 0){
	    //Run filtering of conservative data and dump into new solution container
	    //Alternate first filter direction 
	    if(filterCount%2 == 0){

		double *rhoTemp = new double[Nx*Ny];
		double *rhoUTemp = new double[Nx*Ny];
		double *rhoVTemp = new double[Nx*Ny];
		double *rhoETemp = new double[Nx*Ny];
		FilterPeriodicX(rho2, diagFx, offFx, alphaF, a0, a1, a2, a3, a4, a5, Nx, Ny, rhoTemp);
		FilterPeriodicX(rhoU2, diagFx, offFx, alphaF, a0, a1, a2, a3, a4, a5, Nx, Ny, rhoUTemp);
		FilterPeriodicX(rhoV2, diagFx, offFx, alphaF, a0, a1, a2, a3, a4, a5, Nx, Ny, rhoVTemp);
		FilterPeriodicX(rhoE2, diagFx, offFx, alphaF, a0, a1, a2, a3, a4, a5, Nx, Ny, rhoETemp);

		FilterPeriodicY(rhoTemp, diagFy, offFy, alphaF, a0, a1, a2, a3, a4, a5, Nx, Ny, rho1);
		FilterPeriodicY(rhoUTemp, diagFy, offFy, alphaF, a0, a1, a2, a3, a4, a5, Nx, Ny, rhoU1);
		FilterPeriodicY(rhoVTemp, diagFy, offFy, alphaF, a0, a1, a2, a3, a4, a5, Nx, Ny, rhoV1);
		FilterPeriodicY(rhoETemp, diagFy, offFy, alphaF, a0, a1, a2, a3, a4, a5, Nx, Ny, rhoE1);
		delete[] rhoTemp;
		delete[] rhoUTemp;
		delete[] rhoVTemp;
		delete[] rhoETemp;
	    }else{

		double *rhoTemp = new double[Nx*Ny];
		double *rhoUTemp = new double[Nx*Ny];
		double *rhoVTemp = new double[Nx*Ny];
		double *rhoETemp = new double[Nx*Ny];
		FilterPeriodicY(rho2, diagFy, offFy, alphaF, a0, a1, a2, a3, a4, a5, Nx, Ny, rhoTemp);
		FilterPeriodicY(rhoU2, diagFy, offFy, alphaF, a0, a1, a2, a3, a4, a5, Nx, Ny, rhoUTemp);
		FilterPeriodicY(rhoV2, diagFy, offFy, alphaF, a0, a1, a2, a3, a4, a5, Nx, Ny, rhoVTemp);
		FilterPeriodicY(rhoE2, diagFy, offFy, alphaF, a0, a1, a2, a3, a4, a5, Nx, Ny, rhoETemp);

		FilterPeriodicX(rhoTemp, diagFx, offFx, alphaF, a0, a1, a2, a3, a4, a5, Nx, Ny, rho1);
		FilterPeriodicX(rhoUTemp, diagFx, offFx, alphaF, a0, a1, a2, a3, a4, a5, Nx, Ny, rhoU1);
		FilterPeriodicX(rhoVTemp, diagFx, offFx, alphaF, a0, a1, a2, a3, a4, a5, Nx, Ny, rhoV1);
		FilterPeriodicX(rhoETemp, diagFx, offFx, alphaF, a0, a1, a2, a3, a4, a5, Nx, Ny, rhoE1);
		delete[] rhoTemp;
		delete[] rhoUTemp;
		delete[] rhoVTemp;
		delete[] rhoETemp;
	    }
	    filterCount++;
	}else{

	    //Just Update the solutions
	    memcpy(rho1, rho2, Nx*Ny*sizeof(double));	
	    memcpy(rhoU1, rhoU2, Nx*Ny*sizeof(double));	
	    memcpy(rhoV1, rhoV2, Nx*Ny*sizeof(double));	
	    memcpy(rhoE1, rhoE2, Nx*Ny*sizeof(double));	
	}


	//Solve for derived data
	
	//Solve for U and V velocities
	solveU(rho1, rhoU1, Nx, Ny, U);
	solveU(rho1, rhoV1, Nx, Ny, V);

	//Solve for P
	solvep(rho1, rhoE1, U, V, gamma, Nx, Ny, p);

	//Solve for T
	solveT(rho1, p, R_gas, Nx, Ny, T);

	//Solve for mu
	solvemu(T, T_ref, mu_ref, Nx, Ny, mu);

	//Solve for speed of sound
	solveSOS(rho1, p, gamma, Nx, Ny, SOS);
	

	if(ip%checkStep == 0){
	    t2 = high_resolution_clock::now();
	    cout << "-------------------------------------------------" << endl;
	    cout << " Step = "<< ip << ", time = " << time << ", dt = " << dt << endl; 
	    cout << "-------------------------------------------------" << endl;
	    cout << " Time since last timestep = " << duration_cast<nanoseconds>(t2-t1).count()/(double)1000000000 << endl;
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

	    t1 = high_resolution_clock::now();
	}


	if(ip%100 == 0){
	    cout << "DUMPING FIELD" << endl;

    	    ofstream outfile;
	    string outputFileName;
    	    outputFileName = "rho.out.";
	    outputFileName.append(to_string(ip)); 
    	    outfile.open(outputFileName);
    	    for(int jp = 0; jp < Nx; jp++){
		for(int kp = 0; kp < Ny; kp++){
       	 	    outfile << rho1[jp*Ny+kp] << " ";
		}
		outfile << endl;
    	    }
	    outfile.close();


    	    outputFileName = "rhoU.out."; 
	    outputFileName.append(to_string(ip)); 
    	    outfile.open(outputFileName);
    	    for(int jp = 0; jp < Nx; jp++){
		for(int kp = 0; kp < Ny; kp++){
       	 	    outfile << rhoU1[jp*Ny+kp] << " ";
		}
		outfile << endl;
    	    }
	    outfile.close();

    	    outputFileName = "rhoV.out."; 
	    outputFileName.append(to_string(ip)); 
    	    outfile.open(outputFileName);
    	    for(int jp = 0; jp < Nx; jp++){
		for(int kp = 0; kp < Ny; kp++){
       	 	    outfile << rhoV1[jp*Ny+kp] << " ";
		}
		outfile << endl;
    	    }
	    outfile.close();

    	    outputFileName = "rhoE.out."; 
	    outputFileName.append(to_string(ip)); 
    	    outfile.open(outputFileName);
    	    for(int jp = 0; jp < Nx; jp++){
		for(int kp = 0; kp < Ny; kp++){
       	 	    outfile << rhoE1[jp*Ny+kp]*SOS[jp*Ny+kp] << " ";
		}
		outfile << endl;
    	    }
	    outfile.close();
	}


	

	if(time > timeEnd){
	    cout << "HIT END OF TIME" << endl;

    	    ofstream outfile;
	    string outputFileName;
    	    outputFileName = "rho.out"; 
    	    outfile.open(outputFileName);
    	    for(int jp = 0; jp < Nx; jp++){
		for(int kp = 0; kp < Ny; kp++){
       	 	    outfile << rho1[jp*Ny+kp] << " ";
		}
		outfile << endl;
    	    }
	    outfile.close();


    	    outputFileName = "p.out"; 
    	    outfile.open(outputFileName);
    	    for(int jp = 0; jp < Nx; jp++){
		for(int kp = 0; kp < Ny; kp++){
       	 	    outfile << p[jp*Ny+kp] << " ";
		}
		outfile << endl;
    	    }
	    outfile.close();

    	    outputFileName = "mach.out"; 
    	    outfile.open(outputFileName);
    	    for(int jp = 0; jp < Nx; jp++){
		for(int kp = 0; kp < Ny; kp++){
       	 	    outfile << U[jp*Ny+kp]*SOS[jp*Ny+kp] << " ";
		}
		outfile << endl;
    	    }
	    outfile.close();

	    break;
	}

    }
    //Garbage Collection
    delete[] mu;
    delete[] U;
    delete[] V;
    delete[] SOS;
    delete[] p;
    delete[] T;

    delete[] Ux;
    delete[] Uy;
    delete[] Vx;
    delete[] Vy;
    delete[] Tx;
    delete[] Ty;

    delete[] rho1;
    delete[] rhok;
    delete[] rho2;

    delete[] rhoU1;
    delete[] rhoUk;
    delete[] rhoU2;
 
    delete[] rhoV1;
    delete[] rhoVk;
    delete[] rhoV2;

    delete[] rhoE1;
    delete[] rhoEk;
    delete[] rhoE2;

    delete[] rho1k;
    delete[] rhoU1k;
    delete[] rhoV1k;
    delete[] rhoE1k;

    delete[] rhsDxIn;
    delete[] rhsDxOut;
    delete[] rhsDyIn;
    delete[] rhsDyOut;
 
*/
 
    return 0;
}

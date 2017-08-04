#include "Utils.hpp"
#include <chrono>
#include <ctime>
#include <iostream>
#include <fstream>
#include "Filter.hpp"
#include "Derivatives.hpp"
#include "IdealGas.hpp"

class Solver{

  public:

    //Class Objects that we'll need
    IdealGas *idealGas;
    Filter *filter;
    Derivatives *derivatives; 

    //Solver Options
    enum TimeStepping { CONST_DT, CONST_CFL};
    TimeStepping timeStepping;

    enum BoundaryConditions { PERIODIC, SPONGE, WALL, MOVING_WALL, INLET};
    BoundaryConditions bcX0, bcX1, bcY0, bcY1;
 
    //Mesh
    int Nx, Ny;
    double Lx, Ly;
    double *x, *y;
    double dx, dy;

    //Time and step stuff
    double dt, time, timeEnd;
    int step, checkStep, timeStep, filterCount, filterStep, outputStep;
    double CFL; bool endFlag;

    //Timing 
    std::chrono::system_clock::time_point t1, t2;
    
    //Initial Conditions
    double *rho0;
    double *U0;
    double *V0;
    double *p0;

    //Solver Data Arrays
    double *mu;
    double *U;
    double *V;
    double *SOS;
    double *p;
    double *T;
    double *Ux;
    double *Uy;
    double *Vx;
    double *Vy;
    double *Tx;
    double *Ty;    

    //Solver RK Step Data
    double *rho1;
    double *rhok;
    double *rho1k;
    double *rho2;

    double *rhoU1;
    double *rhoUk;
    double *rhoU1k;
    double *rhoU2;
    
    double *rhoV1;
    double *rhoVk;
    double *rhoV1k;
    double *rhoV2;

    double *rhoE1;
    double *rhoEk;
    double *rhoE1k;
    double *rhoE2;

    double *rhsDxIn; double *rhsDxOut;
    double *rhsDyIn; double *rhsDyOut;

    //Sponge Stuff
    double spongeAvgT;
    double spongeEpsP;
    double *spongeAvgRho;
    double *spongeAvgRhoU;
    double *spongeAvgRhoV;
    double *spongeAvgRhoE;
    double *spongeSigma;
    double spongeStrength;
    double spongeLengthX0;
    double spongeLengthX1;
    double spongeLengthY0;
    double spongeLengthY1;

    //Default Constructor
    Solver(){
        //Initial Conditions
        rho0 = NULL; U0   = NULL;
	V0   = NULL; p0   = NULL;

        //Solver Data Arrays
        mu  = NULL;
        U   = NULL; V = NULL;
        SOS = NULL;
        p   = NULL; T = NULL;
        Ux  = NULL; Uy = NULL;
        Vx  = NULL; Vy = NULL;
        Tx  = NULL; Ty = NULL;    

        //Solver RK Step Data
        rho1  = NULL; rhok = NULL;
        rho1k = NULL; rho2 = NULL;

        rhoU1  = NULL; rhoUk = NULL;
        rhoU1k = NULL; rhoU2 = NULL;
    
        rhoV1  = NULL; rhoVk = NULL;
        rhoV1k = NULL; rhoV2 = NULL;

        rhoE1  = NULL;  rhoEk = NULL;
    	rhoE1k = NULL;  rhoE2 = NULL;
	
	rhsDxIn = NULL; rhsDxOut = NULL;
	rhsDyIn = NULL; rhsDyOut = NULL;

	x = NULL; y = NULL;

	Nx = 0; Lx = 0;
	Ny = 0; Ly = 0;
	dx = 0; dy = 0;

	dt = 0; time = 0; timeEnd = 0;
	step = 0; checkStep = 0; timeStep = 0;
	filterCount = 0; filterStep = 0; outputStep = 0; CFL = 0;
	endFlag = 0;

	idealGas = NULL;
	filter   = NULL;
	derivatives = NULL;

	//Default to PERIODIC and CONST_CFL
	timeStepping = CONST_CFL;
	bcX0 = PERIODIC;
	bcX1 = PERIODIC;
	bcY0 = PERIODIC;
	bcY1 = PERIODIC;

        //Sponge Stuff
        spongeAvgT = 0;
        spongeEpsP = 0;
 	spongeStrength = 0;
        spongeLengthX0 = 0;
        spongeLengthX1 = 0;
        spongeLengthY0 = 0;
        spongeLengthY1 = 0;
       
	spongeAvgRho = NULL;
        spongeAvgRhoU = NULL;
        spongeAvgRhoV = NULL;
        spongeAvgRhoE = NULL;
        spongeSigma = NULL;
    }

    //Array initializing constructor
    Solver(int NX, int NY, double LX, double LY){
	
	//Allocate our array size containers
	Nx = NX; Lx = LX;
	Ny = NY; Ly = LY;

	x = new double[Nx];
	y = new double[Ny];

	for(int ip = 0; ip < Nx; ip++){
	    x[ip] = (((double)ip)/((double)Nx - 1.0))*Lx;
	}

	for(int jp = 0; jp < Ny; jp++){
	    y[jp] = (((double)jp)/((double)Ny - 1.0))*Ly;
	}

	dx = x[1]-x[0];
	dy = y[1]-y[0];

        //Initial Condition Arrays
        rho0 = new double[Nx*Ny]; 
        U0   = new double[Nx*Ny]; 
        V0   = new double[Nx*Ny]; 
        p0   = new double[Nx*Ny]; 

        //Solver Data Arrays 
        //Derived Data
        mu   = new double[Nx*Ny]; 
        U    = new double[Nx*Ny]; 
        V    = new double[Nx*Ny]; 
        SOS  = new double[Nx*Ny]; 
        p    = new double[Nx*Ny]; 
        T    = new double[Nx*Ny]; 
        Ux   = new double[Nx*Ny]; 
        Uy   = new double[Nx*Ny]; 
        Vx   = new double[Nx*Ny]; 
        Vy   = new double[Nx*Ny]; 
        Tx   = new double[Nx*Ny]; 
        Ty   = new double[Nx*Ny]; 

        //SolverStepData
        rho1   = new double[Nx*Ny]; 
        rhok   = new double[Nx*Ny]; 
        rho1k  = new double[Nx*Ny]; 
        rho2   = new double[Nx*Ny]; 

        rhoU1  = new double[Nx*Ny]; 
        rhoUk  = new double[Nx*Ny]; 
        rhoU1k = new double[Nx*Ny]; 
        rhoU2  = new double[Nx*Ny]; 

        rhoV1  = new double[Nx*Ny]; 
        rhoVk  = new double[Nx*Ny]; 
        rhoV1k = new double[Nx*Ny]; 
        rhoV2  = new double[Nx*Ny]; 

        rhoE1  = new double[Nx*Ny]; 
        rhoEk  = new double[Nx*Ny]; 
        rhoE1k = new double[Nx*Ny]; 
        rhoE2  = new double[Nx*Ny]; 

        rhsDxIn   = new double[Nx*Ny]; 
        rhsDxOut  = new double[Nx*Ny]; 
        rhsDyIn   = new double[Nx*Ny]; 
        rhsDyOut  = new double[Nx*Ny]; 

	rhsDxIn = new double[Nx*Ny]; 
	rhsDxOut = new double[Nx*Ny]; 
	rhsDyIn = new double[Nx*Ny]; 
	rhsDyOut = new double[Nx*Ny]; 

	dt = 0; time = 0; timeEnd = 0;
	step = 0; checkStep = 0; timeStep = 0;
	outputStep = 0; filterCount = 0; 
	filterStep = 0; CFL = 0;
	endFlag = 0;

	idealGas = new IdealGas(Nx, Ny);
	filter   = new Filter(Nx, Ny);
	derivatives = new Derivatives(Nx, Ny, dx, dy);

	//Default to PERIODIC and CONST_CFL
	timeStepping = CONST_CFL;
	bcX0 = PERIODIC;
	bcX1 = PERIODIC;
	bcY0 = PERIODIC;
	bcY1 = PERIODIC;

        t1 = std::chrono::system_clock::now();
        t2 = std::chrono::system_clock::now();

        //Sponge Stuff
        spongeAvgT = 0.0;
        spongeEpsP = 0.0;
 	spongeStrength = 0.0;
        spongeLengthX0 = 0;
        spongeLengthX1 = 0;
        spongeLengthY0 = 0;
        spongeLengthY1 = 0;
      
 
	spongeAvgRho = NULL;
        spongeAvgRhoU = NULL;
        spongeAvgRhoV = NULL;
        spongeAvgRhoE = NULL;
        spongeSigma = NULL;
        
    }


    void setBCForDerivatives();

    void initSpongeStuff();

    void calcSpongeSource(double *phi, double *spongeAvgPhi, double *spongeSource);

    void applyInitialCondition();
    
    void computeDtFromCFL_advanceTime();

    void computeCompactDY(double *phi, double *dphidy);

    void computeCompactDX(double *phi, double *dphidy);

    void computeVelocityTemperatureGradients();

    void computeContinuity(double *rhoU, double *rhoV);

    void computeXMomentum(double *rhoU, double *rhoV);

    void computeYMomentum(double *rhoU, double *rhoV);

    void computeEnergy(double *rhoE);

    void enforceBCs();

    void updateSolutionRKStep1();

    void updateSolutionRKStep2();

    void updateSolutionRKStep3();

    void updateSolutionRKStep4();

    void filterCompactY(double *phi, double *phiF);

    void filterCompactX(double *phi, double *phiF);

    void filterAndUpdateSolution();

    void updateEndOfStepPrimAndTemp();

    void checkSolution();

    void dumpSolution();

    void checkEnd();
};

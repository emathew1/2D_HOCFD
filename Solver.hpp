#include "Utils.hpp"
#include "Filter.hpp"
#include "Derivatives.hpp"
#include "IdealGas.hpp"

class Solver{

  public:

    //Class Objects that we'll need
    IdealGas *idealGas;
    Filter *filter;
    Derivatives *derivatives; 

    //Mesh
    int Nx, Ny;
    double Lx, Ly;
    double *x, *y;
    double dx, dy;

    //Time and step stuff
    double dt, time, timeEnd;
    int step, checkStep, timeStep, filterStep, outputStep;
    double CFL;
    
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
	filterStep = 0; outputStep = 0; CFL = 0;

	idealGas = NULL;
	filter   = NULL;
	derivatives = NULL;

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
	outputStep = 0; filterStep = 0; CFL = 0;

	idealGas = new IdealGas(Nx, Ny);
	filter   = new Filter(Nx, Ny);
	derivatives = new Derivatives(Nx, Ny, dx, dy);

    }

    void hellotest();

    void applyInitialCondition();
    
    void computeDtFromCFL_advanceTime();

    void computeVelocityTemperatureGradients();

    void computeContinuity();

    void computeXMomentum();

    void computeYMomentum();

    void computeEnergy();

};

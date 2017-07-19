
class Solver{

  public:

    //Size of Mesh
    int Nx, Ny;
    double Lx, Ly;

    //Time and step stuff
    double dt, time, timeEnd;
    int step, checkStep, timeStep;
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
    double *p
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
    double *rho1k
    double *rho2;

    double *rhoU;
    double *rhoU1;
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

        rhoU   = NULL; rhoU1 = NULL;
        rhoU1k = NULL; rhoU2 = NULL;
    
        rhoV1  = NULL; rhoVk = NULL;
        rhoV1k = NULL; rhoV2 = NULL;

        rhoE1  = NULL;  rhoEk = NULL;
    	rhoE1k = NULL;  rhoE2 = NULL;
    }

    //Array initializing constructor
    Solver(int NX, int NY){
	
	//Allocate our array size containers
	Nx = NX;
	Ny = NY;

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

    }


};

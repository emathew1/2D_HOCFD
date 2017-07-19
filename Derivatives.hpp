
class Derivatives{


    public:
	
	double alpha;
	double beta;
	double a, b, c;

	double *diagx, *offx;
	double *diagy, *offy;

    Derivatives(){
	alpha = 1.0/3.0;
	beta  = 0.0;
	a = 14.0/9.0;
	b = 1.0/9.0;
	c = 0.0;

	diagx = NULL; offx = NULL;
	diagy = NULL; offy = NULL;
    }

    Derivatives(int NX, int NY){
	alpha = 1.0/3.0;
	beta  = 0.0;
	a = 14.0/9.0;
	b = 1.0/9.0;
	c = 0.0;

	diagx = new double[NX]; 
	offx  = new double[NX];
	diagy = new double[NY]; 
	offy  = new double[NY];
    }

};

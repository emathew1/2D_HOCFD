
class Filter{

    public:
	double alphaF;
	double betaF;
	double a0, a1, a2, a3, a4, a5;


	double *diagFx, *offFx;
	double *diagFy, *offFy;

    Filter(){
	alphaF = 0.49;
	betaF  = 0.0;

	a0 = (93.0 + 70.0*alphaF)/128.0;
        a1     = ( 7.0 + 18.0*alphaF)/16.0;
        a2     = (-7.0 + 14.0*alphaF)/32.0;
        a3     = ( 1.0 -  2.0*alphaF)/16.0;
        a4     = (-1.0 +  2.0*alphaF)/128.0;
        a5     = 0.0;

	diagFx = NULL; offFx = NULL;
	diagFy = NULL; offFy = NULL;
    }


    Filter(int NX, int NY){
	alphaF = 0.49;
	betaF  = 0.0;

	a0 = (93.0 + 70.0*alphaF)/128.0;
        a1     = ( 7.0 + 18.0*alphaF)/16.0;
        a2     = (-7.0 + 14.0*alphaF)/32.0;
        a3     = ( 1.0 -  2.0*alphaF)/16.0;
        a4     = (-1.0 +  2.0*alphaF)/128.0;
        a5     = 0.0;

	diagFx = new double[NX]; 
	offFx  = new double[NX];
	diagFy = new double[NY]; 
	offFy  = new double[NY];
	
    }

};


class IdealGas{

    public:

	double gamma;
	double Pr;
	double R_gas;
	double cp;
	double p_ref;
	double rho_ref;
	double T_ref;
	double mu_ref;

    IdealGas(){
	gamma = 1.4;
	Pr    = 0.7;
	p_ref = 1.0/gamma;
	rho_ref = 1.0;
	T_ref = 1.0;
	mu_ref = 1.0;
	R_gas = p_ref/rho_ref/T_ref;
	cp = R_gas*gamma/(gamma - 1.0);

    }

};

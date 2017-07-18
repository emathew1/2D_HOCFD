#if !defined(_FUNC_H_)
#define _FUNC_H_

void solvemu(double *T, double T_ref, double mu_ref, int Nx, int Ny, double *mu);
void solveU(double *rho, double *rhoU, int Nx, int Ny, double *U);
void solvep(double *rho, double *rhoE, double *U, double *V, double gamma, int Nx, int Ny, double *p);
void solveT(double *rho, double *p, double R_gas, int Nx, int Ny, double *T);
void solvee(double *rho, double *p, double gamma, int Nx, int Ny, double *e);
void solveSOS(double *rho, double *p, double gamma, int Nx, int Ny, double *SOS);
void solveTri(double a[], double b[], double c[], double d[], double x[], double *work, int n);
void cyclic(double *a, double *b, double *c, double alpha, double beta,
	double *r, int n, double *x);
void multRHSDeriv(double a, double b, double c, double dh, double *phi, int N, double *RHSvec);
void CompactDYPeriodic(double *phi, double *diagy, double *offy, double alpha, double a, double b, double c, double dy, int Nx, int Ny, double *dphidy);
void transposeMatrix(double *in, int Nx, int Ny, double *out);
void transposeMatrix_Fast1(const double *in, int n, int p, double *out, int block);
void transposeMatrix_Fast2(const double *in, int n, int p, double *out, int blocksize);
void CompactDXPeriodic(double *phi, double *diagx, double *offx, double alpha, double a, double b, double c, double dx, int Nx, int Ny, double *dphidx);
void multRHSFilter(double a0, double a1, double a2, double a3, double a4, double a5, double *phi, int N, double *RHSvec);
void FilterPeriodicY(double *phi, double *diagFy, double *offFy, double alphaF, double a0, double a1, double a2, double a3, double a4, double a5, int Nx, int Ny, double *phiF);
void FilterPeriodicX(double *phi, double *diagFy, double *offFy, double alphaF, double a0, double a1, double a2, double a3, double a4, double a5, int Nx, int Ny, double *phiF);
void getRange(double *phi, std::string dataName, int Nx, int Ny);



#endif

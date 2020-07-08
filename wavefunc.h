//For ode subroutine
extern const int neqn;
extern double *y;
extern double *yp;
extern double *x;
extern double xout;
extern const double relerr;
extern const double abserr;
extern int *iflag;
extern double *work;
extern int *iwork;

//Properties of psi
extern double E;
extern double V;
extern double x_0;

extern void f ( double t, double y[], double yp[], double E, double V);
double glomin(double a, double b, double c, double m, double machep,
              double e, double t, double f(double x), double *x);
double local_min(double a, double b, double eps, double t,
                 double f(double x), double *x);
double local_min_rc(double *a, double *b, int *status, double value);
double r8_max(double x, double y);
double r8_sign(double x);
double zero ( double a, double b, double machep, double t, 
  double _f(void (*f)(double x, double *y, double *yp, double E, double V), double *x, double xout, double E, double V, double x_0, const double igrad, double *y, double *yp
        , const int neqn, double relerr, double abserr, int *iflag, double *work, int *iwork), void (*f)(double x, double *y, double *yp, double E, double V), 
        double *x, double xout, double V, double x_0, const double igrad, double *y, double *yp, 
        const int neqn, double relerr, double abserr, int *iflag, double *work, int *iwork);
void zero_rc(double a, double b, double t, double *arg, int *status,
             double value);

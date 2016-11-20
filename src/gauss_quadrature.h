#ifndef GAUSS_QUADRATURE_H
#define GAUSS_QUADRATURE_H

int gauss_quadrature(int order, int kind, double alpha, double beta, double a, double b, double x[], double w[]);
int gauss_quadrature_standard(int order, int kind, double alpha, double beta, double a, double b, double x[], double w[]);
int scale_gausspoints(int order, int kind, double alpha, double beta, double a, double b, double x[], double w[]);

void cdgqf ( int nt, int kind, double alpha, double beta, double t[], 
  double wts[] );
void cgqf_onlyscale ( int nt, int kind, double alpha, double beta, double a, double b, double t[], double wts[] );

void cgqf ( int nt, int kind, double alpha, double beta, double a, double b, 
  double t[], double wts[] );
double class_matrix ( int kind, int m, double alpha, double beta, double aj[], 
  double bj[] );
void imtqlx ( int n, double d[], double e[], double z[] );
void parchk ( int kind, int m, double alpha, double beta );
double r8_abs ( double x );
double r8_epsilon ( );
double r8_sign ( double x );
void scqf ( int nt, double t[], int mlt[], double wts[], int nwts, int ndx[], 
  double swts[], double st[], int kind, double alpha, double beta, double a, 
  double b );
void sgqf ( int nt, double aj[], double bj[], double zemu, double t[], 
  double wts[] );

#endif

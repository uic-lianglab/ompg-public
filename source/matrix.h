// matrix.h

#ifndef _MATRIX_
#define _MATRIX_

#include "atom.h"

//move from rootMeanSquare.cpp to here.
#define NR_END 1
#define NP 3
#define MP 500
#define FREE_ARG char*
#define MAXSTR 80
#define DEBUG 0

// from nr.h
double pythag(double a, double b);

void svdcmp(double **a, int m, int n, double w[], double **v);

// from nrutil.h
static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
static double maxarg1, maxarg2;
#define FMAX(a, b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))
static int iminarg1, iminarg2;
#define IMIN(a, b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))
#define IMAX(a, b) (iminarg1=(a),iminarg2=(b),(iminarg1) > (iminarg2) ? (iminarg1) : (iminarg2))

double *Vector(long nl, long nh);

void free_Vector(double *v, long nl, long nh);

double **matrix(long nrl, long nrh, long ncl, long nch);

void free_matrix(double **m, long nrl, long nrh, long ncl, long nch);

void nrerror(string error_text);

// from matrix.h
double **covarianceMatrix(double **x, double *x_mu, double **y, double *y_mu, int start, int end, int n);

double **subtract_matrix(double **x, double **y, int m, int n);

/* function prototypes */
double **var_matrix(double **x, double **x_mu, int m, int n);

double matrixNorm(double **x, double *x_mean, int start, int end, int n);

/* determines the mean Vector of a matrix */
double *mu_matrix(double **x, int start, int end, int n);

/* creates an identity matrix of size mxn and returns it */
double **identity_matrix(int m, int n);

double determinant_matrix_3x3(double **y);

/* creates an mirror identity matrix of size mxn and returns it */
double **mirror_identity_matrix(int m, int n);

double traceVector(double *x, int m);

int load_matrix(char *, int, double **);

int file_count(char file_name[]);

/* transpose matrix x of size m x n */
double **transpose_matrix(double **x, int m, int n);

double **multiply_matrix(double **x, int m, int n, double **y, int v, int b);

void clean_Vector(double *x, int m);

void print_matrix(double **a, int m, int n);

double rootMeanSquare(double **x, double **y, int start, int end);

void P_PT(Point p0, Point p1, double matrix[4][4]);

void RotMatrix(double **, double, double, double, double);


#endif

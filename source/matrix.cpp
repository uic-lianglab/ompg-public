// matrix.cpp

//---------------------------------------
// rootMeanSquare.cpp
//
#include <iostream>
#include <new>
#include <cstdio>
#include <cmath>
#include <cstdlib>

#define NRANSI

#include "matrix.h"

// from pythag.c
double pythag(double a, double b) {
    double absa, absb;
    absa = fabs(a);
    absb = fabs(b);
    if (absa > absb) return absa * sqrt(1.0 + SQR(absb / absa));
    else return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + SQR(absa / absb)));
}

// from svdcmp.c
void svdcmp(double **a, int m, int n, double w[], double **v) {
//  svd ,  a*diag(w)*v'=svd(a); a: 1:m * 1:n
    double pythag(double a, double b);
    int flag, i, its, j, jj, k, l, nm;
    double anorm, c, f, g, h, s, scale, x, y, z;
    double *rv1;

    rv1 = Vector(1, n);
    g = scale = anorm = 0.0;
    for (i = 1; i <= n; i++) {
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i <= m) {
            for (k = i; k <= m; k++) scale += fabs(a[k][i]);
            if (scale) {
                for (k = i; k <= m; k++) {
                    a[k][i] /= scale;
                    s += a[k][i] * a[k][i];
                }
                f = a[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][i] = f - g;
                for (j = l; j <= n; j++) {
                    for (s = 0.0, k = i; k <= m; k++) s += a[k][i] * a[k][j];
                    f = s / h;
                    for (k = i; k <= m; k++) a[k][j] += f * a[k][i];
                }
                for (k = i; k <= m; k++) a[k][i] *= scale;
            }
        }
        w[i] = scale * g;
        g = s = scale = 0.0;
        if (i <= m && i != n) {
            for (k = l; k <= n; k++) scale += fabs(a[i][k]);
            if (scale) {
                for (k = l; k <= n; k++) {
                    a[i][k] /= scale;
                    s += a[i][k] * a[i][k];
                }
                f = a[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][l] = f - g;
                for (k = l; k <= n; k++) rv1[k] = a[i][k] / h;
                for (j = l; j <= m; j++) {
                    for (s = 0.0, k = l; k <= n; k++) s += a[j][k] * a[i][k];
                    for (k = l; k <= n; k++) a[j][k] += s * rv1[k];
                }
                for (k = l; k <= n; k++) a[i][k] *= scale;
            }
        }
        anorm = FMAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
    }
    for (i = n; i >= 1; i--) {
        if (i < n) {
            if (g) {
                for (j = l; j <= n; j++)
                    v[j][i] = (a[i][j] / a[i][l]) / g;
                for (j = l; j <= n; j++) {
                    for (s = 0.0, k = l; k <= n; k++) s += a[i][k] * v[k][j];
                    for (k = l; k <= n; k++) v[k][j] += s * v[k][i];
                }
            }
            for (j = l; j <= n; j++) v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }
    for (i = IMIN(m, n); i >= 1; i--) {
        l = i + 1;
        g = w[i];
        for (j = l; j <= n; j++) a[i][j] = 0.0;
        if (g) {
            g = 1.0 / g;
            for (j = l; j <= n; j++) {
                for (s = 0.0, k = l; k <= m; k++) s += a[k][i] * a[k][j];
                f = (s / a[i][i]) * g;
                for (k = i; k <= m; k++) a[k][j] += f * a[k][i];
            }
            for (j = i; j <= m; j++) a[j][i] *= g;
        } else for (j = i; j <= m; j++) a[j][i] = 0.0;
        ++a[i][i];
    }
    for (k = n; k >= 1; k--) {
        for (its = 1; its <= 30; its++) {
            flag = 1;
            for (l = k; l >= 1; l--) {
                nm = l - 1;
                if ((double) (fabs(rv1[l]) + anorm) == anorm) {
                    flag = 0;
                    break;
                }
                if ((double) (fabs(w[nm]) + anorm) == anorm) break;
            }
            if (flag) {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++) {
                    f = s * rv1[i];
                    rv1[i] = c * rv1[i];
                    if ((double) (fabs(f) + anorm) == anorm) break;
                    g = w[i];
                    h = pythag(f, g);
                    w[i] = h;
                    h = 1.0 / h;
                    c = g * h;
                    s = -f * h;
                    for (j = 1; j <= m; j++) {
                        y = a[j][nm];
                        z = a[j][i];
                        a[j][nm] = y * c + z * s;
                        a[j][i] = z * c - y * s;
                    }
                }
            }
            z = w[k];
            if (l == k) {
                if (z < 0.0) {
                    w[k] = -z;
                    for (j = 1; j <= n; j++) v[j][k] = -v[j][k];
                }
                break;
            }
            if (its == 20) nrerror("no convergence in 20 svdcmp iterations");


            x = w[l];
            nm = k - 1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = pythag(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
            c = s = 1.0;
            for (j = l; j <= nm; j++) {
                i = j + 1;
                g = rv1[i];
                y = w[i];
                h = s * g;
                g = c * g;
                z = pythag(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y *= c;
                for (jj = 1; jj <= n; jj++) {
                    x = v[jj][j];
                    z = v[jj][i];
                    v[jj][j] = x * c + z * s;
                    v[jj][i] = z * c - x * s;
                }
                z = pythag(f, h);
                w[j] = z;
                if (z) {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = c * g + s * y;
                x = c * y - s * g;
                for (jj = 1; jj <= m; jj++) {
                    y = a[jj][j];
                    z = a[jj][i];
                    a[jj][j] = y * c + z * s;
                    a[jj][i] = z * c - y * s;
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
        }
    }
    free_Vector(rv1, 1, n);
}

// from nrutil.c
void nrerror(string error_text)
/* Numerical Recipes standard error handler */
{
    fprintf(stderr, "Numerical Recipes run-time error...\n");
    fprintf(stderr, "%s\n", error_text.c_str());
    fprintf(stderr, "...now exiting to system...\n");
    exit(1);
}

double *Vector(long nl, long nh)
/* allocate a double Vector with subscript range v[nl..nh] */
{
    double *v;
    v = new double[nh - nl + 1 + NR_END];
//	if (!v) std::cout<<"allocation failure in Vector()\n";
    return v;
}

int **matrixint(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    int **m;
    /* allocate pointers to rows */
    //        m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
    m = new int *[nrow + NR_END];
//	if (m==NULL) std::cout<<"allocation failure 1 in matrix()\n";
    //        m += NR_END;
    //        m -= nrl;

    /* allocate rows and set pointers to them */
    //        m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
    for (i = 0; i <= nrow; i++) {
        m[i] = new int[ncol + NR_END];
//		if (m[i]==NULL) std::cout<<"allocation failure 2 in matrix()\n";
    }
    //       m[nrl] += NR_END;
    //        m[nrl] -= ncl;
    //        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
    /* return pointer to array of pointers to rows */
    return m;
}

extern double **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    double **m;
    /* allocate pointers to rows */
    //        m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
    m = new double *[nrow + NR_END];
//	if (m==NULL) std::cout<<"allocation failure 1 in matrix()\n";
    //        m += NR_END;
    //        m -= nrl;

    /* allocate rows and set pointers to them */
    //        m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
    for (i = 0; i <= nrow; i++) {
        m[i] = new double[ncol + NR_END];
//		if (m[i]==NULL) std::cout<<"allocation failure 2 in matrix()\n";
    }
    //       m[nrl] += NR_END;
    //        m[nrl] -= ncl;
    //        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
    /* return pointer to array of pointers to rows */
    return m;
}

void free_Vector(double *v, long nl, long nh)
/* free a double Vector allocated with Vector() */
{
    delete[] v;

}

void free_matrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by matrix() */
{
    long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    for (i = 0; i <= nrow; ++i) {
        delete[] m[i];
    }
    delete[] m;

}

// from matrix.c
/**********************************************************************/
/* file_count()
/* determines the number of lines in a file                           */
/**********************************************************************/
int file_count(char file_name[MAXSTR]) {

    FILE *fp;           /* file pointer */
    int line_count = 0;
    char line[1500];



    /* Count number of lines in a file */
    if ((fp = fopen(file_name, "r")) == NULL) {
        printf("%s not found\n", file_name);
        nrerror("Data file not found\n");
    }

    while (!feof(fp)) {
        char *bob = fgets(line, 1500, fp);
        delete[] bob;
        line_count++;
    }
    fclose(fp);               /* close file pointer */
    line_count--;

    return line_count;
}

/**********************************************************************/
/* load_matrix()                                                      */
/* load in the matrix from file                                       */
/* returns the m dimension of the matrix, since n=3 for our purposes  */
/**********************************************************************/
int load_matrix(char file_name[MAXSTR], int atoms, double **x) {
    int k, l, n = NP;        /* counter variables */
    FILE *fp;           /* file pointer */

    /* reopen file and load into matrix */
    if ((fp = fopen(file_name, "r")) == NULL) {
        printf("%s not found\n", file_name);
        nrerror("Data file not found\n");
    }

    for (k = 1; k <= atoms; k++) {
        for (l = 1; l <= n; l++) {
            int bob = fscanf(fp, "%lf ", &x[k][l]);
        }
    }

    fclose(fp);               /* close file pointer */

    return 1;
}

double *vector_minus(double *r1, double *r2, int m) {
    int i;
    double *tempdif;
    tempdif = Vector(1, m);

    for (i = 1; i <= m; i++) {
        tempdif[i] = r1[i] - r2[i];
    }
    return tempdif;
}

double *vector_ctimes(double c1, double *r2, int m) {
    int i;
    double *tempdif;
    tempdif = Vector(1, m);

    for (i = 1; i <= m; i++) {
        tempdif[i] = c1 * r2[i];
    }
    return tempdif;
}

double *vector_xtimes(double *r1, double *r2, int m) {
    /*
    (Cross product)
    Vector1(x1,y1,z1) X Vector2(x2,y2,z2)=(ox,oy,oz)
    ox = (y1 * z2) - (y2 * z1)
    oy = (z1 * x2) - (z2 * x1)
    oz = (x1 * y2) - (x2 * y1)
    http://www.imagic3d.com/laohe1.htm
    	*/


    double *tempdif;
    tempdif = Vector(1, m);

    tempdif[1] = r1[2] * r2[3] - r1[3] * r2[2];
    tempdif[2] = r1[3] * r2[1] - r1[1] * r2[3];
    tempdif[3] = r1[1] * r2[2] - r1[2] * r2[1];

    return tempdif;
}


double vector_ptimes(double *r1, double *r2, int m) {
    double ptimes = 0;
    for (int i = 1; i <= m; i++) {
        ptimes += r1[i] * r2[i];
    }
    return ptimes;
}


double vector_dist(double *r1, double *r2, int m) {
    int i;
    double vdist = 0;
    for (i = 1; i <= m; i++) {
        vdist += (r1[i] - r2[i]) * (r1[i] - r2[i]);
    }
    return sqrt(vdist);
}

/**********************************************************************/
/* determinant_matrix_3x3()                                           */
/**********************************************************************/
double determinant_matrix_3x3(double **y) {

    double det = 0;

    det =
            (y[1][1] * y[2][2] * y[3][3]) +
            (y[1][3] * y[2][1] * y[3][2]) +
            (y[1][2] * y[2][3] * y[3][1]) -
            (y[1][3] * y[2][2] * y[3][1]) -
            (y[1][2] * y[2][1] * y[3][3]) -
            (y[1][1] * y[2][3] * y[3][2]);

    return det;
}

/**********************************************************************/
/* covariance_matrix()                                                */
/**********************************************************************/
double **covarianceMatrix(double **x, double *x_mu, double **y, double *y_mu, int start, int end, int n) {
    int k, l;
    double **tmp;
    double **tmpx;
    double **tmpy;
    double **tmpyT;

    int m = end - start + 1;

    tmpy = matrix(1, m, 1, NP);
    tmpx = matrix(1, m, 1, NP);
    /* precompute the mean matrix of x and y for speed */
    for (k = 1; k <= m; k++) {
        for (l = 1; l <= n; l++) {
            tmpy[k][l] = y[k + start - 1][l] - y_mu[l];  //printf("TMP Y\n");  print_matrix(tmpy,3,3);
            tmpx[k][l] = x[k + start - 1][l] - x_mu[l];  //printf("TMP X\n");  print_matrix(tmpx,3,3);
        }
    }
    tmpyT = transpose_matrix(tmpy, m, n);
    tmp = multiply_matrix(tmpyT, n, m, tmpx, m, n);
    for (k = 1; k <= n; k++) {
        for (l = 1; l <= n; l++)
            tmp[k][l] /= m - 1;  // Why m-1????
    }

    /* cleanup */
    free_matrix(tmpx, 1, m, 1, NP);
    free_matrix(tmpy, 1, m, 1, NP);
    free_matrix(tmpyT, 1, n, 1, m);
    //    if (DEBUG)  printf("\t\tcovariance_matrix()::Covariance Matrix\n"); print_matrix(tmp,3,3);

    return tmp;
}

/**********************************************************************/
/* mu_matrix()                                                        */
/* determines the mean Vector of a matrix                             */
/**********************************************************************/
double *mu_matrix(double **x, int start, int end, int n) {
    int k, l;        /* counter variables */
    double *tmp;
    int m = end - start + 1;

    tmp = Vector(1, n);
    clean_Vector(tmp, n);          /* Vector initialized to 2.364208 ? */

    for (k = 1; k < n; k++)
        tmp[k] = 0;


    for (k = 1; k <= m; k++) {
        for (l = 1; l <= n; l++) {
            tmp[l] += x[k + start - 1][l];
        }
    }

    for (l = 1; l <= n; l++) {          /* divide by m columns */
        tmp[l] /= m;
    }

    return tmp;
}
/**********************************************************************/
/* identity matrix()                                                  */
/* creates an identity matrix of size mxn and returns it              */
/**********************************************************************/
double **identity_matrix(int m, int n) {
    int k, l;                             /* counter variables */
    double **tmp;                         /* returned matrix */

    tmp = matrix(1, m, 1, NP);               /* initialize matrix */

    for (k = 1; k <= m; k++) {
        for (l = 1; l <= n; l++)
            if (k == l) {                    /* value 1 on diagonal */
                tmp[k][l] = 1;
            } else {
                tmp[k][l] = 0;
            }
    }
    return tmp;
}

/**********************************************************************/
/* identity matrix()                                                  */
/* creates an identity matrix of size mxn and returns it              */
/**********************************************************************/
double **mirror_identity_matrix(int m, int n) {
    int k, l;                             /* counter variables */
    double **tmp;                         /* returned matrix */

    tmp = matrix(1, m, 1, NP);               /* initialize matrix */

    for (k = 1; k <= m; k++) {
        for (l = 1; l <= n; l++)
            if (k == l) {                    /* value 1 on diagonal */
                tmp[k][l] = 1;
            } else {
                tmp[k][l] = 0;
            }
    }
    tmp[m][n] = -1;
    return tmp;
}
/**********************************************************************/
/* transpose_matrix()                                                 */
/* transposes a matrix of size mxn                                    */
/**********************************************************************/
double **transpose_matrix(double **x, int m, int n) {
    int k, l;
    double **tmp;

    tmp = matrix(1, n, 1, m);
    for (k = 1; k <= m; k++) {
        for (l = 1; l <= n; l++)
            tmp[l][k] = x[k][l];
    }
    return tmp;
}
/**********************************************************************/
/* multiply_matrix()                                                  */
/* multiply one matrix by another                                     */
/* careful when reusing this code, certain exceptions in variable
passing have been made with the knowledge that matrix 'y' will be
a transpose of a mx3 matrix and so it will be of size 3xm
(therefore the dimensions of y donot have to be passed             */
/**********************************************************************/
double **multiply_matrix(double **x, int m, int n, double **y, int v, int b) {
    int k, l, j;
    double **tmp;

    tmp = matrix(1, m, 1, b);
    if (n != v) nrerror("bad dimensions in mulitply_matrix()");
    for (k = 1; k <= m; k++)
        for (l = 1; l <= b; l++)
            tmp[k][l] = 0;

    for (k = 1; k <= m; k++) {     /* row of x */
        for (l = 1; l <= b; l++) {    /* column of x */

            for (j = 1; j <= v; j++) {
                tmp[k][l] += x[k][j] * y[j][l];
            }
        }
    }
    return tmp;
}

/**********************************************************************/
/* traceVector()                                                      */
/* compute trace for a Vector representing the diagonal of a matrix   */
/**********************************************************************/
double traceVector(double *x, int m) {

    int k;
    double tmp = 0;

    for (k = 1; k <= m; k++) {
        tmp += x[k];
    }
    if (DEBUG) {
        printf("\t\ttraceVector() :: sigma sqaure = %lf \n", tmp);
    }

    return tmp;
}
/**********************************************************************/
/* clean_Vector()                                                     */
/* initializes a numerical recipes Vector and initializes it to zero  */
/**********************************************************************/
void clean_Vector(double *x, int m) {
    int k;
    for (k = 1; k <= m; k++)
        x[k] = 0;
}
/**********************************************************************/
/* print_matrix()                                                     */
/* print matrix routine                                               */
/**********************************************************************/
void print_matrix(double **a, int m, int n) {
    int k, l;
    printf("[");
    for (k = 1; k <= m; k++) {
        for (l = 1; l <= n; l++)
            printf("%12.6f", a[k][l]);
        printf(";\n");
    }
    printf("]\n");
}

void print_vector(double *a, int m) {
    int k;

    for (k = 1; k <= m; k++) {
        printf("%12.6f", a[k]);
        printf("\n");
    }
}
/**********************************************************************/
/**********************************************************************/
/* THIS IS WRONG */  // made some  changes 09.05.2003
/**********************************************************************/
/**********************************************************************/
double matrixNorm(double **x, double *x_mu, int start, int end, int n) {
    int k, l;

    int m = end - start + 1;
    double tmp = 0;

    for (k = 1; k <= m; k++) {
        for (l = 1; l <= n; l++) {
            tmp += (x[k + start - 1][l] - x_mu[l]) * (x[k + start - 1][l] - x_mu[l]);
        }
    }

    tmp /= (m - 1);
    /* print out sigma squared */
    if (DEBUG) {
        printf("\t\tmatrixNorm() :: ||x||^2 = %lf \n", tmp);
    }

    return tmp;
}

double rootMeanSquare(double **x, double **y, int start, int end) {
    if (end - start == 1)
        return 0;
    /* counter variables */
    double *x_mu;
    double x_norm, y_norm;
    double *y_mu;

    /* determinatnt of u */
    double ds;                                   /* trace(DS) */
    double rms;                                  /* root mean square */
    double mse;                                  /* mean square error */

    double *w;
    double **u;
    double **v;
    /* svd matrices, Vectors */

    int m = end - start + 1;

    /* initialize least squares Vectors and matrices */

    w = Vector(1, NP);
    v = matrix(1, m, 1, NP);
    /* create mean Vector and matrix norms */
    x_mu = mu_matrix(x, start, end, 3);
    x_norm = matrixNorm(x, x_mu, start, end, 3);
    y_mu = mu_matrix(y, start, end, 3);
    y_norm = matrixNorm(y, y_mu, start, end, 3);
    /* create covariance */
    u = covarianceMatrix(x, x_mu, y, y_mu, start, end, 3);


    /* begin svd stuff */

    svdcmp(u, 3, 3, w, v);

    ds = traceVector(w, 3); // need to mulitply by S


    mse = x_norm + y_norm - (2 * ds);
    /* Rounding error zero rms correction */
    if (mse <= 0.000001) {
        rms = 0.0;
    } else {
        rms = (double) sqrt(mse);
    }
    free_matrix(u, 1, NP, 1, NP);
    free_matrix(v, 1, m, 1, NP);


    free_Vector(w, 1, NP);
    free_Vector(x_mu, 1, NP);
    free_Vector(y_mu, 1, NP);
    return rms;
}

/*---------------------------------------------
//above are from zhangjinfeng.
//following codes are added by chenyu for estructure.
---------------------------------------------*/
void fprintf_vectori(char file_name[MAXSTR], int atoms, int *x) {
    int k, n = NP;        /* counter variables */
    FILE *fp;           /* file pointer */

    /* reopen file and load into matrix */
    if ((fp = fopen(file_name, "w")) == NULL) {
        printf("%s error\n", file_name);
        nrerror("output file error\n");
    }

    for (k = 1; k <= atoms; k++) {
        fprintf(fp, "%d\t%d\n", k, x[k]);
    }
    fclose(fp);               /* close file pointer */

}

void fprintf_vector(char file_name[MAXSTR], int atoms, double *x) {
    int k, n = NP;        /* counter variables */
    FILE *fp;           /* file pointer */

    /* reopen file and load into matrix */
    if ((fp = fopen(file_name, "w")) == NULL) {
        printf("%s error\n", file_name);
        nrerror("output file error\n");
    }

    for (k = 1; k <= atoms; k++) {
        fprintf(fp, "%d\t %lf\n", k, x[k]);
    }
    fclose(fp);               /* close file pointer */

}

void fprintf_matrix(char file_name[MAXSTR], int atoms, double **x) {
    int k, l, n = NP;        /* counter variables */
    FILE *fp;           /* file pointer */

    /* reopen file and load into matrix */
    if ((fp = fopen(file_name, "w")) == NULL) {
        printf("%s fprintf_matrix error\n", file_name);
        nrerror("output file error\n");
    }

    for (k = 1; k <= atoms; k++) {
        for (l = 1; l <= n; l++) {
            fprintf(fp, "%lf\t", x[k][l]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);               /* close file pointer */

}

void fprintf_matrixint(char file_name[MAXSTR], int atoms, int ncol, int **x) {
    int k, l, n;        /* counter variables */
    FILE *fp;           /* file pointer */
    n = ncol;
    /* reopen file and load into matrix */
    if ((fp = fopen(file_name, "w")) == NULL) {
        printf("%s fprintf_matrixint error\n", file_name);
        nrerror("output file error\n");
    }
    for (k = 1; k <= atoms; k++) {
        for (l = 1; l <= n; l++) {
            fprintf(fp, "%d\t", x[k][l]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);               /* close file pointer */

}

void fprintm_matrix(char file_name[MAXSTR], int atoms, double **x) {
    int k, l, n = NP;        /* counter variables */
    FILE *fp;           /* file pointer */

    /* reopen file and load into matrix */
    if ((fp = fopen(file_name, "w")) == NULL) {
        printf("%s error\n", file_name);
        nrerror("output file error\n");
    }
    fprintf(fp, "[");

    for (k = 1; k <= atoms; k++) {
        for (l = 1; l <= n; l++) {
            fprintf(fp, "%lf\t", x[k][l]);
        }
        fprintf(fp, ";\n");
    }
    fprintf(fp, "]");
    fclose(fp);               /* close file pointer */

}


extern void linear_eq(double **A, double *y, double *x, int nume) {
    int i, j, k;
    double **B, temp;
    B = matrix(1, nume, 1, nume + 1);
    int zoline = 0;
    int zonum[5];

    for (k = 1, i = 1; i <= nume; i++) //row
    {
        if (y[i] == 0) {
            zoline += 1;
            zonum[zoline] = i;
            continue;


        }
        for (j = 1; j <= nume; j++) //column
            B[k][j] = A[i][j];
        B[k][nume + 1] = y[i];
        k += 1;
    }

    /*---------------------------------------------
    move all the 0 lines to the end.
    ---------------------------------------------*/

    for (i = 1; i <= zoline; i++) {
        for (j = 1; j <= nume; j++) //column
        {
            B[k][j] = A[zonum[i]][j];
        }
        B[k][nume + 1] = y[zonum[i]];
        k += 1;
    }

    /*---------------------------------------------
    transform B to up-triangle matirx
    ---------------------------------------------*/

    for (i = 2; i <= nume - zoline; i++) {
        for (j = i; j <= nume - zoline; j++) //line
        {
            temp = B[j][i - 1];
            for (k = i - 1; k <= nume + 1; k++) //column
            {
                B[j][k] = B[j][k] - B[i - 1][k] * temp / B[i - 1][i - 1];
            }
        }
    }

    /*---------------------------------------------
    calculate x[]
    ---------------------------------------------*/

    double tcl = 0.0;
    int cxl = nume - zoline;
    if (B[cxl][nume] == 0)
        x[nume] = 0;
    else {
        x[nume] = B[cxl][nume + 1] / B[cxl][nume];
        cxl -= 1;
    }

    for (k = nume - 1, i = nume - 1; i >= 1; i--, k--) //k the column
    {
        if (B[i][k] == 0) {
            x[k] = 0;
        }
        else {
            tcl = 0.0;
            for (j = nume; j >= k + 1; j--)
                tcl += B[i][j] * x[j];
            x[k] = (B[i][nume + 1] - tcl) / B[i][k];
        }
    }
    /*---------------------------------------------
    check x[]
    ---------------------------------------------*/
    double yht[5];
    for (i = 1; i <= nume; i++) //line
    {
        yht[i] = 0;
        for (j = 1; j <= nume; j++) //column
        {
            yht[i] += x[j] * A[i][j];
        }
    }

}

void addextention(char *filename, char extn1[MAXSTR], int extn2) {
    char extention[MAXSTR];

    sprintf(extention, "_%s_%d.txt", extn1, extn2);
    strcat(filename, extention);
}


void fprintpdb_matrix(char file_name[MAXSTR], int atoms, double **x, char cctype[20][10], int *sequ, double trate) {
    int k, l, n = NP;        /* counter variables */
    FILE *fp;           /* file pointer */

    /* reopen file and load into matrix */
    if ((fp = fopen(file_name, "w")) == NULL) {
        printf("%s error\n", file_name);
        nrerror("output file error\n");
    }
    for (k = 1; k <= atoms; k++) {
        fprintf(fp, "ATOM%7d  CA  %s A%4d    ", k, cctype[sequ[k]], k);
        for (l = 1; l <= n; l++) {
            fprintf(fp, "%8.3f", trate * x[k][l]);
        }
        fprintf(fp, "  1.0   0.0\n");
    }
    fclose(fp);               /* close file pointer */

}


void pvector_minus(double *tempdif, double *r1, double *r2, int m) {
    int i;

    for (i = 1; i <= m; i++) {
        tempdif[i] = r1[i] - r2[i];
    }
}

void pvector_ctimes(double *tempdif, double c1, double *r2, int m) {
    int i;


    for (i = 1; i <= m; i++) {
        tempdif[i] = c1 * r2[i];
    }

}

void pvector_xtimes(double *tempdif, double *r1, double *r2, int m) {
    tempdif[1] = r1[2] * r2[3] - r1[3] * r2[2];
    tempdif[2] = r1[3] * r2[1] - r1[1] * r2[3];
    tempdif[3] = r1[1] * r2[2] - r1[2] * r2[1];

}

void pcalculate_u(double *tempu, double *ur1, double *ur2, int m) {
    double determ;
    int i;
    determ = vector_dist(ur1, ur2, m);  // vector distance
    for (i = 1; i <= m; i++) {
        tempu[i] = (ur1[i] - ur2[i]) / determ;
    }

}

void P_PT(Point p0, Point p1, double matrix[4][4]) {
    matrix[1][1] = p0.x * p1.x;
    matrix[1][2] = p0.x * p1.y;
    matrix[1][3] = p0.x * p1.z;
    matrix[2][1] = p0.y * p1.x;
    matrix[2][2] = p0.y * p1.y;
    matrix[2][3] = p0.y * p1.z;
    matrix[3][1] = p0.z * p1.x;
    matrix[3][2] = p0.z * p1.y;
    matrix[3][3] = p0.z * p1.z;
};

void RotMatrix(double **m_matrix, double x, double y, double z, double theta) {
    m_matrix[1][1] = cos(theta) + (1 - cos(theta)) * x * x;
    m_matrix[1][2] = (1 - cos(theta)) * x * y - z * sin(theta);
    m_matrix[1][3] = (1 - cos(theta)) * x * z + y * sin(theta);
    m_matrix[2][1] = (1 - cos(theta)) * x * y + z * sin(theta);
    m_matrix[2][2] = cos(theta) + (1 - cos(theta)) * y * y;
    m_matrix[2][3] = (1 - cos(theta)) * y * z - x * sin(theta);
    m_matrix[3][1] = (1 - cos(theta)) * x * z - y * sin(theta);
    m_matrix[3][2] = (1 - cos(theta)) * y * z + x * sin(theta);
    m_matrix[3][3] = cos(theta) + (1 - cos(theta)) * z * z;
}

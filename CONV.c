#include "mex.h"
#include <stdio.h>
#include <math.h>
#include <complex.h>
typedef struct {
    double *real;
    double *imag;
} Complex;
Complex X,H,Y;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    X.real=mxGetPr(prhs[0]);
    X.imag=mxGetPi(prhs[0]);
    H.real=mxGetPr(prhs[1]);
    H.imag=mxGetPi(prhs[1]);
    int X_length=mxGetN(prhs[0]);
    int H_length=mxGetN(prhs[1]);
    int Y_length=X_length+H_length-1;
    Y=mxCreateDoubleMatrix(1,Y_length,mxCOMPLEX);
    Y.real=mxGetPr(plhs[0]);
    Y.imag=mxGetPi(plhs[0]);

    double complex *x=malloc(1*X_length*sizeof(double));
    double complex *h=malloc(1*H_length*sizeof(double));
    double complex *y=malloc(1*Y_length*sizeof(double));
        
    for(int i=0;i<X_length;i++){
        x[i] = X.real[i] + X.imag[i]*I;
     }
     
    for(int i=0;i<H_length;i++){
        h[i] = H.real[i] + H.imag[i]*I;
    }

    for(int i=0;i<X_length;i++){
        for(int j=0;j<H_length;j++){
            y[i+j]=x[i]*h[j];
        }
    }

    for(int i=0;i<Y_length;i++){
        Y.real[i]=creal(y[i]);
        Y.imag[i]=cimag(y[i]);
    }
    free(x);
    free(y);
    free(h);

}
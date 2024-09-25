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
    plhs[0]=mxCreateDoubleMatrix(1,Y_length,mxCOMPLEX);
    Y.real=mxGetPr(plhs[0]);
    Y.imag=mxGetPi(plhs[0]);

    double complex *x=malloc(X_length*sizeof(double complex));
    double complex *h=malloc(H_length*sizeof(double complex ));
    double complex *y=malloc(Y_length*sizeof(double complex ));
        
    for(int i=0;i<X_length;i++){
        if(mxIsComplex(prhs[0])){
            x[i] = X.real[i] + X.imag[i]*I;
        }
        else{
            x[i] = X.real[i];
        }
    }
     
     for(int i=0;i<H_length;i++){
         if(mxIsComplex(prhs[1])){
              h[i] = H.real[i] + H.imag[i]*I;
         }
         else{
              h[i] = H.real[i];
         }
     }
    for(int i=0;i<H_length;i++){
        for(int j=0;j<X_length;j++){
            y[i+j]+=x[j]*conj(h[H_length-i-1]);
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
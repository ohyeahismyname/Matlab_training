#include "mex.h"
#include <stdio.h>
#include <math.h>
#include <complex.h>

typedef struct{
    double *real;
    double *imag;
}Complex;
Complex A,Ans;

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[] ){
    
    A.real=mxGetPr(prhs[0]);       //取出數值
    A.imag=mxGetPi(prhs[0]);
    int col=mxGetN(prhs[0]);           //幾行幾列
    int row=mxGetM(prhs[0]);
    output=mxCreateDoubleMatrix(row,col,mxCOMPLEX);
    Ans.real=mxGetPr(plhs[0]); 
    Ans.imag=mxGetPi(plhs[0]); 
// 
//     double *decimal=mxGetPr(prhs[0]);
    int *QAM=mxGetPr(prhs[1]);
    x=(0:QAM-1);              // 0到M的數列
    y=qammod(x,QAM);
    Eavg = mean(abs(y).^2);     //average power
    NF = 1/sqrt(Eavg);
    for(int i=0;i<4;i++){
        A.real[i]=A.real[i]*NF;
        A.imag[i]=A.imag[i]*NF;
        Ans.real[i]=round( A.real[i] );
        Ans.imag[i]=round( A.imag[i] );
    }
}
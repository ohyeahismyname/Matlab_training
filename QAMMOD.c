#include "mex.h"
#include <stdio.h>
#include <math.h>
#include <complex.h>

typedef struct{
    double *real;
    double *imag;
}Complex;
Complex A;

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[] ){
//     int num_series = (int)mxGetScalar(prhs[0]);
    int num_QAM = (int)mxGetScalar(prhs[1]);
    int col = mxGetN(prhs[0]);
    int row = mxGetM(prhs[0]);
//     int N_out=nlhs;
//     gray code 輸出
//     double complex gray_real[4];
//     double complex gray_imag[4];

//     plhs[0]=mxCreateDoubleMatrix(1,4,mxCOMPLEX);
    
    A.real=mxGetPr(prhs[0]);
    

    if(num_QAM==4){
        double complex gray_real[col * row];
        double complex gray_imag[col * row];
        for(int i=0;i<col * row;i++){
            int cal=A.real[i];
            if((cal&2)==0) gray_real[i]+=(-1);
            else if((cal&2)==2) gray_real[i]+=(1);

            if((cal&1)==0) gray_imag[i]+=(1);
            else if ((cal&1)==1)gray_imag[i]+=(-1);
            
        }
        plhs[0]=mxCreateDoubleMatrix(col * row,1,mxCOMPLEX);
        double *real_output = mxGetPr(plhs[0]);
        double *imag_output = mxGetPi(plhs[0]);
    
        for (int i = 0; i < col * row; i++) {
            real_output[i] = gray_real[i];
            imag_output[i] = gray_imag[i];
        }
    }
    else if(num_QAM==16){
        double complex gray_real[col * row];
        double complex gray_imag[col * row];
        for(int i=0;i<col * row;i++){
            int cal=A.real[i];
            
            if((cal&12)==0) gray_real[i]+=(-3);
            else if((cal&12)==4) gray_real[i]+=(-1);
            else if((cal&12)==12) gray_real[i]+=(1);
            else if((cal&12)==8) gray_real[i]+=(3);

            if((cal&3)==0) gray_imag[i]+=(3);
            else if((cal&3)==1) gray_imag[i]+=(1);
            else if((cal&3)==3) gray_imag[i]+=(-1);
            else if((cal&3)==2) gray_imag[i]+=(-3);
        }
        plhs[0]=mxCreateDoubleMatrix(col * row,1,mxCOMPLEX);
        double *real_output = mxGetPr(plhs[0]);
        double *imag_output = mxGetPi(plhs[0]);
    
        for (int i = 0; i < col * row; i++) {
            real_output[i] = gray_real[i];
            imag_output[i] = gray_imag[i];
        }
    }
    else if(num_QAM==64){
        double complex gray_real[col * row];
        double complex gray_imag[col * row];
        for(int i=0;i<col * row;i++){
            int cal=A.real[i];
            if((cal&56)==0) gray_real[i]+=(-7);
            else if((cal&56)==8) gray_real[i]+=(-5);
            else if((cal&56)==24) gray_real[i]+=(-3);
            else if((cal&56)==16) gray_real[i]+=(-1);
            else if((cal&56)==48) gray_real[i]+=(1);
            else if((cal&56)==56) gray_real[i]+=(3);
            else if((cal&56)==40) gray_real[i]+=(5);
            else if((cal&56)==32) gray_real[i]+=(7);

            if((cal&7)==0) gray_imag[i]+=(7);
            else if((cal&7)==1) gray_imag[i]+=(5);
            else if((cal&7)==3) gray_imag[i]+=(3);
            else if((cal&7)==2) gray_imag[i]+=(1);
            else if((cal&7)==6) gray_imag[i]+=(-1);
            else if((cal&7)==7) gray_imag[i]+=(-3);
            else if((cal&7)==5) gray_imag[i]+=(-5);
            else if((cal&7)==4) gray_imag[i]+=(-7);
        }
        plhs[0]=mxCreateDoubleMatrix(col * row,1,mxCOMPLEX);
        double *real_output = mxGetPr(plhs[0]);
        double *imag_output = mxGetPi(plhs[0]);
    
        for (int i = 0; i < col * row; i++) {
            real_output[i] = gray_real[i];
            imag_output[i] = gray_imag[i];
        }
    }   
//     double *real_output = mxGetPr(plhs[0]);
//     double *imag_output = mxGetPi(plhs[0]);
// 
//     for (int i = 0; i < 4; i++) {
//         real_output[i] = gray_real[i];
//         imag_output[i] = gray_imag[i];
//     }


}
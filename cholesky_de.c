#include  "mex.h"
#include <complex.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
typedef struct{
    double *real;
    double *imag;
}Complex;
Complex A;
Complex Ainv;

void lowertriangle(int n,double complex M1[n][n],double complex temp[n][n]){//下三角
    for(int i=0;i<n;i++){               //row
        for(int j=0;j<n;j++){           //column
            double complex sum=0.000;
            if(i==j){//對角線
                M1[i][j];
                for(int k=0;k<j;k++){
                    sum += temp[i][k]*conj(temp[i][k]);
                }
                temp[i][j] = csqrt( M1[i][j] - sum );
            }
            else if(i > j){//非對角線元素
                for(int k=0;k<j;k++){
                    sum+=temp[i][k]*conj(temp[j][k]);
                }
                temp[i][j]=(M1[i][j]-sum)/temp[j][j];
            }
            else {
                temp[i][j] = 0.0;
            }
        }
    }
}

void gauss_jordan(int n,double complex trytry[n][n],double complex inverse[n][n]){//輸入的是下三角 所以只要往下做
      
     for (int i = 0; i < n; i++) { //生成單位矩陣
        for (int j = 0; j < n; j++) {
            if (i == j)
                inverse[i][j] = 1.0;
            else
                inverse[i][j] = 0.0;
        }
     }
    for(int i  = 0 ; i < n ; i++){//行
        for (int j = 0 ; j < n ; j++){//列
            if(trytry[j][i]==0){
                continue;
            }
            else if(i==j){
                double complex scale = trytry [j][j];
                for(int k=0;k<n;k++){
                    trytry[j][k]/=scale;
                    inverse[j][k]/=scale;
                }
            }
            else {
                double complex multi=trytry[j][i];
                for(int k=0;k<=i;k++){
                    inverse[j][k]-=inverse[i][k]*multi;
                }            
            }
        }    
    }
}
void Hermitian( int n,double complex inverse[n][n],double complex afterHER[n][n]){
    double complex CH[n][n];
    for (int i = 0; i < n; i++) { //先共軛
        for (int j = 0; j < n; j++) {
            CH[i][j]=conj(inverse[i][j]);
        }   
    }
    for (int i = 0; i < n; i++) { //後轉置
        for (int j = 0; j < n; j++) {
            afterHER[j][i] = CH[i][j];
        }
    }
}



void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[] ){//主程式
    int M,N;
    A.real=mxGetPr(prhs[0]);
    A.imag=mxGetPi(prhs[0]);
    M=mxGetM(prhs[0]);//number of rows
    N=mxGetN(prhs[0]);//number of columns
   
    plhs[0]=mxCreateDoubleMatrix(M,N,mxCOMPLEX);
    Ainv.real=mxGetPr(plhs[0]); 
    Ainv.imag=mxGetPi(plhs[0]);
    
    double complex M1[M][M];   
    for(int i=0;i<M;i++){
        for(int j=0;j<M;j++){
            if(mxIsComplex(prhs[0])){
                M1[i][j]=A.real[j*M+i] ;
                M1[i][j]+=A.imag[j*M+i]*I ;
            }
            else{
                M1[i][j]=A.real[j*M+i] ;
            }
        }
    }

                                                                                               
    double complex temp[M][M];   
    lowertriangle(M,M1,temp);
    
    double complex inverse[M][M];
    gauss_jordan(M,temp,inverse);   //求到反矩陣

    double complex afterHER[M][M];  //經過HERmitian
    Hermitian(M,inverse,afterHER);

    double complex Final[M][M] ;
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            double complex sum = 0.0;
            for (int k = 0; k < M; k++) {
                sum += afterHER[i][k] * inverse[k][j];
            }
            Final[i][j] = sum;
        }
    }
    for(int i=0;i<M;i++){
        for(int j=0;j<M;j++){
            if(mxIsComplex(prhs[0])){
                Ainv.real[j*M+i]=creal(Final[i][j]) ;
                Ainv.imag[j*M+i]=cimag(Final[i][j]) ;
            }
            else {
                Ainv.real[j*M+i]=creal(Final[i][j]) ; 
            }
        }
    }
    
}
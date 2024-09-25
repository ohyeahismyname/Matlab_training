#include "mex.h"
#include <complex.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "matrix.h"

typedef struct{
    double *real;
    double *imag;
}Complex;
Complex ZF,Y,H;

int Tx_num,M;

void Hermitian(double complex input[Tx_num][Tx_num],double complex output[Tx_num][Tx_num]){  //做賀密特
    double complex CH[Tx_num][Tx_num];
    for (int i = 0; i < Tx_num; i++) { //先共軛
        for (int j = 0; j < Tx_num; j++) {
            CH[i][j]=conj(input[i][j]);
        }
    }
    for (int i = 0; i < Tx_num; i++) { //後轉置
        for (int j = 0; j < Tx_num; j++) {
            output[j][i] = CH[i][j];
        }
    }
}
void multisingle(double complex input_h[Tx_num][Tx_num],double complex input_y[Tx_num],double complex output_zf[Tx_num]){//一維乘法
    for(int i=0;i<Tx_num;i++){
        output_zf[i] = 0;
    }
    for(int i=0;i<Tx_num;i++){
        for(int j= 0;j<Tx_num;j++){
            output_zf[i] += input_h[i][j]*input_y[j];
        }
    }
}

void multidouble(double complex input1[Tx_num][Tx_num],double complex input2[Tx_num][Tx_num],double complex output[Tx_num][Tx_num]){//二維乘法
    for(int i=0;i<Tx_num;i++){
        for(int j=0;j<Tx_num;j++){
            output[i][j] =0;
        }
    }
    for(int i=0;i<Tx_num;i++){
        for(int j=0;j<Tx_num;j++){
            for(int k=0;k<Tx_num;k++){
                output[i][j]+=input1[i][k]*input2[k][j];
            }
        }
    }
}

void lower(double complex l[M][M] , double complex mat[M][M]){
    double complex sum;
	for(int i = 0;i < M;i++){
		for(int j = 0;j < M;j++){
            sum = 0;
		    if( j == i){
				for(int k = 0; k < i; k++)
					sum += l[j][k] * conj( l[j][k]);
				l[j][i] = csqrt( mat[j][i] - sum );
		    }
		    else if(i < j){
		        for(int k = 0; k < i; k++)
					sum += l[j][k] * conj(l[i][k]);
				l[j][i] = (1/l[i][i]) * ( mat[j][i] - sum );
		    }
		}
    }
}
void invL(double complex input[M][M] , double complex output[M][M]){
    double complex temp[M][M];
    for(int i=0;i<M;i++){
        for(int j=0;j<M;j++){
            temp[i][j] = input[i][j];
            if(i == j)
                output[i][j] = 1;
            else
                output[i][j] = 0;
        }
    }
    for(int i=0;i<M;i++){
        for(int j=0;j<=i;j++){
            if( i != j){
                for(int k=0;k<=j;k++)
                    output[i][k] -= temp[i][j] * output[j][k];
                temp[i][j] = 0;
            }
            else{
                for(int k=0;k<=j;k++)
                    output[i][k] /= temp[i][j];
                temp[i][j] = 1;
            }
        }
    }
}
void M_conj(double complex input[M][M] , double complex output[M][M]){
    double complex H[M][M];
    for(int i=0;i<M;i++){
        for(int j=0;j<M;j++){
            output[i][j]=0;
			H[i][j]	= conj(input[j][i]);
        }
    }
    for(int i=0;i<M;i++)
        for(int j=0;j<M;j++)
            for(int k=0;k<M;k++)
                output[i][j] += H[i][k] * input[k][j];
}
void inv   (double complex input[Tx_num][Tx_num],double complex output[Tx_num][Tx_num]){
	double complex L[M][M];
    double complex l[M][M];
	lower(input,L);
	invL(L,l);
	M_conj(l,output);
}

int Ytrans(int carrier,int slot,int i){    //換算y的位置
    return carrier+1644*slot+1644*560*i;
}

int Htrans(int carrier,int slot,int i,int j){    //換算h的位置
    return carrier+1644*slot+1644*560*i+1644*560*4*j;
}


void finall(){
    double complex y[Tx_num];
	double complex h[Tx_num][Tx_num];
	double complex hH[Tx_num][Tx_num];
	double complex multi[Tx_num][Tx_num];
	double complex invh[Tx_num][Tx_num];
	double complex x[Tx_num];
	double complex p[Tx_num];
	
	for(int carrier=0;carrier<1644;carrier++){
		for(int slot=0;slot<560;slot++){
			for(int i=0;i<Tx_num;i++){      //寫入y
				y[i]= Y.real[Ytrans(carrier,slot,i)] + Y.imag[Ytrans(carrier,slot,i)]*I;
				for(int j=0;j<Tx_num;j++){  //寫入h
					h[i][j]	= H.real[Htrans(carrier,slot,i,j)] + H.imag[Htrans(carrier,slot,i,j)]*I;
                }
			}
			Hermitian(h,hH);
			multidouble(hH,h,multi);
			inv(invh,multi);
			multisingle(hH,y,p);
			multisingle(invh,p,x);
            
			for(int i=0;i<Tx_num;i++){
				ZF.real[Ytrans(carrier,slot,i)] = creal(x[i]);
				ZF.imag[Ytrans(carrier,slot,i)] = cimag(x[i]);
			}
		}
	}
}
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]){
	
	Y.real= mxGetPr( prhs[1] );     //Y接收的實部,虛部
	Y.imag= mxGetPi( prhs[1] );
	H.real= mxGetPr( prhs[0] );     //H通道的實部,虛部
	H.imag= mxGetPi( prhs[0] );
	Tx_num  = *mxGetPr(prhs[2]);    //輸出端的數量
	M = *mxGetPr(prhs[2]);
    
	mwSize M[3] = {1644,560,Tx_num};
    ZF.real=mxGetPr(plhs[0]);
	ZF.imag=mxGetPi(plhs[0]);
	plhs[0]=mxCreateNumericArray(3,M,mxDOUBLE_CLASS,mxCOMPLEX);
	
	finall();
    free(Y.real);
    free(Y.imag);
    free(H.real);
    free(H.imag);

}
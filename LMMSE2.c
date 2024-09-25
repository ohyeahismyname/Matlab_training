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
Complex LMMSE,Y,H;

int i,j,Tx_num,row,col,ram;
void herm   (double complex output[Tx_num][Tx_num],double complex input[Tx_num][Tx_num]){
    double complex temp[Tx_num][Tx_num];
    for ( i = 0; i < Tx_num; i++) { //先共軛
        for ( j = 0; j < Tx_num; j++) {
            temp[i][j]=conj(input[i][j]);
        }
    }
    for ( i = 0; i < Tx_num; i++) { //後轉置
        for ( j = 0; j < Tx_num; j++) {
            output[j][i] = temp[i][j];
        }
    }
}

void matrix_muilt(int Arow,int Acol,double complex A[Arow][Acol],int Brow,int Bcol,double complex B[Brow][Bcol],double complex Ans[Arow][Bcol]){
	int i,j,k;
	double complex ram;
	for(i=0;i<Arow;++i)
		for(k=0;k<Acol;++k){
			ram = A[i][k];
			for(j=0;j<Bcol;++j)
				Ans[i][j] += ram*B[k][j];
		}
}

void multip1(double complex output[Tx_num]    ,double complex input1[Tx_num][Tx_num],double complex input2[Tx_num]){
	for(i=0;i<Tx_num;i++)
	    output[i] = 0;
	for(i=0;i<Tx_num;i++)
	    for(ram = 0;ram<Tx_num;ram++)
            output[i] += input1[i][ram] * input2[ram];
}
void multip2(double complex output[Tx_num][Tx_num],double complex input1[Tx_num][Tx_num],double complex input2[Tx_num][Tx_num]){
    for(i=0;i<Tx_num;i++)
        for(j=0;j<Tx_num;j++)
			output[i][j] = 0;
	for(i=0;i<Tx_num;i++)
        for(j=0;j<Tx_num;j++)
            for(ram=0;ram<Tx_num;ram++)
                output[i][j] += input1[i][ram] * input2[ram][j];
}
void matadd (double complex output[Tx_num][Tx_num],double complex input1[Tx_num][Tx_num],double complex input2[Tx_num][Tx_num]){
	for(i=0;i<Tx_num;i++)
		for(j=0;j<Tx_num;j++)
			output[i][j] = input1[i][j] + input2[i][j];
}
void lower(double complex l[Tx_num][Tx_num] , double complex mat[Tx_num][Tx_num]){
    double complex sum_2conj;
	for(col = 0;col < Tx_num;col++){
		for(row = 0;row < Tx_num;row++){
		    if(row == col){
		        sum_2conj = 0;
				for(ram = 0; ram < col; ram++)
					sum_2conj += l[row][ram] * conj( l[row][ram]);
				l[row][col] = csqrt( mat[row][col] - sum_2conj );
		    }
		    else if(col < row){
		        sum_2conj = 0;
		        for(ram = 0; ram < col; ram++)
					sum_2conj += l[row][ram] * conj(l[col][ram]);
				l[row][col] = (1/l[col][col]) * ( mat[row][col] - sum_2conj );
		    }
		}
    }
}
void invt_L(double complex l[Tx_num][Tx_num] , double complex L[Tx_num][Tx_num]){
    double complex R[Tx_num][Tx_num];
    for(row=0;row<Tx_num;row++){
        for(col=0;col<Tx_num;col++){
            R[row][col] = L[row][col];
            if(row == col)
                l[row][col] = 1;
            else
                l[row][col] = 0;
        }
    }
    for(row=0;row<Tx_num;row++){
        for(col=0;col<=row;col++){
            if(row != col){
                for(ram=0;ram<=col;ram++)
                    l[row][ram] -= R[row][col] * l[col][ram];
                R[row][col] = 0;
            }
            else{
                for(ram=0;ram<=col;ram++)
                    l[row][ram] /= R[row][col];
                R[row][col] = 1;
            }
        }
    }
}
void M_conj(double complex a[Tx_num][Tx_num] , double complex l[Tx_num][Tx_num]){
    double complex H[Tx_num][Tx_num];
    for(row=0;row<Tx_num;row++){
        for(col=0;col<Tx_num;col++){
            a[row][col]=0;
			H[row][col]	= conj(l[col][row]);
        }
    }
    for(row=0;row<Tx_num;row++)
        for(col=0;col<Tx_num;col++)
            for(ram=0;ram<Tx_num;ram++)
                a[row][col] += H[row][ram] * l[ram][col];
}
void inv   (double complex output[Tx_num][Tx_num],double complex input[Tx_num][Tx_num]){
	double complex L[Tx_num][Tx_num];
    double complex l[Tx_num][Tx_num];
	lower(L,input);
	invt_L(l,L);
	M_conj(output,l);
}
int YPos(int carrier,int slot,int row){
	return carrier + 1644*slot + 1644*560*row;
}
int HPos(int carrier,int slot,int row,int col){
	return carrier + 1644*slot + 1644*560*row + 1644*560*Tx_num*col;
}
void total(){
	int carrier,slot;
	int i,j;
	double complex eye[Tx_num][Tx_num];
	double complex h[Tx_num][Tx_num];
	double complex hH[Tx_num][Tx_num];
	double complex mul[Tx_num][Tx_num];
	double complex inh[Tx_num][Tx_num];
	double complex y[Tx_num]	  ;
	double complex x[Tx_num]	  ;
	double complex temp[Tx_num]	  ;
	double complex te2[Tx_num][Tx_num];
	//init No
	for(i=0;i<Tx_num;i++)
		for(j=0;j<Tx_num;j++){
			if(i==j)
				eye[i][j] = 1;
			else
				eye[i][j] = 0;
		}
			
	//start calcul
	for(carrier=0;carrier<1644;carrier++){
		for(slot=0;slot<560;slot++){
			for(i=0;i<Tx_num;i++){
				y[i]		= Y.real[YPos(carrier,slot,i  )] + I*Y.imag[YPos(carrier,slot,i)];
				for(j=0;j<Tx_num;j++)
					h[i][j]	= H.real[HPos(carrier,slot,i,j)] + I*H.imag[HPos(carrier,slot,i,j)];
			}
			herm(hH	,h);
			multip2	(mul,hH	,h);
			matadd(te2,mul,eye);
			inv	(inh,te2);
			multip1 (temp,hH,y);
			multip1 ( x ,inh,temp);
			for(i=0;i<Tx_num;i++){
				LMMSE.real[YPos(carrier,slot,i)] = creal(x[i]);
				LMMSE.imag[YPos(carrier,slot,i)] = cimag(x[i]);
			}
		}
	}
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]){
	Y.real= mxGetPr(  prhs[0]  );
	Y.imag= mxGetPi(  prhs[0]  );
	H.real= mxGetPr(  prhs[1]  );
	H.imag= mxGetPi(  prhs[1]  );
	Tx_num= *mxGetPr( prhs[2]  );
	mwSize Matrix[3] = {1644,560,Tx_num};
    LMMSE.real = mxGetPr( plhs[0] );
	LMMSE.imag = mxGetPi( plhs[0] );
	prhs[2]= mxCreateNumericArray(3,Matrix,mxDOUBLE_CLASS,mxCOMPLEX);
	
	total();
}
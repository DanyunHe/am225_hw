#include <cmath>
#include <cstdio>

#include "blas.h"


/** Standard matrix matrix multiplication.
 *\ param[in] n matrix dimension.
 *\ Calculate C=AB. */
void std_mul(int n,double* A, double* B, double* C){
	
	double temp;
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			temp=0.;
			for(int k=0;k<n;k++){
				temp+=A[k*n+i]*B[j*n+k];
			}
			C[j*n+i]=temp;
		}
	}

}

/** Strassen's algorithm matrix matrix multiplication.
 *\ param[in] n matrix dimension.
 *\ Calculate C=AB. */
void strassen_mul(int n, double* A,double* B,double* C){

	if(n==1){
		C[0]=A[0]*B[0];
	}

	// Construct submatrix
	int nn=n/2;
	const int nn2=nn*nn;
	double *A00=new double[nn2],*A01=new double[nn2],*A10=new double[nn2],*A11=new double[nn2];
	double *B00=new double[nn2],*B01=new double[nn2],*B10=new double[nn2],*B11=new double[nn2];
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			if(i<nn){
				if(j<nn){
					A00[j*nn+i]=A[j*n+i];
					B00[j*nn+i]=B[j*n+i];
				}
				else{
					A01[(j-nn)*nn+i]=A[j*n+i];
					B01[(j-nn)*nn+i]=B[j*n+i];
				}
			}
			else if(j<nn){
				A10[j*nn+i-nn]=A[i*n+i];
				B10[j*nn+i-nn]=B[i*n+i];

			}
			else{
				A11[(j-nn)*nn+i-nn]=A[j*n+i];
				B11[(j-nn)*nn+i-nn]=B[j*n+i];
			}
		}
	}

	// Calculate Q
	double *Q0=new double[nn2],*Q1=new double[nn2],*Q2=new double[nn2],*Q3=new double[nn2],\
		*Q4=new double[nn2],*Q5=new double[nn2],*Q6=new double[nn2];

	double *temp1=new double[nn2],*temp2=new double[nn2];
	
	// Q0
	for(int i=0;i<nn;i++){
		for(int j=0;j<nn;j++){
			temp1[j*nn+i]=A00[j*nn+i]+A11[j*nn+i];
			temp2[j*nn+i]=B00[j*nn+i]+B11[j*nn+i];
		}
	}
	strassen_mul(nn,temp1,temp2,Q0);

	// Q1
	for(int i=0;i<nn;i++){
		for(int j=0;j<nn;j++){
			temp1[j*nn+i]=A10[j*nn+i]+A11[j*nn+i];
		}
	}
	strassen_mul(nn,temp1,B00,Q1);

	// Q2
	for(int i=0;i<nn;i++){
		for(int j=0;j<nn;j++){
			temp1[j*nn+i]=B01[j*nn+i]-B11[j*nn+i];
		}
	}
	strassen_mul(nn,A00,temp1,Q2);

	// Q3
	for(int i=0;i<nn;i++){
		for(int j=0;j<nn;j++){
			temp1[j*nn+i]=B10[j*nn+i]-B00[j*nn+i];
		}
	}
	strassen_mul(nn,A11,temp1,Q3);


	// Q4
	for(int i=0;i<nn;i++){
		for(int j=0;j<nn;j++){
			temp1[j*nn+i]=A00[j*nn+i]+A01[j*nn+i];
		}
	}
	strassen_mul(nn,temp1,B11,Q4);

	// Q5
	for(int i=0;i<nn;i++){
		for(int j=0;j<nn;j++){
			temp1[j*nn+i]=A10[j*nn+i]-A00[j*nn+i];
			temp2[j*nn+i]=B00[j*nn+i]+B01[j*nn+i];
		}
	}
	strassen_mul(nn,temp1,temp2,Q5);

	// Q6
	for(int i=0;i<nn;i++){
		for(int j=0;j<nn;j++){
			temp1[j*nn+i]=A01[j*nn+i]-A11[j*nn+i];
			temp1[j*nn+i]=B10[j*nn+i]+B11[j*nn+i];
		}
	}
	strassen_mul(nn,temp1,temp2,Q6);

	// Calculate C
	int idx;
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			if(i<nn){
				if(j<nn){
					idx=j*nn+i;
					C[j*n+i]=Q0[idx]+Q3[idx]-Q4[idx]+Q6[idx];
				}
				else{
					idx=(j-nn)*nn+i;
					C[j*n+i]=Q1[idx]+Q3[idx];
				}
			}
			else if(j<nn){
				idx=j*nn+i-nn;
				C[j+n+i]=Q2[idx]+Q4[idx];
			}
			else{
				idx=(j-nn)*nn+i-nn;
				C[j*n+i]=Q0[idx]+Q2[idx]-Q1[idx]+Q5[idx];
			}
		}
	}

}

/** BLAS matrix matrix multiplication.
 *\ param[in] n matrix dimension.
 *\ Calculate C=AB. */
void blas_mul(int n,double* A,double* B,double* C){
	char trans='n';
	double alpha=1.,beta=0.;

	dgemm_(&trans,&trans,&n,&n,&n,&alpha,A,&n,B,&n,&beta,C,&n);

}









#include <TF1.h>            // questo serve per la funzione cout
#include "algebra_func.h"

using namespace std;

void num_random_multivariata_3par(double* mu, double* sigma, double* cov_matrix, int N, double* rnd0, double* rnd1, double* rnd2, int Nrand){

  const double PI = 3.141592653589793;

  double Q2 = 0, fR2, fR1, max_normal;
  double det_cov, sqrt_det;

  double *inv_cov_matrix=(double*)malloc(sizeof(double)*N*N);
  double *min=(double*)malloc(sizeof(double)*N), *max=(double*)malloc(sizeof(double)*N);
  double *R=(double*)malloc(sizeof(double)*N);

  LU_invert(inv_cov_matrix,cov_matrix,N);
  
  det_cov = determinant(cov_matrix, N);
  printf("DET = %.20f\n", det_cov );
  sqrt_det = sqrt(det_cov);
  
  for(int i = 0; i < N; i++){
    
    min[i] = mu[i]-10*sigma[i];
    max[i] = mu[i]+10*sigma[i];
  }

  max_normal = 1/( pow( sqrt(2*PI), N)*sqrt_det);
  
  srand(time(NULL));
  
  for(int n = 0; n < Nrand; n++){
    
    do{
      
      for(int i = 0; i < N; i++){
	
	R[i] = min[i] + (max[i] - min[i])*( ( (double) rand() )/RAND_MAX );	
	
      }// i
      
      fR2 = max_normal*( ( (double) rand() )/RAND_MAX );
      
      for(int i = 0; i < N; i++){
	for(int j = 0; j < N; j++){      
	  
	  Q2 = Q2 + ( R[i] - mu[i])*inv_cov_matrix[i*N+j]*( R[j] - mu[j]);
	  
	}// j
      }// i
      
      fR1 = max_normal*exp( -Q2/2. );
      
      Q2 = 0;
      
    }while(fR2 > fR1);
    
    rnd0[n] = R[0];
    rnd1[n] = R[1];
    rnd2[n] = R[2];
    
  }// n
  
  free(inv_cov_matrix);
  free(max);
  free(min);
  free(R);

}// num_random_multivariata_3par


void num_random_multivariata_4par(double* mu, double* sigma, double* cov_matrix, int N, double* rnd0, double* rnd1, double* rnd2, double* rnd3, int Nrand){

  const double PI = 3.141592653589793;

  double Q2 = 0, fR2, fR1, max_normal;
  double det_cov, sqrt_det;

  double *inv_cov_matrix=(double*)malloc(sizeof(double)*N*N);
  double *min=(double*)malloc(sizeof(double)*N), *max=(double*)malloc(sizeof(double)*N);
  double *R=(double*)malloc(sizeof(double)*N);

  LU_invert(inv_cov_matrix,cov_matrix,N);
  
  det_cov = determinant(cov_matrix, N);
  printf("DET = %.20f\n", det_cov );
  sqrt_det = sqrt(det_cov);
  
  for(int i = 0; i < N; i++){
    
    min[i] = mu[i]-10*sigma[i];
    max[i] = mu[i]+10*sigma[i];
  }

  max_normal = 1/( pow( sqrt(2*PI), N)*sqrt_det);
  
  srand(time(NULL));
  
  for(int n = 0; n < Nrand; n++){
    
    do{
      
      for(int i = 0; i < N; i++){
	
	R[i] = min[i] + (max[i] - min[i])*( ( (double) rand() )/RAND_MAX );	
	
      }// i
      
      fR2 = max_normal*( ( (double) rand() )/RAND_MAX );
      
      for(int i = 0; i < N; i++){
	for(int j = 0; j < N; j++){      
	  
	  Q2 = Q2 + ( R[i] - mu[i])*inv_cov_matrix[i*N+j]*( R[j] - mu[j]);
	  
	}// j
      }// i
      
      fR1 = max_normal*exp( -Q2/2. );
      
      Q2 = 0;
      
    }while(fR2 > fR1);
    
    rnd0[n] = R[0];
    rnd1[n] = R[1];
    rnd2[n] = R[2];
    rnd3[n] = R[3];
    
  }// n
  
  free(inv_cov_matrix);
  free(max);
  free(min);
  free(R);

}// num_random_multivariata_4par


void num_random_multivariata_5par(double* mu, double* sigma, double* cov_matrix, int N, double* rnd0, double* rnd1, double* rnd2, double* rnd3, double* rnd4, int Nrand){

  const double PI = 3.141592653589793;

  double Q2 = 0, fR2, fR1, max_normal;
  double det_cov, sqrt_det;

  double *inv_cov_matrix=(double*)malloc(sizeof(double)*N*N);
  double *min=(double*)malloc(sizeof(double)*N), *max=(double*)malloc(sizeof(double)*N);
  double *R=(double*)malloc(sizeof(double)*N);

  LU_invert(inv_cov_matrix,cov_matrix,N);
  
  det_cov = determinant(cov_matrix, N);
  printf("DET = %.20f\n", det_cov );
  sqrt_det = sqrt(det_cov);
  
  for(int i = 0; i < N; i++){
    
    min[i] = mu[i]-10*sigma[i];
    max[i] = mu[i]+10*sigma[i];
  }

  max_normal = 1/( pow( sqrt(2*PI), N)*sqrt_det);
  
  srand(time(NULL));
  
  for(int n = 0; n < Nrand; n++){
    
    do{
      
      for(int i = 0; i < N; i++){
	
	R[i] = min[i] + (max[i] - min[i])*( ( (double) rand() )/RAND_MAX );	
	
      }// i
      
      fR2 = max_normal*( ( (double) rand() )/RAND_MAX );
      
      for(int i = 0; i < N; i++){
	for(int j = 0; j < N; j++){      
	  
	  Q2 = Q2 + ( R[i] - mu[i])*inv_cov_matrix[i*N+j]*( R[j] - mu[j]);
	  
	}// j
      }// i
      
      fR1 = max_normal*exp( -Q2/2. );
      
      Q2 = 0;
      
    }while(fR2 > fR1);
    
    rnd0[n] = R[0];
    rnd1[n] = R[1];
    rnd2[n] = R[2];
    rnd3[n] = R[3];
    rnd4[n] = R[4];
    
  }// n
  
  free(inv_cov_matrix);
  free(max);
  free(min);
  free(R);

}// num_random_multivariata_5par


void num_random_multivariata_3par_bound(double* mu, double* sigma, double* cov_matrix, int N, double* rnd0, double* rnd1, double* rnd2, double* bound_min, double* bound_max, int Nrand){

  const double PI = 3.141592653589793;

  double Q2 = 0, fR2, fR1, max_normal;
  double det_cov, sqrt_det;

  double *inv_cov_matrix=(double*)malloc(sizeof(double)*N*N);
  double *min=(double*)malloc(sizeof(double)*N), *max=(double*)malloc(sizeof(double)*N);
  double *R=(double*)malloc(sizeof(double)*N);

  LU_invert(inv_cov_matrix,cov_matrix,N);
  
  det_cov = determinant(cov_matrix, N);
  printf("DET = %.20f\n", det_cov );
  sqrt_det = sqrt(det_cov);
  
  for(int i = 0; i < N; i++){
    
    min[i] = mu[i]-10*sigma[i];
    max[i] = mu[i]+10*sigma[i];
  }

  max_normal = 1/( pow( sqrt(2*PI), N)*sqrt_det);
  
  srand(time(NULL));
  
  for(int n = 0; n < Nrand; n++){
    
    do{
      
      for(int i = 0; i < N; i++){
	
	R[i] = min[i] + (max[i] - min[i])*( ( (double) rand() )/RAND_MAX );	
	
      }// i
      
      fR2 = max_normal*( ( (double) rand() )/RAND_MAX );
      
      for(int i = 0; i < N; i++){
	for(int j = 0; j < N; j++){      
	  
	  Q2 = Q2 + ( R[i] - mu[i])*inv_cov_matrix[i*N+j]*( R[j] - mu[j]);
	  
	}// j
      }// i
      
      fR1 = max_normal*exp( -Q2/2. );
      
      Q2 = 0;
      
    }while( (fR2 > fR1) || (R[0] > bound_max[0]) || (R[1] > bound_max[1]) || (R[2] > bound_max[2]) || (R[0] < bound_min[0]) || (R[1] < bound_min[1]) || (R[2] < bound_min[2]) );
    
    rnd0[n] = R[0];
    rnd1[n] = R[1];
    rnd2[n] = R[2];
    
  }// n
  
  free(inv_cov_matrix);
  free(max);
  free(min);
  free(R);

}// num_random_multivariata_3par_bound

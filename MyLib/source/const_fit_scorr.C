#include <stdio.h>
#include <TF1.h>

void constant_fit_scorrelated( double *fitted_variable, double *variable, double *sigma_variable, int clusterfile, int Nt, int ti_fit, int tf_fit){

  double *num=(double*)malloc(sizeof(double)), *den=(double*)malloc(sizeof(double));
  
  num[0] = 0;
  den[0] = 0;
    
  for(int ijk = 0; ijk < clusterfile; ijk++){
    for(int t = ti_fit; t <= tf_fit; t++){
      
      num[0] = num[0] + ( variable[ijk*Nt+t]*(1/pow(sigma_variable[t], 2)) );
      den[0] = den[0] + ( 1/pow(sigma_variable[t], 2) );

    }// t
    fitted_variable[ijk] = num[0]/den[0];
    
    num[0] = 0;
    den[0] = 0;
    
  }// ijk

  free(num);
  free(den);
  
}// matrix_el_constant_fit

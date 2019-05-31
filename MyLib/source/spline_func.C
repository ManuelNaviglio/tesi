#include <TF1.h>            // questo serve per la funzione cout
#include "lettura_ult_inp.h"
#include "definitions.h"

using namespace std;

void spline2(int N, double x[], double y[], double B[], double C[]){
  
  double delta_r, delta_l, delta_1, delta_2, slope_r, slope_l, slope_1, slope_2;
  
  if(N == 1){
    B[0] = 0;
    C[0] = 0;
  }
  
  if(N == 2){
    B[0] = (y[1]-y[0])/(x[1]-x[0]);
    B[1] = B[0];
    C[0] = 0;
    C[1] = 0;
  }
  
  if(N > 2){
    
    for(int i = 1; i < N-1; i++){
      delta_r = x[i+1] - x[i];
      delta_l = x[i] - x[i-1];
      
      slope_r = (y[i+1] - y[i])/delta_r;
      slope_l = (y[i] - y[i-1])/delta_l;
      
      B[i] = ( slope_r*delta_l + slope_l*delta_r)/(delta_l + delta_r);
      C[i] = ( slope_r-slope_l )/(delta_l + delta_r);
    }
    
    //// punto iniziale
    
    delta_1 = x[1] - x[0];
    delta_2 = x[2] - x[0];
    
    slope_1 = (y[1] - y[0])/delta_1;
    slope_2 = (y[2] - y[0])/delta_2;

    B[0] = (slope_1*delta_2 - slope_2*delta_1)/(delta_2 - delta_1);
    C[0] = (slope_2 - slope_1)/(delta_2 - delta_1);
    
    //// punto finale
    
    delta_1 = x[N-2] - x[N-1];
    delta_2 = x[N-3] - x[N-1];
    
    slope_1 = (y[N-2] - y[N-1])/delta_1;
    slope_2 = (y[N-3] - y[N-1])/delta_2;
    
    B[N-1] = (slope_1*delta_2 - slope_2*delta_1)/(delta_2 - delta_1);
    C[N-1] = (slope_2 - slope_1)/(delta_2 - delta_1);
    
  }
  
}


double yspline2(int N, double x[], double y[], double B[], double C[], double x0){
  
  int NL = 0, NU = N-1;
  double dx, yspline_2;
  
  for(int i = 0; i < N; i++){
    
    if( x0 < x[N-1-i]){
      NU = N-1-i;
    }
    if( x0 > x[i] ){
      NL = i;
    }
    
  }
  
  if( NU == 0 ){
    dx = x0 - x[0];
    yspline_2 = y[0] + dx*( B[0] + dx*C[0] );
  }
  if( NL == N-1 ){
    dx = x0 - x[N-1];
    yspline_2 = y[N-1] + dx*( B[N-1] + dx*C[N-1] );
  }
  if( (NU != 0) && (NL != N-1) && ((x0-x[NL]) <= (x[NU]-x0)) ){  // nota che ho messo <= per tenere in conto la possibilità che x0 stia esattamente in mezzo
    
    dx = x0 - x[NL];
    yspline_2 = y[NL] + dx*( B[NL] + dx*C[NL] );
    
  }
  if( (NU != 0) && (NL != N-1) && ((x0-x[NL]) > (x[NU]-x0)) ){
    
    dx = x0 - x[NU];
    yspline_2 = y[NU] + dx*( B[NU] + dx*C[NU] );
    
  }
  
  return yspline_2;
  
}


void bootstrap_sampling( double* array, char* file_open, int ibeta, int imusea, int clusterfile){					       
  
  ////////  LETTURA Ultimate Input
  
  double mlight[Nev+1], mstrange[Nev+1], mcharm[Nev+1], a[Nbeta][Nev+1], ainv[Nbeta][Nev+1], r0[Nev+1], Zev[Nbeta][Nev+1], ZTev[Nbeta][Nev+1], f0[Nev+1], B0[Nev+1], fkfpi[Nev+1];
  int  iboot[Nbeta][Nmusea][Nev+1];
  
  lettura_ultimate_input( mlight, mstrange, mcharm, a, ainv, r0, Zev, iboot, f0, B0, fkfpi, ZTev);
  
  ////////  FINE LETTURA Ultimate Input
  
  double* temp = (double*)malloc(sizeof(double)*(clusterfile));
  FILE *fin;
  
  if ((fin = fopen(file_open, "r")) == NULL ){
    printf("Error opening the input file %s!!\n",file_open);
    exit(EXIT_FAILURE);
  }
  
  for(int ijk = 0; ijk < clusterfile; ijk++ ){
    fscanf(fin, "%lf\n", &temp[ijk]);
  }
  fclose(fin);
  
  for(int iev = 0; iev <= Nev; iev++ ){
    
    if(iev != 0){
      array[iev] = temp[iboot[ibeta][imusea][iev]-1];
    }
    if(iev == 0){
      array[iev] = temp[15];
    }
    
  }// iev
    
  free(temp);
  
}

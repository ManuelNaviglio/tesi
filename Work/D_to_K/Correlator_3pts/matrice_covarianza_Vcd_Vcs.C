#include <TF1.h>            // questo serve per la funzione cout
#include <fstream>
#include "stat_analysis_func.h"
#include "definitions.h"

#define LEN_NAME 1024
#define PI 3.141592653589793
#define GF 0.000011664

const int Nq = 1;
const int Nmatrix = 2;
const int Nanalisys = 64;
const int Nev_tot = 6400;
int analisys_in = 1, analisys_fin = 64;

//////////////////////////////

int energia_da_sinhDR = 0;
int energia_da_stdDR = 1;
int energia_dal_fit = 0;

//////////////////////////////

//////////////////////////////

int f0_with_S = 1;
int f0_without_S = 0;

//////////////////////////////

int main(){

  //////////////////////////////
  
  string dir_E[3] = {"E_from_sinhDR", "E_from_stdDR", "E_from_fit"};
  int iE;
  
  if(energia_da_sinhDR == 1){
    iE = 0;
  }
  if(energia_da_stdDR == 1){
    iE = 1;
  }
  if(energia_dal_fit == 1){
    iE = 2;
  }
  //////////////////////////////
  
  //////////////////////////////
  
  string dir_S[2] = {"with_S", "without_S"};
  int iS;
  
  if(f0_with_S == 1){
    iS = 0;
  }
  if(f0_without_S == 1){
    iS = 1;
  }
  
  //////////////////////////////

  
  //////////////// LEGGO Vcd e Vcs EVENTO PER EVENTO
  
  FILE *fin_Vcd, *fin_Vcs;   
  char file_Vcd[LEN_NAME], file_Vcs[LEN_NAME];
  
  double Vcd_boot[Nev_tot+1], Vcs_boot[Nev_tot+1], av_Vcd=0, av_Vcs=0;
  
  sprintf(file_Vcd, "OUTPUT_SMEAR/%s/%s/CKM/Vcd_corr.out", dir_S[iS].c_str(), dir_E[iE].c_str());
  sprintf(file_Vcs, "OUTPUT_SMEAR/%s/%s/CKM/Vcs_corr.out", dir_S[iS].c_str(), dir_E[iE].c_str());
  
  if ((fin_Vcd = fopen( file_Vcd, "r")) == NULL ){
    printf("Error opening the output file: %s\n",file_Vcd);
    exit(EXIT_FAILURE);
  }
  if ((fin_Vcs = fopen( file_Vcs, "r")) == NULL ){
    printf("Error opening the output file: %s\n",file_Vcs);
    exit(EXIT_FAILURE);
  }
  
  for(int iev = 1; iev <= Nev_tot; iev++){
    
    fscanf(fin_Vcd,"%lf\n", &Vcd_boot[iev]);
    fscanf(fin_Vcs,"%lf\n", &Vcs_boot[iev]);
    
    av_Vcd = av_Vcd + Vcd_boot[iev]/Nev_tot;
    av_Vcs = av_Vcs + Vcs_boot[iev]/Nev_tot;
    
  }// iev
  
  Vcd_boot[0] = av_Vcd;
  Vcs_boot[0] = av_Vcs;
  
  fclose(fin_Vcd);
  fclose(fin_Vcs);

  //////////////// FINE LEGGO Vcd e Vcs EVENTO PER EVENTO

  /////////// CALCOLO LA MATRICE DI COVARIANZA 

  double cov_matrix[Nmatrix][Nmatrix];
  
  for(int i = 0; i < Nmatrix; i++){
    for(int j = i; j < Nmatrix; j++){

     
      if( (i < Nq) && (j < Nq) ){

	cov_matrix[i][j] = cov_bootstrap( Vcd_boot, Vcd_boot, analisys_in, analisys_fin, Nev_an, clusterfile);
	
      }// if( (i < Nq) && (j < Nq) )


      if( (i < Nq) && (j >= Nq) ){

	cov_matrix[i][j] = cov_bootstrap( Vcd_boot, Vcs_boot, analisys_in, analisys_fin, Nev_an, clusterfile);
	
      }// if( (i < Nq) && (j >= Nq) )


      if( i >= Nq ){

	cov_matrix[i][j] = cov_bootstrap( Vcs_boot, Vcs_boot, analisys_in, analisys_fin, Nev_an, clusterfile);

      }// if( i >= Nq )      

      cov_matrix[j][i] = cov_matrix[i][j];
      
    }// j
  }// i

  /////////// FINE CALCOLO LA MATRICE DI COVARIANZA

  ////////////// CALCOLO LA VARIANZA
  
  double sigma[Nmatrix];
  
  for(int i = 0; i < Nmatrix; i++){
    
    if( i < Nq ){
      
      sigma[i] = sigma_bootstrap( Vcd_boot, analisys_in, analisys_fin, Nev_an, clusterfile);
      
    }// if( i < Nq )

    
    if( i >= Nq ){
      
      sigma[i] = sigma_bootstrap( Vcs_boot, analisys_in, analisys_fin, Nev_an, clusterfile);
      
    }// if( i < Nq )
    

    ///// CHECK COVARIANZA = VARIANZA SULLA DIAGONALE
    
    printf("VVVV i=%d cov=%f var=%f diff=%f\n", i, cov_matrix[i][i], pow( sigma[i], 2), cov_matrix[i][i]-pow( sigma[i], 2));
    
    ///// FINE CHECK COVARIANZA = VARIANZA SULLA DIAGONALE
    
  }// i
  
  //////////// FINE CALCOLO LA VARIANZA
  

  ////////// PRINT OUTPUT PER MATHEMATICA
  
  printf("OOOO\n\n");
  
  printf("{");
  
  for(int i = 0; i < Nmatrix; i++){
    
    printf("{");
    
    for(int j = 0; j < Nmatrix; j++){
      
      printf("%.6e", cov_matrix[i][j]);
      
      if( j != Nmatrix-1 ){
	
	printf(",");
	
      }
      
    }// j
    
    if(i != Nmatrix-1){
      
      printf("},\n");
      
    }else if(i == Nmatrix-1){
      
      printf("}");
	
    }
    
  }// i
  
  printf("}");
  
  
  ////////// CHECK
  printf("\nCCCC\n");
  printf("%.10f %.10f %.10f", cov_matrix[0][0], cov_matrix[0][1], cov_matrix[1][1]);
  
  
  ////////// FINE CHECK
  
  ////////// FINE PRINT OUTPUT PER MATHEMATICA  
  
  return 0;
  
}


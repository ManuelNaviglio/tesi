#include <TF1.h>            // questo serve per la funzione cout
#include <fstream>
#include "stat_analysis_func.h"
#include "algebra_func.h"

using namespace std;

void read_synt_form_factor( double *ff_synt, double *sigma_ff_synt, int Nq, int Nev_synt, int Nev_an, int analysis_in, int analysis_fin, double q2_max, string *dir_S, string *dir_E, string *dir_ff, int clusterfile, int iS, int iE, int iff){
  
  FILE *f_input;
  int LEN_NAME = 1024;
  char *name_file_synt=(char*)malloc(sizeof(char)*(LEN_NAME));
  
  double *q2_step=(double*)malloc(sizeof(double)), *av_ff=(double*)calloc( 1, sizeof(double)), *temp_ff_synt=(double*)malloc(sizeof(double)*(Nev_synt+1));
    
  for(int iq = 0; iq < Nq; iq++){
    
    sprintf(name_file_synt, "OUTPUT_SMEAR/%s/%s/dati_sintetici/%s/%s_sint_q2_%d.out", dir_S[iS].c_str(), dir_E[iE].c_str(), dir_ff[iff].c_str(), dir_ff[iff].c_str(), iq);
    
    if ((f_input = fopen( name_file_synt, "r")) == NULL ){
      printf("Error opening the output file %s\n",name_file_synt);
      exit(EXIT_FAILURE);
    }
    
    for(int iev = 1; iev <= Nev_synt; iev++){
      
      fscanf(f_input,"%lf\n", &ff_synt[iq*(Nev_synt+1)+iev]);

      *av_ff += ff_synt[iq*(Nev_synt+1)+iev]/Nev_synt;

      temp_ff_synt[iev] = ff_synt[iq*(Nev_synt+1)+iev];
      
    }// iev

    ff_synt[iq*(Nev_synt+1)+0] = *av_ff;

    sigma_ff_synt[iq] = sigma_bootstrap( temp_ff_synt, analysis_in, analysis_fin, Nev_an, clusterfile);
    
    fclose(f_input);

    *av_ff = 0;
    
  } // iq

  printf("#%s\n", dir_ff[iff].c_str());
  printf("@type xydy\n");

  for(int iq = 0; iq < Nq; iq++){

    *q2_step = (q2_max/( (double)  (Nq - 1)))*iq;
    
    printf("%f  %f  %f\n", *q2_step, ff_synt[iq*(Nev_synt+1)+0], sigma_ff_synt[iq]);
    
  }// iq

  printf("\n");

  free(temp_ff_synt);
  free(name_file_synt);
  free(av_ff);
  free(q2_step);
  
} //read_synt_form_factor


void build_ffs_ratio(double* Ratio, double* sigma_Ratio, double* ff_1, double* ff_2, double q2_max, int Nq, int Nev_synt, int Nev_ffs, int Nev_an, string *dir_ff, int analysis_in, int analysis_fin, int clusterfile, int iff1, int iff2, int check_partial_blocks){
  
  double *q2_step=(double*)malloc(sizeof(double)), *av_Ratio=(double*)calloc( 1, sizeof(double)), *av_Ratio_test=(double*)calloc( Nq, sizeof(double)), *temp_Ratio=(double*)malloc(sizeof(double)*(Nev_synt+1));
  
  int num_analysis = (analysis_fin - analysis_in) + 1;
  int start_ev = (analysis_in - 1)*Nev_an + 1, end_ev = analysis_fin*Nev_an;
  int Nev_real = num_analysis*Nev_an;
  
  
  for(int iq = 0; iq < Nq; iq++){
    
    for(int iev = 1; iev <= Nev_synt; iev++){
      
      Ratio[iq*(Nev_synt+1)+iev] = ff_1[iq*(Nev_ffs+1)+iev]/ff_2[iq*(Nev_ffs+1)+iev];
      
      *av_Ratio += Ratio[iq*(Nev_synt+1)+iev]/Nev_synt;

      temp_Ratio[iev] = Ratio[iq*(Nev_synt+1)+iev];

    }// iev

    ////// Check per blocchi di analisi generici
    for(int iev = start_ev; iev <= end_ev; iev++){
      
      av_Ratio_test[iq] += Ratio[iq*(Nev_synt+1)+iev]/Nev_real;

    }// iev
    /////////////

    Ratio[iq*(Nev_synt+1)+0] = *av_Ratio;
    temp_Ratio[0] = *av_Ratio;
    
    sigma_Ratio[iq] = sigma_bootstrap( temp_Ratio, analysis_in, analysis_fin, Nev_an, clusterfile);

    *av_Ratio = 0;
    
  }// iq
  
  if( ((iff1 == 0) && (iff2 == 1)) || ((iff1 == 1) && (iff2 == 0))){ // f0(0)/f+(0) = 1 esattamente
    sigma_Ratio[0] = 0.0;
  }
  
  
  if(check_partial_blocks == 0){
  
    printf("#Ratio %s/%s\n", dir_ff[iff1].c_str(), dir_ff[iff2].c_str());
    printf("@type xydy\n");

    for(int iq = 0; iq < Nq; iq++){
      
      *q2_step = (q2_max/( (double)  (Nq - 1)))*iq;
      
      printf("%f  %f  %f\n", *q2_step, Ratio[iq*(Nev_synt+1)+0], sigma_Ratio[iq]);
      
    }// iq
  
  }else if(check_partial_blocks == 1){

    printf("#Check for Ratio %s/%s, analysis %d and %d\n", dir_ff[iff1].c_str(), dir_ff[iff2].c_str(), analysis_in, analysis_fin);
    printf("@type xydy\n");

    for(int iq = 0; iq < Nq; iq++){
      
      *q2_step = (q2_max/( (double)  (Nq - 1)))*iq;
      
      printf("%f  %f  %f\n", *q2_step, av_Ratio_test[iq], sigma_Ratio[iq]);
      
    }// iq

  }// check_partial_blocks

  printf("\n");
  
  free(temp_Ratio);
  free(av_Ratio);
  free(av_Ratio_test);
  free(q2_step);
  
} // build_ffs_ratio


void covariance_syntetic_data(double* cov_, int Nq, int Nev_synt, int Nev_an, double var_, double* synt_data, int analysis_in, int analysis_fin, int clusterfile, int check){

  double *inv_cov_=(double*)malloc(sizeof(double)*(Nq)*(Nq));
  double *identity_=(double*)calloc( (Nq)*(Nq), sizeof(double));

  double **temp_synt_data=(double**)malloc(sizeof(double)*(Nq));
  for(int iq = 0; iq < Nq; iq++) temp_synt_data[iq]=(double*)malloc(sizeof(double)*(Nev_synt+1));

  for(int iq = 0; iq < Nq; iq++){
    for(int iev = 0; iev <= Nev_synt; iev++){
      temp_synt_data[iq][iev] = synt_data[iq*(Nev_synt+1)+iev];
    }// iev
  }// iq

  for(int i = 0; i < Nq; i++){
    for(int j = i; j < Nq; j++){
      
      cov_[i*Nq+j] = cov_bootstrap( temp_synt_data[i], temp_synt_data[j], analysis_in, analysis_fin, Nev_an, clusterfile);
      
      cov_[j*Nq+i] = cov_[i*Nq+j];
      
    }// j
  }// i

  for(int i = 0; i < Nq; i++){
    for(int j = 0; j < Nq; j++){
      if( i == j ){
				cov_[i*Nq+j] += var_;
      }// if      
    }// j
  }// i

  printf("\n------------------------------------------------------------\n");

  printf("Covariance matrix\n");
  print_matrix( Nq, Nq, cov_, 1);
  
  if(check == 1){

		LU_invert( inv_cov_, cov_, Nq);

		for(int i = 0; i < Nq; i++){
			for(int j = 0; j < Nq; j++){
				for(int k = 0; k < Nq; k++){
					identity_[i*Nq+k] += cov_[i*Nq+j]*inv_cov_[j*Nq+k];
				}// k
			}// j
		}// i
    
    printf("\nInverse of Covariance matrix\n");
    print_matrix_for_mathematica( Nq, Nq, inv_cov_, 0);
    
    printf("\nCheck identity after the inversion\n");
    print_matrix( Nq, Nq, identity_, 0);
    
  }// check
    
  printf("------------------------------------------------------------\n");
  
  free(temp_synt_data);
  free(identity_);
  free(inv_cov_);
  
}// covariance_syntetic_data

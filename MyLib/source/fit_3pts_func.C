#include <TF1.h>            // questo serve per la funzione cout
#include <fstream>
#include "fit_2pts_func.h"
#include "clean_str_cl.h"
#include "stat_analysis_func.h"
#include "definitions.h"

using namespace std;

void read_corr_3pts( double *corr_3pts, int Nt, int clusterfile, double parity, string *av, const string *beta_V, const string mu_sea_1[Nbeta][Nmusea], string *reim, const string n_conf[Nbeta][Nmusea], const string *mu_sea_2, const string *m_, const string *th_, const string *sme_, string *op, int ibeta, int imusea, int im1, int im2, int ith1, int ith2, int ismear, int iav, int iop, int ireim){

  int LEN_NAME = 1024;
  
  FILE *fr;
  char *open_file_name=(char*)malloc(sizeof(char)*LEN_NAME);
  double *temp_corr=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  
  if( (iop != 1) && (iop != 3) ){ // OCCHIO CHE QUESTO PREVEDE CHE op[1]="Vi" E op[3]="Ti"
    parity = 1;
}

sprintf(open_file_name, "%s/%s/%s/3pts.%s.%s.%s.m1m2_%s%s.th_%s%s.sme_30_%s0.%sP5.out", av[iav].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), reim[ireim].c_str(), n_conf[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith1].c_str(), th_[ith2].c_str(), sme_[ismear].c_str(), op[iop].c_str());

if ((fr = fopen(open_file_name, "r")) == NULL ){
  printf("Error opening the file to read: %s\n", open_file_name);
  exit(EXIT_FAILURE);
}

clean_string_cluster(fr, temp_corr, Nt, clusterfile);

fclose(fr);

for(int ijk = 0; ijk < clusterfile; ijk++){
  for(int t = 0; t < Nt; t++){

    corr_3pts[ijk*Nt+t] = parity*temp_corr[ijk*Nt+t];

    }// t
  }// ijk

  free(open_file_name);
  free(temp_corr);
  
}// read_corr_3pts


void read_corr_2pts( double *corr_2pts, int Nt, int clusterfile, const string *beta_V, const string mu_sea_1[Nbeta][Nmusea], const string n_conf[Nbeta][Nmusea], const string *mu_sea_2, const string *m_, const string *th_, const string *sme_, int ibeta, int imusea, int im1, int im2, int ith, int ismear){

  int LEN_NAME = 1024;
  
  FILE *fr_2pts_th_r1r2_00, *fr_2pts_th_r1r2_11;  
  char *open_2pts_th_r1r2_00=(char*)malloc(sizeof(char)*LEN_NAME), *open_2pts_th_r1r2_11=(char*)malloc(sizeof(char)*LEN_NAME);
  
  double *temp_2pts_th_r1r2_00=(double*)malloc(sizeof(double)*(clusterfile)*(Nt)), *temp_2pts_th_r1r2_11=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  double *temp_corr_2pts=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  
  if(ith != 3){

    FILE *fr_2pts_th_rev_r1r2_00, *fr_2pts_th_rev_r1r2_11;
    char *open_2pts_th_rev_r1r2_00=(char*)malloc(sizeof(char)*LEN_NAME), *open_2pts_th_rev_r1r2_11=(char*)malloc(sizeof(char)*LEN_NAME);
    double *temp_2pts_th_rev_r1r2_00=(double*)malloc(sizeof(double)*(clusterfile)*(Nt)), *temp_2pts_th_rev_r1r2_11=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
    
    // CASO ith != 3
    sprintf(open_2pts_th_r1r2_00, "../Correlator_2pts/DATA_SMEAR/%s/%s/2pts.%s.%s.m1m2_%s%s.th_%s.sme_30_%s0.r1r2_00.P5P5.out", beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), n_conf[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith].c_str(), sme_[ismear].c_str());
    sprintf(open_2pts_th_r1r2_11, "../Correlator_2pts/DATA_SMEAR/%s/%s/2pts.%s.%s.m1m2_%s%s.th_%s.sme_30_%s0.r1r2_11.P5P5.out", beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), n_conf[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith].c_str(), sme_[ismear].c_str());
    sprintf(open_2pts_th_rev_r1r2_00, "../Correlator_2pts/DATA_SMEAR/%s/%s/2pts.%s.%s.m1m2_%s%s.th_%s.sme_30_%s0.r1r2_00.P5P5.out", beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), n_conf[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[6-ith].c_str(), sme_[ismear].c_str());
    sprintf(open_2pts_th_rev_r1r2_11, "../Correlator_2pts/DATA_SMEAR/%s/%s/2pts.%s.%s.m1m2_%s%s.th_%s.sme_30_%s0.r1r2_11.P5P5.out", beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), n_conf[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[6-ith].c_str(), sme_[ismear].c_str());
    
    printf("READING: DATA_SMEAR/%s/%s/2pts.%s.%s.m1m2_%s%s.th_%s.sme_30_%s0.r1r2_00.P5P5.out\n", beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), n_conf[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith].c_str(), sme_[ismear].c_str());
    printf("READING: DATA_SMEAR/%s/%s/2pts.%s.%s.m1m2_%s%s.th_%s.sme_30_%s0.r1r2_11.P5P5.out\n", beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), n_conf[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith].c_str(), sme_[ismear].c_str());
    printf("READING: DATA_SMEAR/%s/%s/2pts.%s.%s.m1m2_%s%s.th_%s.sme_30_%s0.r1r2_00.P5P5.out\n", beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), n_conf[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[6-ith].c_str(), sme_[ismear].c_str());
    printf("READING: DATA_SMEAR/%s/%s/2pts.%s.%s.m1m2_%s%s.th_%s.sme_30_%s0.r1r2_11.P5P5.out\n", beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), n_conf[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[6-ith].c_str(), sme_[ismear].c_str());

    if ((fr_2pts_th_r1r2_00 = fopen(open_2pts_th_r1r2_00, "r")) == NULL ){
      printf("Error opening the file to read: 2pts_0m1_th1_r1r2_00\n");
      exit(EXIT_FAILURE);
    }
    if ((fr_2pts_th_r1r2_11 = fopen(open_2pts_th_r1r2_11, "r")) == NULL ){
      printf("Error opening the file to  read: %s\n",open_2pts_th_r1r2_11);
      exit(EXIT_FAILURE);
    }
    if ((fr_2pts_th_rev_r1r2_00 = fopen(open_2pts_th_rev_r1r2_00, "r")) == NULL ){
      printf("Error opening the file to read: 2pts_0m1_th1_rev_r1r2_00 %s\n",open_2pts_th_rev_r1r2_00);
      exit(EXIT_FAILURE);
    }
    if ((fr_2pts_th_rev_r1r2_11 = fopen(open_2pts_th_rev_r1r2_11, "r")) == NULL ){
      printf("Error opening the file to read: 2pts_0m1_th1_rev_r1r2_11 %s\n",open_2pts_th_rev_r1r2_11);
      exit(EXIT_FAILURE);
    }
    free(open_2pts_th_rev_r1r2_00);
    free(open_2pts_th_rev_r1r2_11);

    clean_string_cluster(fr_2pts_th_r1r2_00, temp_2pts_th_r1r2_00, Nt, clusterfile);
    clean_string_cluster(fr_2pts_th_r1r2_11, temp_2pts_th_r1r2_11, Nt, clusterfile);
    clean_string_cluster(fr_2pts_th_rev_r1r2_00, temp_2pts_th_rev_r1r2_00, Nt, clusterfile);
    clean_string_cluster(fr_2pts_th_rev_r1r2_11, temp_2pts_th_rev_r1r2_11, Nt, clusterfile);
    
    fclose(fr_2pts_th_r1r2_00);
    fclose(fr_2pts_th_r1r2_11);
    fclose(fr_2pts_th_rev_r1r2_00);
    fclose(fr_2pts_th_rev_r1r2_11);
    
    th_r_average(temp_2pts_th_r1r2_00, temp_2pts_th_r1r2_11, temp_2pts_th_rev_r1r2_00, temp_2pts_th_rev_r1r2_11, temp_corr_2pts, clusterfile*Nt);
    
    simmetrize(temp_corr_2pts, corr_2pts, Nt, clusterfile);
    
    free(temp_2pts_th_rev_r1r2_00);
    free(temp_2pts_th_rev_r1r2_11);
    
  }// ith != 3

  if(ith == 3){

    // CASO ith1 == 3
    sprintf(open_2pts_th_r1r2_00, "../Correlator_2pts/DATA_SMEAR/%s/%s/2pts.%s.%s.m1m2_%s%s.th_%s.sme_30_%s0.r1r2_00.P5P5.out", beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), n_conf[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith].c_str(), sme_[ismear].c_str());
    sprintf(open_2pts_th_r1r2_11, "../Correlator_2pts/DATA_SMEAR/%s/%s/2pts.%s.%s.m1m2_%s%s.th_%s.sme_30_%s0.r1r2_11.P5P5.out", beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), n_conf[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith].c_str(), sme_[ismear].c_str());
    
    if ((fr_2pts_th_r1r2_00 = fopen(open_2pts_th_r1r2_00, "r")) == NULL ){
      printf("Error opening the file to read: 2pts_0m1_th1_r1r2_00 %s\n",open_2pts_th_r1r2_00);
      exit(EXIT_FAILURE);
    }
    if ((fr_2pts_th_r1r2_11 = fopen(open_2pts_th_r1r2_11, "r")) == NULL ){
      printf("Error opening the file to read: 2pts_0m1_th1_r1r2_11 %s\n",open_2pts_th_r1r2_11);
      exit(EXIT_FAILURE);
    }
    
    clean_string_cluster(fr_2pts_th_r1r2_00, temp_2pts_th_r1r2_00, Nt, clusterfile);
    clean_string_cluster(fr_2pts_th_r1r2_11, temp_2pts_th_r1r2_11, Nt, clusterfile);
    
    fclose(fr_2pts_th_r1r2_00);
    fclose(fr_2pts_th_r1r2_11);
    
    r_average(temp_2pts_th_r1r2_00, temp_2pts_th_r1r2_11, temp_corr_2pts, clusterfile*Nt);
    
    simmetrize(temp_corr_2pts, corr_2pts, Nt, clusterfile);
    
  }// ith == 3
  
  free(open_2pts_th_r1r2_00);
  free(open_2pts_th_r1r2_11);
  
  free(temp_2pts_th_r1r2_00);
  free(temp_2pts_th_r1r2_11);
  free(temp_corr_2pts);

} // read_corr_2pts


void read_2pts_energy( double *energy_2pts, int clusterfile, const string *beta_V, const string mu_sea_1[Nbeta][Nmusea], const string *mu_sea_2, const string *m_, const string *th_, const string *sme_, int ibeta, int imusea, int im1, int im2, int ith, int ismear, int iE){

  int LEN_NAME = 1024;
  
  FILE *fr_E_from_2pts;
  char *open_E_from_2pts=(char*)malloc(sizeof(char)*LEN_NAME);
  
  if(ith != 3){

    //// ENERGY FROM 2pts sinhDR
    if(iE == 0){
      sprintf(open_E_from_2pts, "../Correlator_2pts/OUTPUT_SMEAR/M/fit_data/%s/Fit_%s/fit_2pts_E_from_sinh_DR.%s.m1m2_%s%s.th_%s.sme_30_%s0.r1r2_Av_00_11.P5P5.out", beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith].c_str(), sme_[ismear].c_str());
    }// iE == 0
    
    //// ENERGY FROM 2pts stdDR
    if(iE == 1){
      sprintf(open_E_from_2pts, "../Correlator_2pts/OUTPUT_SMEAR/M/fit_data/%s/Fit_%s/fit_2pts_E_from_std_DR.%s.m1m2_%s%s.th_%s.sme_30_%s0.r1r2_Av_00_11.P5P5.out", beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith].c_str(), sme_[ismear].c_str());
    }// iE == 1
    
    //// ENERGY FROM 2pts FIT
    if(iE == 2){     
      sprintf(open_E_from_2pts, "../Correlator_2pts/OUTPUT_SMEAR/M/fit_data/%s/Fit_%s/fit_2pts_E.%s.m1m2_%s%s.th_%s.sme_30_%s0.r1r2_Av_00_11.P5P5.out", beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith].c_str(), sme_[ismear].c_str()); 
    }// iE == 2

  }// ith != 3
  
  if(ith == 3){

    // MASS FROM 2pts FIT    
    sprintf(open_E_from_2pts, "../Correlator_2pts/OUTPUT_SMEAR/M/fit_data/%s/Fit_%s/fit_2pts_M.%s.m1m2_%s%s.th_%s.sme_30_%s0.r1r2_Av_00_11.P5P5.out", beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith].c_str(), sme_[ismear].c_str());
    
  }// ith == 3    


  if ((fr_E_from_2pts = fopen(open_E_from_2pts, "r")) == NULL ){
    printf("Error opening the file to read: %s\n", open_E_from_2pts);
    exit(EXIT_FAILURE);
  }
  
  for(int ijk = 0; ijk < clusterfile; ijk++){
    fscanf(fr_E_from_2pts,"%lf", &energy_2pts[ijk]);
  }
  fclose(fr_E_from_2pts);
  
  free(open_E_from_2pts);
  
}// read_2pts_energy


void read_2pts_zeta( double *zeta_2pts, int clusterfile, const string *beta_V, const string mu_sea_1[Nbeta][Nmusea], const string *mu_sea_2, const string *m_, const string *th_, const string *sme_, int ibeta, int imusea, int im1, int im2, int ith, int ismear){

  int LEN_NAME = 1024;
  
  FILE *fr_Z_from_2pts;
  char *open_Z_from_2pts=(char*)malloc(sizeof(char)*LEN_NAME);
  
  sprintf(open_Z_from_2pts, "../Correlator_2pts/OUTPUT_SMEAR/M/fit_data/%s/Fit_%s/fit_2pts_Z.%s.m1m2_%s%s.th_%s.sme_30_%s0.r1r2_Av_00_11.P5P5.out", beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith].c_str(), sme_[ismear].c_str());
  
  if ((fr_Z_from_2pts = fopen(open_Z_from_2pts, "r")) == NULL ){
    printf("Error opening the file to read: %s\n", open_Z_from_2pts);
    exit(EXIT_FAILURE);
  }
  
  for(int ijk = 0; ijk < clusterfile; ijk++){
    fscanf(fr_Z_from_2pts,"%lf", &zeta_2pts[ijk]);    
  }
  fclose(fr_Z_from_2pts);

  free(open_Z_from_2pts);
  
}// read_2pts_zeta


void plot_plateau( double *fit_operatore, double sigma_fit_operatore, double *operatore, double *sigma_operatore, int *fit_interval, int clusterfile, int Nt, int *Tsep_num, string *strategy, string *dir_S, string *dir_E, const string *beta_V, const string mu_sea_1[Nbeta][Nmusea], string *op, const string *mu_sea_2, const string *m_, const string *th_, int iS, int iE, int ibeta, int imusea, int im1, int im2, int ith1, int ith2, int istrategy, int iop){

  int LEN_NAME = 1024;
  
  FILE *fout_plateaux;
  char *file_out_plateaux=(char*)malloc(sizeof(char)*LEN_NAME);
  
  sprintf(file_out_plateaux, "OUTPUT_SMEAR/%s/%s/Plateaux/%s/%s/%s/%s_plateaux.%s.m1m2_%s%s.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), op[iop].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  
  if ((fout_plateaux = fopen(file_out_plateaux, "w")) == NULL ){
    printf("Error opening the output file: %s\n", file_out_plateaux);
    exit(EXIT_FAILURE);
  }
  free(file_out_plateaux);
  
  fprintf(fout_plateaux,"@type xydy\n");
  for(int t = 0; t <= Tsep_num[ibeta]; t++){
    fprintf(fout_plateaux,"%d %f %f\n", t, operatore[(clusterfile-1)*Nt+t], sigma_operatore[t]);
  }
  
  fprintf(fout_plateaux,"@type xy\n");
  fprintf(fout_plateaux,"%d %f\n", Tsep_num[ibeta]/2 - fit_interval[ibeta], fit_operatore[clusterfile-1] + sigma_fit_operatore);
  fprintf(fout_plateaux,"%d %f\n", Tsep_num[ibeta]/2 + fit_interval[ibeta], fit_operatore[clusterfile-1] + sigma_fit_operatore);
  
  fprintf(fout_plateaux,"&\n");
  fprintf(fout_plateaux,"%d %f\n", Tsep_num[ibeta]/2 - fit_interval[ibeta], fit_operatore[clusterfile-1] - sigma_fit_operatore);
  fprintf(fout_plateaux,"%d %f\n", Tsep_num[ibeta]/2 + fit_interval[ibeta], fit_operatore[clusterfile-1] - sigma_fit_operatore);
  
  fprintf(fout_plateaux,"&\n");
  fprintf(fout_plateaux,"%d %f\n", Tsep_num[ibeta]/2 - fit_interval[ibeta], fit_operatore[clusterfile-1]);
  fprintf(fout_plateaux,"%d %f\n", Tsep_num[ibeta]/2 + fit_interval[ibeta], fit_operatore[clusterfile-1]);
  
  fclose(fout_plateaux);
  
}// plot_plateau


void write_output_matrix_el_fit( double *out_array, string *dir_S, string *dir_E, string *op, string *strategy, const string *beta_V, const string mu_sea_1[Nbeta][Nmusea], const string *mu_sea_2, const string *m_, const string *th_, int clusterfile, int iS, int iE, int ibeta, int imusea, int im1, int im2, int ith1, int ith2, int iop, int istrategy){

  int LEN_NAME = 1024;
  
  FILE *fout;
  char *file_output=(char*)malloc(sizeof(char)*LEN_NAME);

  sprintf(file_output, "OUTPUT_SMEAR/%s/%s/%s/%s/%s/%s/%s_fit.%s.m1m2_%s%s.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), op[iop].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), op[iop].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith1].c_str(), th_[ith2].c_str());

  if ((fout = fopen(file_output, "w")) == NULL ){
    printf("Error opening the output file: %s\n", file_output);
    exit(EXIT_FAILURE);
  }
  free(file_output);
  
  for(int ijk = 0; ijk < clusterfile; ijk++){
    fprintf(fout, "%f\n", out_array[ijk]);
  }
  fclose(fout);

}// write_output_matrix_el_fit


void write_output_form_factors_fit( double *out_array, double *q2, string *dir_S, string *dir_E, string *ffs, string *strategy, const string *beta_V, const string mu_sea_1[Nbeta][Nmusea], const string *mu_sea_2, const string *m_, const string *th_, int clusterfile, int iS, int iE, int ibeta, int imusea, int im1, int im2, int ith1, int ith2, int iffs, int istrategy){

  int LEN_NAME = 1024;
  
  FILE *fout;
  char *file_output=(char*)malloc(sizeof(char)*LEN_NAME);

  sprintf(file_output, "OUTPUT_SMEAR/%s/%s/%s/%s/%s/%s/%s.%s.m1m2_%s%s.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), ffs[iffs].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), ffs[iffs].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith1].c_str(), th_[ith2].c_str());

  if ((fout = fopen(file_output, "w")) == NULL ){
    printf("Error opening the output file: %s\n", file_output);
    exit(EXIT_FAILURE);
  }
  free(file_output);
  
//  printf("%s\n",file_output);

  for(int ijk = 0; ijk < clusterfile; ijk++){
    fprintf(fout, "%f\t%f\n", out_array[ijk], q2[ijk]);
  }
  fclose(fout);

}// write_output_form_factors_fit


void read_ffs_and_take_the_error( double form_factor[Nth1][Nth2][clusterfile], double q2[Nth1][Nth2][clusterfile], double sigma_form_factor[Nth1][Nth2], string *dir_S, string *dir_E, const string *beta_V, const string mu_sea_1[Nbeta][Nmusea], const string *mu_sea_2, const string *m_, const string *th_, int clusterfile, int Nth1, int Nth2, int istrategy, int iffs, int iS, int iE, int ibeta, int imusea, int im1, int im2){

  int LEN_NAME = 1024;
  
  FILE *fr_form_factor;
  char *open_form_factor=(char*)malloc(sizeof(char)*LEN_NAME);
  
  string strategy[5] = {"light_to_Heavy", "Heavy_to_light", "Hl_sym", "Hl_sym_av", "Hl_sym_av_S_ratio"};
  string ffs[6] = {"fzero", "fplus", "fminus", "ftens", "hminus", "hplus"};

  double *temp_form_factor=(double*)malloc(sizeof(double)*clusterfile);

  if( (iffs == 1) || (iffs == 2) || (iffs == 4)){            // Per fzero considero anche (ith1=3 & ith2=3), negli altri casi devo escluderlo. 
                                               // Per ftens nel caso elastico devo escludere ith1=ith2
    for(int ith1 = 0; ith1 < Nth1; ith1++){
      for(int ith2 = 0; ith2 < Nth2; ith2++){

       if(ith1 != 3 || ith2 != 3){

         sprintf(open_form_factor, "OUTPUT_SMEAR/%s/%s/%s/%s/%s/%s/%s.%s.m1m2_%s%s.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), ffs[iffs].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), ffs[iffs].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith1].c_str(), th_[ith2].c_str());

         if ((fr_form_factor = fopen(open_form_factor, "r")) == NULL ){
           printf("Error opening the input file: %s\n", open_form_factor);
           exit(EXIT_FAILURE);
         }

         for(int ijk = 0; ijk < clusterfile; ijk++){

           fscanf(fr_form_factor,"%lf %lf", &form_factor[ith1][ith2][ijk], &q2[ith1][ith2][ijk]);

           temp_form_factor[ijk] = form_factor[ith1][ith2][ijk];

	       }// ijk

         fclose(fr_form_factor);

         sigma_form_factor[ith1][ith2] = sigma_JK(temp_form_factor, clusterfile);

	     }// if(ith1 != 3 || ith2 != 3){

      }// ith2
    }// ith1
    
  }else if( iffs == 3 ){ // Per ftens nel caso elastico devo escludere ith1=ith2

    if(im1 != im2){      // ftens: Caso non elastico

      for(int ith1 = 0; ith1 < Nth1; ith1++){
       for(int ith2 = 0; ith2 < Nth2; ith2++){

         if(ith1 != 3 || ith2 != 3){

           sprintf(open_form_factor, "OUTPUT_SMEAR/%s/%s/%s/%s/%s/%s/%s.%s.m1m2_%s%s.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), ffs[iffs].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), ffs[iffs].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith1].c_str(), th_[ith2].c_str());

           if ((fr_form_factor = fopen(open_form_factor, "r")) == NULL ){
             printf("Error opening the input file: %s\n", open_form_factor);
             exit(EXIT_FAILURE);
           }

           for(int ijk = 0; ijk < clusterfile; ijk++){

             fscanf(fr_form_factor,"%lf %lf", &form_factor[ith1][ith2][ijk], &q2[ith1][ith2][ijk]);

             temp_form_factor[ijk] = form_factor[ith1][ith2][ijk];

	         }// ijk

           fclose(fr_form_factor);

           sigma_form_factor[ith1][ith2] = sigma_JK(temp_form_factor, clusterfile);

	       }// if(ith1 != 3 || ith2 != 3){

	     }// ith2
      }// ith1
      
    }else if(im1 == im2){  // ftens: Caso elastico

      for(int ith1 = 0; ith1 < Nth1; ith1++){
       for(int ith2 = 0; ith2 < Nth2; ith2++){

         if(ith1 != ith2){

           sprintf(open_form_factor, "OUTPUT_SMEAR/%s/%s/%s/%s/%s/%s/%s.%s.m1m2_%s%s.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), ffs[iffs].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), ffs[iffs].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith1].c_str(), th_[ith2].c_str());

           if ((fr_form_factor = fopen(open_form_factor, "r")) == NULL ){
             printf("Error opening the input file: %s\n", open_form_factor);
             exit(EXIT_FAILURE);
           }

           for(int ijk = 0; ijk < clusterfile; ijk++){

             fscanf(fr_form_factor,"%lf %lf", &form_factor[ith1][ith2][ijk], &q2[ith1][ith2][ijk]);

             temp_form_factor[ijk] = form_factor[ith1][ith2][ijk];

	          }// ijk

           fclose(fr_form_factor);

           sigma_form_factor[ith1][ith2] = sigma_JK(temp_form_factor, clusterfile);

	         }// if(ith1 != ith2)

	       }// ith2
      }// ith1

    }

  }else if((iffs == 0) || (iffs == 5)){   // Per fzero considero anche (ith1=3 & ith2=3)

    for(int ith1 = 0; ith1 < Nth1; ith1++){
      for(int ith2 = 0; ith2 < Nth2; ith2++){

       sprintf(open_form_factor, "OUTPUT_SMEAR/%s/%s/%s/%s/%s/%s/%s.%s.m1m2_%s%s.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), ffs[iffs].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), ffs[iffs].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith1].c_str(), th_[ith2].c_str());

       if ((fr_form_factor = fopen(open_form_factor, "r")) == NULL ){
         printf("Error opening the input file: %s\n", open_form_factor);
         exit(EXIT_FAILURE);
       }

       for(int ijk = 0; ijk < clusterfile; ijk++){

         fscanf(fr_form_factor,"%lf %lf", &form_factor[ith1][ith2][ijk], &q2[ith1][ith2][ijk]);

         temp_form_factor[ijk] = form_factor[ith1][ith2][ijk];

	     }// ijk

       fclose(fr_form_factor);

       sigma_form_factor[ith1][ith2] = sigma_JK(temp_form_factor, clusterfile);

      }// ith2
    }// ith1
    
  }// if(iffs)

  free(open_form_factor);	  
  free(temp_form_factor);
  
}// read_ffs_and_take_the_error


void write_plot_grace_ffs( double form_factor[Nth1][Nth2][clusterfile], double q2[Nth1][Nth2][clusterfile], double sigma_form_factor[Nth1][Nth2], string *dir_S, string *dir_E, const string *beta_V_2, const string *mu_sea_2, const string *m_, int clusterfile, int Nth1, int Nth2, int istrategy, int iffs, int iS, int iE, int ibeta, int imusea, int im1, int im2){

  int LEN_NAME = 1024;
  
  FILE *fout_form_factor;
  char *file_out_form_factor=(char*)malloc(sizeof(char)*LEN_NAME);
  
  string strategy[5] = {"light_to_Heavy", "Heavy_to_light", "Hl_sym", "Hl_sym_av", "Hl_sym_av_S_ratio"};
  string ffs[6] = {"fzero", "fplus", "fminus", "ftens", "hminus", "hplus"};

  sprintf(file_out_form_factor, "OUTPUT_SMEAR/%s/%s/%s/%s/%s_vs_q2/%s_vs_q2.%s.%s.m1m2_%s%s.xmg", dir_S[iS].c_str(), dir_E[iE].c_str(), ffs[iffs].c_str(), strategy[istrategy].c_str(), ffs[iffs].c_str(), ffs[iffs].c_str(), beta_V_2[ibeta].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str());
  
  if ((fout_form_factor = fopen(file_out_form_factor, "w")) == NULL ){
    printf("Error opening the output file: %s\n", file_out_form_factor);
    exit(EXIT_FAILURE);
  }
  free(file_out_form_factor);
  
  fprintf(fout_form_factor, "@type xydy\n");
  
  if( (iffs == 1) || (iffs == 2) || (iffs == 4)){            // Per fzero considero anche (ith1=3 & ith2=3), negli altri casi devo escluderlo.
                                              // Per ftens nel caso elastico devo escludere ith1=ith2
    for(int ith1 = 0; ith1 < Nth1; ith1++){
      for(int ith2 = 0; ith2 < Nth2; ith2++){

       if(ith1 != 3 || ith2 != 3){

         fprintf(fout_form_factor, "%f %f %f\n", q2[ith1][ith2][15], form_factor[ith1][ith2][15], sigma_form_factor[ith1][ith2]);

	     }// 	if(ith1 != 3 || ith2 != 3){

      }// ith2

      if(ith1 != 3){
       fprintf(fout_form_factor, "&\n");
      }// if(ith1 != 3)
      
    }// ith1

  }else if( iffs == 3 ){ // Per ftens nel caso elastico devo escludere ith1=ith2

    if(im1 != im2){      // ftens: Caso non elastico

      for(int ith1 = 0; ith1 < Nth1; ith1++){
       for(int ith2 = 0; ith2 < Nth2; ith2++){

         if(ith1 != 3 || ith2 != 3){

           fprintf(fout_form_factor, "%f %f %f\n", q2[ith1][ith2][15], form_factor[ith1][ith2][15], sigma_form_factor[ith1][ith2]);

	  }// 	if(ith1 != 3 || ith2 != 3){

	}// ith2
	
	if(ith1 != 3){
   fprintf(fout_form_factor, "&\n");
	}// if(ith1 != 3)
	
      }// ith1

    }else if(im1 == im2){  // ftens: Caso elastico

      for(int ith1 = 0; ith1 < Nth1; ith1++){
       for(int ith2 = 0; ith2 < Nth2; ith2++){

         if(ith1 != ith2 ){

           fprintf(fout_form_factor, "%f %f %f\n", q2[ith1][ith2][15], form_factor[ith1][ith2][15], sigma_form_factor[ith1][ith2]);

	  }// 	if(ith1 != ith2 )
	  
	}// ith2
	
	if(ith1 != 3){
   fprintf(fout_form_factor, "&\n");
	}// if(ith1 != 3)
	
      }// ith1

    }

  }else if((iffs == 0) || (iffs == 5)){   // Per fzero considero anche (ith1=3 & ith2=3)

    for(int ith1 = 0; ith1 < Nth1; ith1++){
      for(int ith2 = 0; ith2 < Nth2; ith2++){

       fprintf(fout_form_factor, "%f %f %f\n", q2[ith1][ith2][15], form_factor[ith1][ith2][15], sigma_form_factor[ith1][ith2]);

      }// ith2
      
      if(ith1 != 3){
       fprintf(fout_form_factor, "&\n");
      }// if(ith1 != 3)
      
    }// ith1
    
  }// if(iffs)
  
  fclose(fout_form_factor);
  
}// write_plot_grace_ffs

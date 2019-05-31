#include <TF1.h>
#include <fstream>
#include "stat_analysis_func.h"
#include "fit_3pts_func.h"
#include "clean_str_cl.h"
#include "fit_2pts_func.h"
#include "definitions.h"

using namespace std;

#define LEN_NAME 1024
#define PI 3.141592653589793

//////////////////////////////

int f0_with_S = 1;
int f0_without_S = 0;

//////////////////////////////

//////////////////////////////
   
int energia_da_sinhDR = 0;  
int energia_da_stdDR = 1;
int energia_dal_fit = 0;

//////////////////////////////

int ibeta, imusea, im1, im2;

int main(){
  
  //////////////////// LEGGO IL FILE DI INPUT
  
  FILE *fi;
  
  if ((fi = fopen("Input_corr_3pts/file_input_corr_3pts.out", "r")) == NULL ){
    printf("Error opening the input file!!\n");
    exit(EXIT_FAILURE);
  }
  fscanf(fi, "%d %d %d %d", &ibeta, &imusea, &im1, &im2);
  fclose(fi);
  
  //////////////////// FINE LEGGO IL FILE DI INPUT
    
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
  
  /////// LEGGO IL FATT. DI FORMA PER OGNI CINEMATICA, ED IL q2 RELATIVO, E NE CALCOLO L'ERRORE JK  
  double q2[Nth1][Nth2][clusterfile];

  double fzero_ml_to_mH[Nth1][Nth2][clusterfile], fzero_mH_to_ml[Nth1][Nth2][clusterfile], fzero_Hl_sym[Nth1][Nth2][clusterfile], fzero_Hl_sym_av[Nth1][Nth2][clusterfile], fzero_Hl_sym_av_S_ratio[Nth1][Nth2][clusterfile];
  double fplus_ml_to_mH[Nth1][Nth2][clusterfile], fplus_mH_to_ml[Nth1][Nth2][clusterfile], fplus_Hl_sym[Nth1][Nth2][clusterfile], fplus_Hl_sym_av[Nth1][Nth2][clusterfile], fplus_Hl_sym_av_S_ratio[Nth1][Nth2][clusterfile];
  double fminus_ml_to_mH[Nth1][Nth2][clusterfile], fminus_mH_to_ml[Nth1][Nth2][clusterfile], fminus_Hl_sym[Nth1][Nth2][clusterfile], fminus_Hl_sym_av[Nth1][Nth2][clusterfile], fminus_Hl_sym_av_S_ratio[Nth1][Nth2][clusterfile];
  
  double sigma_fzero_ml_to_mH[Nth1][Nth2], sigma_fzero_mH_to_ml[Nth1][Nth2], sigma_fzero_Hl_sym[Nth1][Nth2], sigma_fzero_Hl_sym_av[Nth1][Nth2], sigma_fzero_Hl_sym_av_S_ratio[Nth1][Nth2];
  double sigma_fplus_ml_to_mH[Nth1][Nth2], sigma_fplus_mH_to_ml[Nth1][Nth2], sigma_fplus_Hl_sym[Nth1][Nth2], sigma_fplus_Hl_sym_av[Nth1][Nth2], sigma_fplus_Hl_sym_av_S_ratio[Nth1][Nth2];
  double sigma_fminus_ml_to_mH[Nth1][Nth2], sigma_fminus_mH_to_ml[Nth1][Nth2], sigma_fminus_Hl_sym[Nth1][Nth2], sigma_fminus_Hl_sym_av[Nth1][Nth2], sigma_fminus_Hl_sym_av_S_ratio[Nth1][Nth2];
  
  //// fzero ////
  read_ffs_and_take_the_error( fzero_ml_to_mH, q2, sigma_fzero_ml_to_mH, dir_S, dir_E, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, Nth1, Nth2, 0, 0, iS, iE, ibeta, imusea, im1, im2);
  read_ffs_and_take_the_error( fzero_mH_to_ml, q2, sigma_fzero_mH_to_ml, dir_S, dir_E, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, Nth1, Nth2, 1, 0, iS, iE, ibeta, imusea, im1, im2);
  read_ffs_and_take_the_error( fzero_Hl_sym, q2, sigma_fzero_Hl_sym, dir_S, dir_E, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, Nth1, Nth2, 2, 0, iS, iE, ibeta, imusea, im1, im2);
  read_ffs_and_take_the_error( fzero_Hl_sym_av, q2, sigma_fzero_Hl_sym_av, dir_S, dir_E, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, Nth1, Nth2, 3, 0, iS, iE, ibeta, imusea, im1, im2);
  read_ffs_and_take_the_error( fzero_Hl_sym_av_S_ratio, q2, sigma_fzero_Hl_sym_av_S_ratio, dir_S, dir_E, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, Nth1, Nth2, 4, 0, iS, iE, ibeta, imusea, im1, im2);
  
  //// fplus ////
  read_ffs_and_take_the_error( fplus_ml_to_mH, q2, sigma_fplus_ml_to_mH, dir_S, dir_E, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, Nth1, Nth2, 0, 1, iS, iE, ibeta, imusea, im1, im2);
  read_ffs_and_take_the_error( fplus_mH_to_ml, q2, sigma_fplus_mH_to_ml, dir_S, dir_E, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, Nth1, Nth2, 1, 1, iS, iE, ibeta, imusea, im1, im2);
  read_ffs_and_take_the_error( fplus_Hl_sym, q2, sigma_fplus_Hl_sym, dir_S, dir_E, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, Nth1, Nth2, 2, 1, iS, iE, ibeta, imusea, im1, im2);
  read_ffs_and_take_the_error( fplus_Hl_sym_av, q2, sigma_fplus_Hl_sym_av, dir_S, dir_E, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, Nth1, Nth2, 3, 1, iS, iE, ibeta, imusea, im1, im2);
  read_ffs_and_take_the_error( fplus_Hl_sym_av_S_ratio, q2, sigma_fplus_Hl_sym_av_S_ratio, dir_S, dir_E, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, Nth1, Nth2, 4, 1, iS, iE, ibeta, imusea, im1, im2);


  //// fminus ////
  read_ffs_and_take_the_error( fminus_ml_to_mH, q2, sigma_fminus_ml_to_mH, dir_S, dir_E, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, Nth1, Nth2, 0, 2, iS, iE, ibeta, imusea, im1, im2);
  read_ffs_and_take_the_error( fminus_mH_to_ml, q2, sigma_fminus_mH_to_ml, dir_S, dir_E, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, Nth1, Nth2, 1, 2, iS, iE, ibeta, imusea, im1, im2);
  read_ffs_and_take_the_error( fminus_Hl_sym, q2, sigma_fminus_Hl_sym, dir_S, dir_E, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, Nth1, Nth2, 2, 2, iS, iE, ibeta, imusea, im1, im2);
  read_ffs_and_take_the_error( fminus_Hl_sym_av, q2, sigma_fminus_Hl_sym_av, dir_S, dir_E, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, Nth1, Nth2, 3, 2, iS, iE, ibeta, imusea, im1, im2);
  read_ffs_and_take_the_error( fminus_Hl_sym_av_S_ratio, q2, sigma_fminus_Hl_sym_av_S_ratio, dir_S, dir_E, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, Nth1, Nth2, 4, 2, iS, iE, ibeta, imusea, im1, im2);

  /////// FINE LEGGO IL FATT. DI FORMA PER OGNI CINEMATICA E NE CALCOLO L'ERRORE JK

  
  /////// SCRIVO L'OUTPUT GRACE CON IL PLOT DEL FATTORE DI FORMA

  //// fzero ////
  write_plot_grace_ffs( fzero_ml_to_mH, q2, sigma_fzero_ml_to_mH, dir_S, dir_E, beta_V_2, mu_sea_2, m_, clusterfile, Nth1, Nth2, 0, 0, iS, iE, ibeta, imusea, im1, im2);
  write_plot_grace_ffs( fzero_mH_to_ml, q2, sigma_fzero_mH_to_ml, dir_S, dir_E, beta_V_2, mu_sea_2, m_, clusterfile, Nth1, Nth2, 1, 0, iS, iE, ibeta, imusea, im1, im2);
  write_plot_grace_ffs( fzero_Hl_sym, q2, sigma_fzero_Hl_sym, dir_S, dir_E, beta_V_2, mu_sea_2, m_, clusterfile, Nth1, Nth2, 2, 0, iS, iE, ibeta, imusea, im1, im2);
  write_plot_grace_ffs( fzero_Hl_sym_av, q2, sigma_fzero_Hl_sym_av, dir_S, dir_E, beta_V_2, mu_sea_2, m_, clusterfile, Nth1, Nth2, 3, 0, iS, iE, ibeta, imusea, im1, im2);
  write_plot_grace_ffs( fzero_Hl_sym_av_S_ratio, q2, sigma_fzero_Hl_sym_av_S_ratio, dir_S, dir_E, beta_V_2, mu_sea_2, m_, clusterfile, Nth1, Nth2, 4, 0, iS, iE, ibeta, imusea, im1, im2);

  //// fplus ////
  write_plot_grace_ffs( fplus_ml_to_mH, q2, sigma_fplus_ml_to_mH, dir_S, dir_E, beta_V_2, mu_sea_2, m_, clusterfile, Nth1, Nth2, 0, 1, iS, iE, ibeta, imusea, im1, im2);
  write_plot_grace_ffs( fplus_mH_to_ml, q2, sigma_fplus_mH_to_ml, dir_S, dir_E, beta_V_2, mu_sea_2, m_, clusterfile, Nth1, Nth2, 1, 1, iS, iE, ibeta, imusea, im1, im2);
  write_plot_grace_ffs( fplus_Hl_sym, q2, sigma_fplus_Hl_sym, dir_S, dir_E, beta_V_2, mu_sea_2, m_, clusterfile, Nth1, Nth2, 2, 1, iS, iE, ibeta, imusea, im1, im2);
  write_plot_grace_ffs( fplus_Hl_sym_av, q2, sigma_fplus_Hl_sym_av, dir_S, dir_E, beta_V_2, mu_sea_2, m_, clusterfile, Nth1, Nth2, 3, 1, iS, iE, ibeta, imusea, im1, im2);
  write_plot_grace_ffs( fplus_Hl_sym_av_S_ratio, q2, sigma_fplus_Hl_sym_av_S_ratio, dir_S, dir_E, beta_V_2, mu_sea_2, m_, clusterfile, Nth1, Nth2, 4, 1, iS, iE, ibeta, imusea, im1, im2);
  
  //// fminus ////
  write_plot_grace_ffs( fminus_ml_to_mH, q2, sigma_fminus_ml_to_mH, dir_S, dir_E, beta_V_2, mu_sea_2, m_, clusterfile, Nth1, Nth2, 0, 2, iS, iE, ibeta, imusea, im1, im2);
  write_plot_grace_ffs( fminus_mH_to_ml, q2, sigma_fminus_mH_to_ml, dir_S, dir_E, beta_V_2, mu_sea_2, m_, clusterfile, Nth1, Nth2, 1, 2, iS, iE, ibeta, imusea, im1, im2);
  write_plot_grace_ffs( fminus_Hl_sym, q2, sigma_fminus_Hl_sym, dir_S, dir_E, beta_V_2, mu_sea_2, m_, clusterfile, Nth1, Nth2, 2, 2, iS, iE, ibeta, imusea, im1, im2);
  write_plot_grace_ffs( fminus_Hl_sym_av, q2, sigma_fminus_Hl_sym_av, dir_S, dir_E, beta_V_2, mu_sea_2, m_, clusterfile, Nth1, Nth2, 3, 2, iS, iE, ibeta, imusea, im1, im2);
  write_plot_grace_ffs( fminus_Hl_sym_av_S_ratio, q2, sigma_fminus_Hl_sym_av_S_ratio, dir_S, dir_E, beta_V_2, mu_sea_2, m_, clusterfile, Nth1, Nth2, 4, 2, iS, iE, ibeta, imusea, im1, im2);
  
  /////// FINE SCRIVO L'OUTPUT GRACE CON IL PLOT DEL FATTORE DI FORMA
  
  return 0;
  
}


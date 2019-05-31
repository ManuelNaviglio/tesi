#include <TMinuit.h>
#include <TF1.h>
#include <fstream>
#include "fit_2pts_func.h"
#include "clean_str_cl.h"
#include "stat_analysis_func.h"
#include "fit_3pts_func.h"
#include "const_fit.h"
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

//////////////////////////////

int f0_scalar_no_ratio = 1;
int f0_scalar_ratio = 0;

//////////////////////////////

string op[3] = {"V0", "Vi", "S0"};
string reim[2] = {"RE", "IM"};
string av[2] = {"SPATIAL_THETA_AV", "TOTAL_AVERAGE"};
string strategy[4] = {"light_to_Heavy", "Heavy_to_light", "Hl_sym", "Hl_sym_av"};
string form_factors[3] = {"fzero", "fplus", "fminus"};

int ibeta, imusea, im1, im2, ith1, ith2, ismear, rev_ith1, rev_ith2;

const int npar = 2;

int T_max[5] = { 64, 48, 64, 48, 96};
int Tsep_num[5] = {18, 18, 20, 20, 26};
int fit_interval[5] = {2, 2, 2, 2, 2};

/////////////// VARIABILI GLOBALI PER IL FIT
double V0_fit_glb, Vi_fit_glb, f0_S_fit_glb;
double sigma_V0_fit_glb, sigma_Vi_fit_glb, sigma_f0_S_fit_glb;
double q2_glb, E_H_glb, E_l_glb, pl_glb, pH_glb, M_l_glb, M_H_glb;
//////////////////////////////////////////////////////////////////

double chi2_f1, chi2_f2, chi2_f3, chi2_f;

void chi2( int &, double *, double &, double *, int);
double V0_function( double *);
double Vi_function( double *);

void sigma_matrix_el( double *, double *, int, int, int);
  
int main(){

  //////////////////// LEGGO IL FILE DI INPUT
  
  FILE *fi;
  
  if ((fi = fopen("Input_corr_3pts/file_input_corr_3pts.out", "r")) == NULL ){
    printf("Error opening the input file!!\n");
    exit(EXIT_FAILURE);
  }
  fscanf(fi, "%d %d %d %d %d %d %d", &ibeta, &imusea, &im1, &im2, &ith1, &ith2, &ismear);
  fclose(fi);

  int Nt = T_max[ibeta];
  
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
  
  string dir_f0_Hl_sym_av[2] = {"Hl_sym_av", "Hl_sym_av_S_ratio"};
  int if0;
  
  if(f0_scalar_no_ratio == 1){
    if0 = 0;
  }
  if(f0_scalar_ratio == 1){
    if0 = 1;
  }
  
  //////////////////////////////
  
  //////////////////////////////
  
  string dir_S[2] = {"with_S", "without_S"};
  int iS, sigma_S_enhancement;
  
  if(f0_with_S == 1){
    sigma_S_enhancement = 1;
    iS = 0;
  }
  if(f0_without_S == 1){
    sigma_S_enhancement = 100000;
    iS = 1;
  }
  
  //////////////////////////////
  
  //////////////////// COSTRUISCO I MOMENTI
  double parity=1;
  double th1_val, th2_val, momentum_1, momentum_2;
  
  if(ith2 <= 3){
    rev_ith1 = ith1;
    rev_ith2 = ith2;
    parity = +1;
  }else{
    rev_ith1 = (Nth-1)-ith1;
    rev_ith2 = (Nth-1)-ith2;
    parity = -1;
  }
  
  th1_val = theta_value[ibeta][ith1];
  th2_val = theta_value[ibeta][ith2];
  
  momentum_1 = (th1_val*PI)/L[ibeta];
  momentum_2 = (th2_val*PI)/L[ibeta];
  
  //////////////////// COSTRUISCO I MOMENTI

  /////////////// COSTRUISCO L'ARRAY mq
  
  double mq[Nbeta][Nmasses] = {{0.0030, 0.01800, 0.02200, 0.02600, 0.21256, 0.25000, 0.29404, 0.34583, 0.40675, 0.47840, 0.56267, 0.66178, 0.77836, 0.91546, 1.07672},
	                             {0.0040, 0.01800, 0.02200, 0.02600, 0.21256, 0.25000, 0.29404, 0.34583, 0.40675, 0.47840, 0.56267, 0.66178, 0.77836, 0.91546, 1.07672},
			                         {0.0025, 0.01550, 0.01900, 0.02250, 0.18705, 0.22000, 0.25875, 0.30433, 0.35794, 0.42099, 0.49515, 0.58237, 0.68495, 0.80561, 0.94752},
			                         {0.0085, 0.01550, 0.01900, 0.02250, 0.18705, 0.22000, 0.25875, 0.30433, 0.35794, 0.42099, 0.49515, 0.58237, 0.68495, 0.80561, 0.94752},
			                         {0.0015, 0.01230, 0.01500, 0.01770, 0.14454, 0.17000, 0.19995, 0.23517, 0.27659, 0.32531, 0.38262, 0.45001, 0.52928, 0.62252, 0.73217}};
    
  mq[ibeta][0] = mq_l[ibeta][imusea]; // questo è il motivo per cui non posso inserirlo in "definitions.h"
  
  /////////////// COSTRUISCO L'ARRAY mq
  
  
  /////////////////////////////////////
  //                                 //
  // READ 3pts CORRELATION FUNCTIONS //
  //                                 //
  /////////////////////////////////////

  ////////////// DICHIARAZIONE DEI CORRELATORI
  double *corr_3pts_ml_to_mH_V0=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  double *corr_3pts_mH_to_ml_V0=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  double *corr_3pts_ml_to_ml_V0=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  double *corr_3pts_mH_to_mH_V0=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));

  double *corr_3pts_ml_to_mH_V0_av=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  double *corr_3pts_mH_to_ml_V0_av=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  double *corr_3pts_ml_to_ml_V0_av=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  double *corr_3pts_mH_to_mH_V0_av=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));

  double *corr_3pts_ml_to_mH_Vi=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  double *corr_3pts_mH_to_ml_Vi=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  double *corr_3pts_ml_to_ml_Vi=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  double *corr_3pts_mH_to_mH_Vi=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));

  double *corr_3pts_ml_to_mH_Vi_av=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  double *corr_3pts_mH_to_ml_Vi_av=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  double *corr_3pts_ml_to_ml_Vi_av=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  double *corr_3pts_mH_to_mH_Vi_av=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  
  double *corr_3pts_ml_to_mH_S0=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  double *corr_3pts_mH_to_ml_S0=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));

  double *corr_3pts_ml_to_mH_S0_av=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  double *corr_3pts_mH_to_ml_S0_av=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));

  double *corr_3pts_rest_ml_to_mH_V0_av=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  double *corr_3pts_rest_mH_to_ml_V0_av=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  double *corr_3pts_rest_ml_to_ml_V0_av=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  double *corr_3pts_rest_mH_to_mH_V0_av=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  
  double *corr_3pts_rest_ml_to_mH_S0_av=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  double *corr_3pts_rest_mH_to_ml_S0_av=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  ////////////// FINE DICHIARAZIONE DEI CORRELATORI
  
  //// SPATIAL_THETA_AV Vi
  read_corr_3pts( corr_3pts_ml_to_mH_Vi, T_max[ibeta], clusterfile, parity, av, beta_V, mu_sea_1, reim, n_conf, mu_sea_2, m_, th_, sme_, op, ibeta, imusea, im1, im2, ith1, ith2, ismear, 0, 1, 1);
  read_corr_3pts( corr_3pts_mH_to_ml_Vi, T_max[ibeta], clusterfile, parity, av, beta_V, mu_sea_1, reim, n_conf, mu_sea_2, m_, th_, sme_, op, ibeta, imusea, im2, im1, rev_ith2, rev_ith1, ismear, 0, 1, 1);
  read_corr_3pts( corr_3pts_ml_to_ml_Vi, T_max[ibeta], clusterfile, parity, av, beta_V, mu_sea_1, reim, n_conf, mu_sea_2, m_, th_, sme_, op, ibeta, imusea, im1, im1, ith1, ith1, ismear, 0, 1, 1);
  read_corr_3pts( corr_3pts_mH_to_mH_Vi, T_max[ibeta], clusterfile, parity, av, beta_V, mu_sea_1, reim, n_conf, mu_sea_2, m_, th_, sme_, op, ibeta, imusea, im2, im2, rev_ith2, rev_ith2, ismear, 0, 1, 1);

  //// TOTAL_AVERAGE Vi
  read_corr_3pts( corr_3pts_ml_to_mH_Vi_av, T_max[ibeta], clusterfile, parity, av, beta_V, mu_sea_1, reim, n_conf, mu_sea_2, m_, th_, sme_, op, ibeta, imusea, im1, im2, ith1, ith2, ismear, 1, 1, 1);
  read_corr_3pts( corr_3pts_mH_to_ml_Vi_av, T_max[ibeta], clusterfile, parity, av, beta_V, mu_sea_1, reim, n_conf, mu_sea_2, m_, th_, sme_, op, ibeta, imusea, im2, im1, rev_ith2, rev_ith1, ismear, 1, 1, 1);
  read_corr_3pts( corr_3pts_ml_to_ml_Vi_av, T_max[ibeta], clusterfile, parity, av, beta_V, mu_sea_1, reim, n_conf, mu_sea_2, m_, th_, sme_, op, ibeta, imusea, im1, im1, ith1, ith1, ismear, 1, 1, 1);
  read_corr_3pts( corr_3pts_mH_to_mH_Vi_av, T_max[ibeta], clusterfile, parity, av, beta_V, mu_sea_1, reim, n_conf, mu_sea_2, m_, th_, sme_, op, ibeta, imusea, im2, im2, rev_ith2, rev_ith2, ismear, 1, 1, 1);

  //// SPATIAL_THETA_AV V0
  read_corr_3pts( corr_3pts_ml_to_mH_V0, T_max[ibeta], clusterfile, parity, av, beta_V, mu_sea_1, reim, n_conf, mu_sea_2, m_, th_, sme_, op, ibeta, imusea, im1, im2, ith1, ith2, ismear, 0, 0, 0);
  read_corr_3pts( corr_3pts_mH_to_ml_V0, T_max[ibeta], clusterfile, parity, av, beta_V, mu_sea_1, reim, n_conf, mu_sea_2, m_, th_, sme_, op, ibeta, imusea, im2, im1, rev_ith2, rev_ith1, ismear, 0, 0, 0);
  read_corr_3pts( corr_3pts_ml_to_ml_V0, T_max[ibeta], clusterfile, parity, av, beta_V, mu_sea_1, reim, n_conf, mu_sea_2, m_, th_, sme_, op, ibeta, imusea, im1, im1, ith1, ith1, ismear, 0, 0, 0);
  read_corr_3pts( corr_3pts_mH_to_mH_V0, T_max[ibeta], clusterfile, parity, av, beta_V, mu_sea_1, reim, n_conf, mu_sea_2, m_, th_, sme_, op, ibeta, imusea, im2, im2, rev_ith2, rev_ith2, ismear, 0, 0, 0);

  //// TOTAL_AVERAGE V0
  read_corr_3pts( corr_3pts_ml_to_mH_V0_av, T_max[ibeta], clusterfile, parity, av, beta_V, mu_sea_1, reim, n_conf, mu_sea_2, m_, th_, sme_, op, ibeta, imusea, im1, im2, ith1, ith2, ismear, 1, 0, 0);
  read_corr_3pts( corr_3pts_mH_to_ml_V0_av, T_max[ibeta], clusterfile, parity, av, beta_V, mu_sea_1, reim, n_conf, mu_sea_2, m_, th_, sme_, op, ibeta, imusea, im2, im1, rev_ith2, rev_ith1, ismear, 1, 0, 0);
  read_corr_3pts( corr_3pts_ml_to_ml_V0_av, T_max[ibeta], clusterfile, parity, av, beta_V, mu_sea_1, reim, n_conf, mu_sea_2, m_, th_, sme_, op, ibeta, imusea, im1, im1, ith1, ith1, ismear, 1, 0, 0);
  read_corr_3pts( corr_3pts_mH_to_mH_V0_av, T_max[ibeta], clusterfile, parity, av, beta_V, mu_sea_1, reim, n_conf, mu_sea_2, m_, th_, sme_, op, ibeta, imusea, im2, im2, rev_ith2, rev_ith2, ismear, 1, 0, 0);

  //// TOTAL_AVERAGE V0 AT REST
  read_corr_3pts( corr_3pts_rest_ml_to_mH_V0_av, T_max[ibeta], clusterfile, parity, av, beta_V, mu_sea_1, reim, n_conf, mu_sea_2, m_, th_, sme_, op, ibeta, imusea, im1, im2, 3, 3, ismear, 1, 0, 0);
  read_corr_3pts( corr_3pts_rest_mH_to_ml_V0_av, T_max[ibeta], clusterfile, parity, av, beta_V, mu_sea_1, reim, n_conf, mu_sea_2, m_, th_, sme_, op, ibeta, imusea, im2, im1, 3, 3, ismear, 1, 0, 0);
  read_corr_3pts( corr_3pts_rest_ml_to_ml_V0_av, T_max[ibeta], clusterfile, parity, av, beta_V, mu_sea_1, reim, n_conf, mu_sea_2, m_, th_, sme_, op, ibeta, imusea, im1, im1, 3, 3, ismear, 1, 0, 0);
  read_corr_3pts( corr_3pts_rest_mH_to_mH_V0_av, T_max[ibeta], clusterfile, parity, av, beta_V, mu_sea_1, reim, n_conf, mu_sea_2, m_, th_, sme_, op, ibeta, imusea, im2, im2, 3, 3, ismear, 1, 0, 0);
  
  //// SPATIAL_THETA_AV S0
  read_corr_3pts( corr_3pts_ml_to_mH_S0, T_max[ibeta], clusterfile, parity, av, beta_V, mu_sea_1, reim, n_conf, mu_sea_2, m_, th_, sme_, op, ibeta, imusea, im1, im2, ith1, ith2, ismear, 0, 2, 0);
  read_corr_3pts( corr_3pts_mH_to_ml_S0, T_max[ibeta], clusterfile, parity, av, beta_V, mu_sea_1, reim, n_conf, mu_sea_2, m_, th_, sme_, op, ibeta, imusea, im2, im1, rev_ith2, rev_ith1, ismear, 0, 2, 0);

  //// TOTAL_AVERAGE S0
  read_corr_3pts( corr_3pts_ml_to_mH_S0_av, T_max[ibeta], clusterfile, parity, av, beta_V, mu_sea_1, reim, n_conf, mu_sea_2, m_, th_, sme_, op, ibeta, imusea, im1, im2, ith1, ith2, ismear, 1, 2, 0);
  read_corr_3pts( corr_3pts_mH_to_ml_S0_av, T_max[ibeta], clusterfile, parity, av, beta_V, mu_sea_1, reim, n_conf, mu_sea_2, m_, th_, sme_, op, ibeta, imusea, im2, im1, rev_ith2, rev_ith1, ismear, 1, 2, 0);

  //// TOTAL_AVERAGE S0 AT REST
  read_corr_3pts( corr_3pts_rest_ml_to_mH_S0_av, T_max[ibeta], clusterfile, parity, av, beta_V, mu_sea_1, reim, n_conf, mu_sea_2, m_, th_, sme_, op, ibeta, imusea, im1, im2, 3, 3, ismear, 1, 2, 0);
  read_corr_3pts( corr_3pts_rest_mH_to_ml_S0_av, T_max[ibeta], clusterfile, parity, av, beta_V, mu_sea_1, reim, n_conf, mu_sea_2, m_, th_, sme_, op, ibeta, imusea, im2, im1, 3, 3, ismear, 1, 2, 0);

  ///////// FINE READ 3pts CORRELATION FUNCTIONS
  

  ////////////////////////////////////////////////////////////////////////////////
  //                                                                            //
  // READ 2pts CORRELATION FUNCTIONS  &  2pts ENERGIES FROM DISPERSION RELATION //
  //                                                                            //
  ////////////////////////////////////////////////////////////////////////////////

  /////////////////
  //   0m1_th1   //
  /////////////////
  
  ////// CORRELATORE A 2 PUNTI A RIPOSO
  double *corr_rest_2pts_0ml=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  read_corr_2pts( corr_rest_2pts_0ml, Nt, clusterfile, beta_V, mu_sea_1, n_conf, mu_sea_2, m_, th_, sme_, ibeta, imusea, 0, im1, 3, ismear);

  // CORRELATORE A 2 PUNTI
  double *corr_2pts_0m1_th1=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  read_corr_2pts( corr_2pts_0m1_th1, Nt, clusterfile, beta_V, mu_sea_1, n_conf, mu_sea_2, m_, th_, sme_, ibeta, imusea, 0, im1, ith1, ismear);
  
  // MASS light
  double *mass_2pts_ml=(double*)malloc(sizeof(double)*clusterfile);
  read_2pts_energy( mass_2pts_ml, clusterfile, beta_V, mu_sea_1, mu_sea_2, m_, th_, sme_, ibeta, imusea, 0, im1, 3, ismear, iE);

  // ZETA FROM 2pts FIT
  double *zeta_2pts_m1_from_fit=(double*)malloc(sizeof(double)*clusterfile);
  read_2pts_zeta( zeta_2pts_m1_from_fit, clusterfile, beta_V, mu_sea_1, mu_sea_2, m_, th_, sme_, ibeta, imusea, 0, im1, ith1, ismear);

  // ENERGY
  double *energy_2pts_m1=(double*)malloc(sizeof(double)*clusterfile);
  read_2pts_energy( energy_2pts_m1, clusterfile, beta_V, mu_sea_1, mu_sea_2, m_, th_, sme_, ibeta, imusea, 0, im1, ith1, ismear, iE);
  
  ///////////// FINE 0m1_th1
  
  
  /////////////////
  //   0m2_th2   //  
  /////////////////
  
  ////// CORRELATORE A 2 PUNTI A RIPOSO
  double *corr_rest_2pts_0mH=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  read_corr_2pts( corr_rest_2pts_0mH, Nt, clusterfile, beta_V, mu_sea_1, n_conf, mu_sea_2, m_, th_, sme_, ibeta, imusea, 0, im2, 3, ismear);

  // CORRELATORE A 2 PUNTI
  double *corr_2pts_0m2_th2=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  read_corr_2pts( corr_2pts_0m2_th2, Nt, clusterfile, beta_V, mu_sea_1, n_conf, mu_sea_2, m_, th_, sme_, ibeta, imusea, 0, im2, rev_ith2, ismear);

  // MASS heavy
  double *mass_2pts_mH=(double*)malloc(sizeof(double)*clusterfile);
  read_2pts_energy( mass_2pts_mH, clusterfile, beta_V, mu_sea_1, mu_sea_2, m_, th_, sme_, ibeta, imusea, 0, im2, 3, ismear, iE);

  // ZETA FROM 2pts FIT
  double *zeta_2pts_m2_from_fit=(double*)malloc(sizeof(double)*clusterfile);
  read_2pts_zeta( zeta_2pts_m2_from_fit, clusterfile, beta_V, mu_sea_1, mu_sea_2, m_, th_, sme_, ibeta, imusea, 0, im2, rev_ith2, ismear);
  
  // ENERGY
  double *energy_2pts_m2=(double*)malloc(sizeof(double)*clusterfile);
  read_2pts_energy( energy_2pts_m2, clusterfile, beta_V, mu_sea_1, mu_sea_2, m_, th_, sme_, ibeta, imusea, 0, im2, rev_ith2, ismear, iE);
  
  ///////////// FINE 0m2_th2
  
  ////////////// FINE READ 2pts CORRELATION FUNCTIONS  &  2pts ENERGIES, ZETA, AND MASSES

  
  //////////////////////////////
  //                          //
  // 2pts correlators_at_Tsep //
  //                          //
  //////////////////////////////
  
  double *tilde_corr_2pts_Tsep_0m1_th1=(double*)malloc(sizeof(double)*clusterfile), *tilde_corr_2pts_Tsep_0m2_th2=(double*)malloc(sizeof(double)*clusterfile);
  double *tilde_corr_rest_2pts_Tsep_0ml=(double*)malloc(sizeof(double)*clusterfile), *tilde_corr_rest_2pts_Tsep_0mH=(double*)malloc(sizeof(double)*clusterfile);
  
  for(int ijk = 0; ijk < clusterfile; ijk++){
    tilde_corr_2pts_Tsep_0m1_th1[ijk]  = 0.5*( corr_2pts_0m1_th1[ijk*Nt+Tsep_num[ibeta]]  + sqrt( pow(corr_2pts_0m1_th1[ijk*Nt+Tsep_num[ibeta]],2)  - pow(corr_2pts_0m1_th1[ijk*Nt+(Nt/2)],2)  ) );
    tilde_corr_2pts_Tsep_0m2_th2[ijk]  = 0.5*( corr_2pts_0m2_th2[ijk*Nt+Tsep_num[ibeta]]  + sqrt( pow(corr_2pts_0m2_th2[ijk*Nt+Tsep_num[ibeta]],2)  - pow(corr_2pts_0m2_th2[ijk*Nt+(Nt/2)],2)  ) );
    tilde_corr_rest_2pts_Tsep_0ml[ijk] = 0.5*( corr_rest_2pts_0ml[ijk*Nt+Tsep_num[ibeta]] + sqrt( pow(corr_rest_2pts_0ml[ijk*Nt+Tsep_num[ibeta]],2) - pow(corr_rest_2pts_0ml[ijk*Nt+(Nt/2)],2) ) );
    tilde_corr_rest_2pts_Tsep_0mH[ijk] = 0.5*( corr_rest_2pts_0mH[ijk*Nt+Tsep_num[ibeta]] + sqrt( pow(corr_rest_2pts_0mH[ijk*Nt+Tsep_num[ibeta]],2) - pow(corr_rest_2pts_0mH[ijk*Nt+(Nt/2)],2) ) );
  }// ijk
  free(corr_2pts_0m1_th1);
  free(corr_2pts_0m2_th2);
  free(corr_rest_2pts_0ml);
  free(corr_rest_2pts_0mH);
  /////////////////// FINE 2pts CORRELATORS AT Tsep

  
  /////////////
  //         //
  // Ratios  //
  //         //
  /////////////

  double **ratio_ml_to_mH_Vi=(double**)malloc(sizeof(double)*(clusterfile));
  for(int ijk = 0; ijk < clusterfile; ijk++) ratio_ml_to_mH_Vi[ijk]=(double*)malloc(sizeof(double)*Nt);
  double **ratio_mH_to_ml_Vi=(double**)malloc(sizeof(double)*(clusterfile));
  for(int ijk = 0; ijk < clusterfile; ijk++) ratio_mH_to_ml_Vi[ijk]=(double*)malloc(sizeof(double)*Nt);
  double **ratio_Hl_sym_Vi=(double**)malloc(sizeof(double)*(clusterfile));
  for(int ijk = 0; ijk < clusterfile; ijk++) ratio_Hl_sym_Vi[ijk]=(double*)malloc(sizeof(double)*Nt);
  double **ratio_Hl_sym_Vi_av=(double**)malloc(sizeof(double)*(clusterfile));
  for(int ijk = 0; ijk < clusterfile; ijk++) ratio_Hl_sym_Vi_av[ijk]=(double*)malloc(sizeof(double)*Nt);

  double **ratio_ml_to_mH_V0=(double**)malloc(sizeof(double)*(clusterfile));
  for(int ijk = 0; ijk < clusterfile; ijk++) ratio_ml_to_mH_V0[ijk]=(double*)malloc(sizeof(double)*Nt);
  double **ratio_mH_to_ml_V0=(double**)malloc(sizeof(double)*(clusterfile));
  for(int ijk = 0; ijk < clusterfile; ijk++) ratio_mH_to_ml_V0[ijk]=(double*)malloc(sizeof(double)*Nt);
  double **ratio_Hl_sym_V0=(double**)malloc(sizeof(double)*(clusterfile));
  for(int ijk = 0; ijk < clusterfile; ijk++) ratio_Hl_sym_V0[ijk]=(double*)malloc(sizeof(double)*Nt);
  double **ratio_Hl_sym_V0_av=(double**)malloc(sizeof(double)*(clusterfile));
  for(int ijk = 0; ijk < clusterfile; ijk++) ratio_Hl_sym_V0_av[ijk]=(double*)malloc(sizeof(double)*Nt);

  double **ratio_ml_to_mH_S0=(double**)malloc(sizeof(double)*(clusterfile));
  for(int ijk = 0; ijk < clusterfile; ijk++) ratio_ml_to_mH_S0[ijk]=(double*)malloc(sizeof(double)*Nt);
  double **ratio_mH_to_ml_S0=(double**)malloc(sizeof(double)*(clusterfile));
  for(int ijk = 0; ijk < clusterfile; ijk++) ratio_mH_to_ml_S0[ijk]=(double*)malloc(sizeof(double)*Nt);
  double **ratio_Hl_sym_S0=(double**)malloc(sizeof(double)*(clusterfile));
  for(int ijk = 0; ijk < clusterfile; ijk++) ratio_Hl_sym_S0[ijk]=(double*)malloc(sizeof(double)*Nt);
  double **ratio_Hl_sym_S0_av=(double**)malloc(sizeof(double)*(clusterfile));
  for(int ijk = 0; ijk < clusterfile; ijk++) ratio_Hl_sym_S0_av[ijk]=(double*)malloc(sizeof(double)*Nt);
  
  double **ratio_Hl_rest_sym_V0_av=(double**)malloc(sizeof(double)*(clusterfile));
  for(int ijk = 0; ijk < clusterfile; ijk++) ratio_Hl_rest_sym_V0_av[ijk]=(double*)malloc(sizeof(double)*Nt);
  double **ratio_Hl_rest_sym_S0_av=(double**)malloc(sizeof(double)*(clusterfile));
  for(int ijk = 0; ijk < clusterfile; ijk++) ratio_Hl_rest_sym_S0_av[ijk]=(double*)malloc(sizeof(double)*Nt);
  
  
  for(int ijk = 0; ijk < clusterfile; ijk++){
    for(int t = 0; t <= Tsep_num[ibeta]; t++){

      ratio_ml_to_mH_Vi[ijk][t]  = ( corr_3pts_ml_to_mH_Vi[(ijk*Nt+t)]*corr_3pts_ml_to_mH_Vi[(ijk*Nt+Tsep_num[ibeta]-t)] )/( corr_3pts_ml_to_ml_Vi[(ijk*Nt+t)]*corr_3pts_mH_to_mH_Vi[(ijk*Nt+t)]);
      ratio_mH_to_ml_Vi[ijk][t]  = ( corr_3pts_mH_to_ml_Vi[(ijk*Nt+t)]*corr_3pts_mH_to_ml_Vi[(ijk*Nt+Tsep_num[ibeta]-t)] )/( corr_3pts_ml_to_ml_Vi[(ijk*Nt+t)]*corr_3pts_mH_to_mH_Vi[(ijk*Nt+t)]);
      ratio_Hl_sym_Vi[ijk][t]    = ( corr_3pts_mH_to_ml_Vi[(ijk*Nt+t)]*corr_3pts_ml_to_mH_Vi[(ijk*Nt+t)] )/( corr_3pts_ml_to_ml_Vi[(ijk*Nt+t)]*corr_3pts_mH_to_mH_Vi[(ijk*Nt+t)]);

      ratio_Hl_sym_Vi_av[ijk][t] = ( corr_3pts_mH_to_ml_Vi_av[(ijk*Nt+t)]*corr_3pts_ml_to_mH_Vi_av[(ijk*Nt+t)] )/( corr_3pts_ml_to_ml_Vi_av[(ijk*Nt+t)]*corr_3pts_mH_to_mH_Vi_av[(ijk*Nt+t)]);
      
      if(ith1 == 3){

        ratio_ml_to_mH_Vi[ijk][t]  = ( corr_3pts_ml_to_mH_Vi[(ijk*Nt+t)]*corr_3pts_ml_to_mH_Vi[(ijk*Nt+Tsep_num[ibeta]-t)] )/( corr_3pts_ml_to_ml_V0[(ijk*Nt+t)]*corr_3pts_mH_to_mH_Vi[(ijk*Nt+t)]);
        ratio_mH_to_ml_Vi[ijk][t]  = ( corr_3pts_mH_to_ml_Vi[(ijk*Nt+t)]*corr_3pts_mH_to_ml_Vi[(ijk*Nt+Tsep_num[ibeta]-t)] )/( corr_3pts_ml_to_ml_V0[(ijk*Nt+t)]*corr_3pts_mH_to_mH_Vi[(ijk*Nt+t)]);
        ratio_Hl_sym_Vi[ijk][t]    = ( corr_3pts_mH_to_ml_Vi[(ijk*Nt+t)]*corr_3pts_ml_to_mH_Vi[(ijk*Nt+t)] )/( corr_3pts_ml_to_ml_V0[(ijk*Nt+t)]*corr_3pts_mH_to_mH_Vi[(ijk*Nt+t)]);

        ratio_Hl_sym_Vi_av[ijk][t] = ( corr_3pts_mH_to_ml_Vi_av[(ijk*Nt+t)]*corr_3pts_ml_to_mH_Vi_av[(ijk*Nt+t)] )/( corr_3pts_ml_to_ml_V0_av[(ijk*Nt+t)]*corr_3pts_mH_to_mH_Vi_av[(ijk*Nt+t)]);
	
      }// if
	
      if(ith2 == 3){

        ratio_ml_to_mH_Vi[ijk][t]  = ( corr_3pts_ml_to_mH_Vi[(ijk*Nt+t)]*corr_3pts_ml_to_mH_Vi[(ijk*Nt+Tsep_num[ibeta]-t)] )/( corr_3pts_ml_to_ml_Vi[(ijk*Nt+t)]*corr_3pts_mH_to_mH_V0[(ijk*Nt+t)]);
        ratio_mH_to_ml_Vi[ijk][t]  = ( corr_3pts_mH_to_ml_Vi[(ijk*Nt+t)]*corr_3pts_mH_to_ml_Vi[(ijk*Nt+Tsep_num[ibeta]-t)] )/( corr_3pts_ml_to_ml_Vi[(ijk*Nt+t)]*corr_3pts_mH_to_mH_V0[(ijk*Nt+t)]);
        ratio_Hl_sym_Vi[ijk][t]    = ( corr_3pts_mH_to_ml_Vi[(ijk*Nt+t)]*corr_3pts_ml_to_mH_Vi[(ijk*Nt+t)] )/( corr_3pts_ml_to_ml_Vi[(ijk*Nt+t)]*corr_3pts_mH_to_mH_V0[(ijk*Nt+t)]);

        ratio_Hl_sym_Vi_av[ijk][t] = ( corr_3pts_mH_to_ml_Vi_av[(ijk*Nt+t)]*corr_3pts_ml_to_mH_Vi_av[(ijk*Nt+t)] )/( corr_3pts_ml_to_ml_Vi_av[(ijk*Nt+t)]*corr_3pts_mH_to_mH_V0_av[(ijk*Nt+t)]);
	
      }// if
      
      if(ith1 == 3 && ith2 == 3){

        ratio_ml_to_mH_Vi[ijk][t]  = ( corr_3pts_ml_to_mH_Vi[(ijk*Nt+t)]*corr_3pts_ml_to_mH_Vi[(ijk*Nt+Tsep_num[ibeta]-t)] )/( corr_3pts_ml_to_ml_V0[(ijk*Nt+t)]*corr_3pts_mH_to_mH_V0[(ijk*Nt+t)]);
        ratio_mH_to_ml_Vi[ijk][t]  = ( corr_3pts_mH_to_ml_Vi[(ijk*Nt+t)]*corr_3pts_mH_to_ml_Vi[(ijk*Nt+Tsep_num[ibeta]-t)] )/( corr_3pts_ml_to_ml_V0[(ijk*Nt+t)]*corr_3pts_mH_to_mH_V0[(ijk*Nt+t)]);
        ratio_Hl_sym_Vi[ijk][t]    = ( corr_3pts_mH_to_ml_Vi[(ijk*Nt+t)]*corr_3pts_ml_to_mH_Vi[(ijk*Nt+t)] )/( corr_3pts_ml_to_ml_V0[(ijk*Nt+t)]*corr_3pts_mH_to_mH_V0[(ijk*Nt+t)]);

        ratio_Hl_sym_Vi_av[ijk][t] = ( corr_3pts_mH_to_ml_Vi_av[(ijk*Nt+t)]*corr_3pts_ml_to_mH_Vi_av[(ijk*Nt+t)] )/( corr_3pts_ml_to_ml_V0_av[(ijk*Nt+t)]*corr_3pts_mH_to_mH_V0_av[(ijk*Nt+t)]);
	
      }// if
      
      ratio_ml_to_mH_V0[ijk][t]  = ( corr_3pts_ml_to_mH_V0[(ijk*Nt+t)]*corr_3pts_ml_to_mH_V0[(ijk*Nt+Tsep_num[ibeta]-t)] )/( corr_3pts_ml_to_ml_V0[(ijk*Nt+t)]*corr_3pts_mH_to_mH_V0[(ijk*Nt+t)]);
      ratio_mH_to_ml_V0[ijk][t]  = ( corr_3pts_mH_to_ml_V0[(ijk*Nt+t)]*corr_3pts_mH_to_ml_V0[(ijk*Nt+Tsep_num[ibeta]-t)] )/( corr_3pts_ml_to_ml_V0[(ijk*Nt+t)]*corr_3pts_mH_to_mH_V0[(ijk*Nt+t)]);
      ratio_Hl_sym_V0[ijk][t]    = ( corr_3pts_mH_to_ml_V0[(ijk*Nt+t)]*corr_3pts_ml_to_mH_V0[(ijk*Nt+t)] )/( corr_3pts_ml_to_ml_V0[(ijk*Nt+t)]*corr_3pts_mH_to_mH_V0[(ijk*Nt+t)]);

      ratio_Hl_sym_V0_av[ijk][t] = ( corr_3pts_mH_to_ml_V0_av[(ijk*Nt+t)]*corr_3pts_ml_to_mH_V0_av[(ijk*Nt+t)] )/( corr_3pts_ml_to_ml_V0_av[(ijk*Nt+t)]*corr_3pts_mH_to_mH_V0_av[(ijk*Nt+t)]);
      
      ratio_ml_to_mH_S0[ijk][t]  = (corr_3pts_ml_to_mH_S0[ijk*Nt+t]*corr_3pts_ml_to_mH_S0[ijk*Nt+Tsep_num[ibeta]-t])/(tilde_corr_2pts_Tsep_0m1_th1[ijk]*tilde_corr_2pts_Tsep_0m2_th2[ijk]);
      ratio_mH_to_ml_S0[ijk][t]  = (corr_3pts_mH_to_ml_S0[ijk*Nt+t]*corr_3pts_mH_to_ml_S0[ijk*Nt+Tsep_num[ibeta]-t])/(tilde_corr_2pts_Tsep_0m1_th1[ijk]*tilde_corr_2pts_Tsep_0m2_th2[ijk]);
      ratio_Hl_sym_S0[ijk][t]    = (corr_3pts_mH_to_ml_S0[ijk*Nt+t]*corr_3pts_ml_to_mH_S0[ijk*Nt+t])/(tilde_corr_2pts_Tsep_0m1_th1[ijk]*tilde_corr_2pts_Tsep_0m2_th2[ijk]);

      ratio_Hl_sym_S0_av[ijk][t] = (corr_3pts_mH_to_ml_S0_av[ijk*Nt+t]*corr_3pts_ml_to_mH_S0_av[ijk*Nt+t])/(tilde_corr_2pts_Tsep_0m1_th1[ijk]*tilde_corr_2pts_Tsep_0m2_th2[ijk]);
      
      ratio_Hl_rest_sym_V0_av[ijk][t] = ( corr_3pts_rest_ml_to_mH_V0_av[ijk*Nt+t]*corr_3pts_rest_mH_to_ml_V0_av[ijk*Nt+t]  )/(corr_3pts_rest_ml_to_ml_V0_av[ijk*Nt+t]*corr_3pts_rest_mH_to_mH_V0_av[ijk*Nt+t]  );
      ratio_Hl_rest_sym_S0_av[ijk][t] = ( corr_3pts_rest_ml_to_mH_S0_av[ijk*Nt+t]*corr_3pts_rest_mH_to_ml_S0_av[ijk*Nt+t]  )/(tilde_corr_rest_2pts_Tsep_0ml[ijk]*tilde_corr_rest_2pts_Tsep_0mH[ijk]);
      
    }// t
  }// ijk
  free(corr_3pts_ml_to_mH_V0);  
  free(corr_3pts_ml_to_ml_V0);
  free(corr_3pts_mH_to_mH_V0);
  free(corr_3pts_ml_to_mH_V0_av);
  free(corr_3pts_ml_to_ml_V0_av);
  free(corr_3pts_mH_to_mH_V0_av);

  free(corr_3pts_ml_to_mH_Vi);  
  free(corr_3pts_ml_to_ml_Vi);
  free(corr_3pts_mH_to_mH_Vi);
  free(corr_3pts_ml_to_mH_Vi_av);
  free(corr_3pts_ml_to_ml_Vi_av);
  free(corr_3pts_mH_to_mH_Vi_av);
  
  free(corr_3pts_mH_to_ml_S0_av);
  free(corr_3pts_rest_ml_to_mH_V0_av);
  free(corr_3pts_rest_ml_to_ml_V0_av);
  free(corr_3pts_rest_mH_to_mH_V0_av);
  free(corr_3pts_rest_ml_to_mH_S0_av);  
  //////////// FINE RATIOS
  
  
  /////////////////////////////
  //                         //
  // COSTRUISCO Vi, V0 ED S0 // 
  //                         //
  /////////////////////////////
  
  double coeff_V0, coeff_Vi, coeff_S0, coeff_V0_rest, coeff_S0_rest;

  ///////// DICHIARAZIONE ELEMENTI DI MATRICE
  double *Vi_ml_to_mH=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  double *Vi_mH_to_ml=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  double *Vi_Hl_sym=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  double *Vi_Hl_sym_av=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  
  double *V0_ml_to_mH=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  double *V0_mH_to_ml=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  double *V0_Hl_sym=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  double *V0_Hl_sym_av=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  
  double *S0_ml_to_mH=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  double *S0_mH_to_ml=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  double *S0_Hl_sym=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  double *S0_Hl_sym_av=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));

  double *V0_Hl_rest_sym_av=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  
  double *S0_Hl_rest_sym_av=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  ///////// FINE DICHIARAZIONE ELEMENTI DI MATRICE

  
  for(int ijk=0; ijk < clusterfile; ijk++){
    
    coeff_V0 = 4*energy_2pts_m1[ijk]*energy_2pts_m2[ijk];
    coeff_S0 = coeff_V0;
    coeff_Vi = 4*momentum_1*momentum_2;
    
    coeff_V0_rest = 4*mass_2pts_ml[ijk]*mass_2pts_mH[ijk];
    coeff_S0_rest = coeff_V0_rest;
    
    if(ith1 == 3){
      coeff_Vi = 4*energy_2pts_m1[ijk]*momentum_2;
    }
    if(ith2 == 3){
      coeff_Vi = 4*energy_2pts_m2[ijk]*momentum_1;
    }
    if(ith1 == 3 && ith2 == 3){
      coeff_Vi = 4*energy_2pts_m1[ijk]*energy_2pts_m2[ijk];
    }
    
    for(int t = 0; t <= Tsep_num[ibeta]; t++){
      
      //   OCCHIO AL SEGNO - DAVANTI AI VETTORIALI
      // IL SEGNO ASSOCIATO AGLI EL. DI MATRICE VETTORIALI È SEMPRE QUELLO DI HEAVY_TO_LIGHT !!!
      
      Vi_ml_to_mH[ijk*Nt+t]  = -sqrt(fabs(ratio_ml_to_mH_Vi[ijk][t]*coeff_Vi))*( corr_3pts_mH_to_ml_Vi[(ijk*Nt)+(Tsep_num[ibeta]/2)]/fabs(corr_3pts_mH_to_ml_Vi[(ijk*Nt)+(Tsep_num[ibeta]/2)]) );
      Vi_mH_to_ml[ijk*Nt+t]  = -sqrt(fabs(ratio_mH_to_ml_Vi[ijk][t]*coeff_Vi))*( corr_3pts_mH_to_ml_Vi[(ijk*Nt)+(Tsep_num[ibeta]/2)]/fabs(corr_3pts_mH_to_ml_Vi[(ijk*Nt)+(Tsep_num[ibeta]/2)]) );
      Vi_Hl_sym[ijk*Nt+t]    = -sqrt(fabs(ratio_Hl_sym_Vi[ijk][t]*coeff_Vi))*( corr_3pts_mH_to_ml_Vi[(ijk*Nt)+(Tsep_num[ibeta]/2)]/fabs(corr_3pts_mH_to_ml_Vi[(ijk*Nt)+(Tsep_num[ibeta]/2)]) );
      Vi_Hl_sym_av[ijk*Nt+t] = -sqrt(fabs(ratio_Hl_sym_Vi_av[ijk][t]*coeff_Vi))*( corr_3pts_mH_to_ml_Vi_av[(ijk*Nt)+(Tsep_num[ibeta]/2)]/fabs(corr_3pts_mH_to_ml_Vi_av[(ijk*Nt)+(Tsep_num[ibeta]/2)]) );
      
      V0_ml_to_mH[ijk*Nt+t]  = -sqrt(fabs(ratio_ml_to_mH_V0[ijk][t]*coeff_V0))*( corr_3pts_mH_to_ml_V0[(ijk*Nt)+(Tsep_num[ibeta]/2)]/fabs(corr_3pts_mH_to_ml_V0[(ijk*Nt)+(Tsep_num[ibeta]/2)]) );
      V0_mH_to_ml[ijk*Nt+t]  = -sqrt(fabs(ratio_mH_to_ml_V0[ijk][t]*coeff_V0))*( corr_3pts_mH_to_ml_V0[(ijk*Nt)+(Tsep_num[ibeta]/2)]/fabs(corr_3pts_mH_to_ml_V0[(ijk*Nt)+(Tsep_num[ibeta]/2)]) );
      V0_Hl_sym[ijk*Nt+t]    = -sqrt(fabs(ratio_Hl_sym_V0[ijk][t]*coeff_V0))*( corr_3pts_mH_to_ml_V0[(ijk*Nt)+(Tsep_num[ibeta]/2)]/fabs(corr_3pts_mH_to_ml_V0[(ijk*Nt)+(Tsep_num[ibeta]/2)]) );
      V0_Hl_sym_av[ijk*Nt+t] = -sqrt(fabs(ratio_Hl_sym_V0_av[ijk][t]*coeff_V0))*( corr_3pts_mH_to_ml_V0_av[(ijk*Nt)+(Tsep_num[ibeta]/2)]/fabs(corr_3pts_mH_to_ml_V0_av[(ijk*Nt)+(Tsep_num[ibeta]/2)]) );
      
      V0_Hl_rest_sym_av[ijk*Nt+t] = -sqrt(fabs(ratio_Hl_rest_sym_V0_av[ijk][t]*coeff_V0_rest))*(corr_3pts_rest_mH_to_ml_V0_av[(ijk*Nt)+(Tsep_num[ibeta]/2)]/fabs(corr_3pts_rest_mH_to_ml_V0_av[(ijk*Nt)+(Tsep_num[ibeta]/2)]) );
      
      
      S0_ml_to_mH[ijk*Nt+t]  = sqrt(fabs(ratio_ml_to_mH_S0[ijk][t]*coeff_S0))*( corr_3pts_ml_to_mH_S0[(ijk*Nt)+(Tsep_num[ibeta]/2)]/fabs(corr_3pts_ml_to_mH_S0[(ijk*Nt)+(Tsep_num[ibeta]/2)]) );
      S0_mH_to_ml[ijk*Nt+t]  = sqrt(fabs(ratio_mH_to_ml_S0[ijk][t]*coeff_S0))*( corr_3pts_mH_to_ml_S0[(ijk*Nt)+(Tsep_num[ibeta]/2)]/fabs(corr_3pts_mH_to_ml_S0[(ijk*Nt)+(Tsep_num[ibeta]/2)]) );
      S0_Hl_sym[ijk*Nt+t]    = sqrt(fabs(ratio_Hl_sym_S0[ijk][t]*coeff_S0))*( corr_3pts_ml_to_mH_S0[(ijk*Nt)+(Tsep_num[ibeta]/2)]/fabs(corr_3pts_ml_to_mH_S0[(ijk*Nt)+(Tsep_num[ibeta]/2)]) );
      S0_Hl_sym_av[ijk*Nt+t] = sqrt(fabs(ratio_Hl_sym_S0_av[ijk][t]*coeff_S0))*( corr_3pts_ml_to_mH_S0_av[(ijk*Nt)+(Tsep_num[ibeta]/2)]/fabs(corr_3pts_ml_to_mH_S0_av[(ijk*Nt)+(Tsep_num[ibeta]/2)]) );
      
      S0_Hl_rest_sym_av[ijk*Nt+t] = sqrt(fabs(ratio_Hl_rest_sym_S0_av[ijk][t]*coeff_S0_rest))*(corr_3pts_rest_mH_to_ml_S0_av[(ijk*Nt)+(Tsep_num[ibeta]/2)]/fabs(corr_3pts_rest_mH_to_ml_S0_av[(ijk*Nt)+(Tsep_num[ibeta]/2)]) );
      
    }// t
  }// ijk
  free(corr_3pts_mH_to_ml_V0);
  free(corr_3pts_mH_to_ml_V0_av);
  free(corr_3pts_mH_to_ml_Vi);
  free(corr_3pts_mH_to_ml_Vi_av);
  free(corr_3pts_ml_to_mH_S0);
  free(corr_3pts_mH_to_ml_S0);
  free(corr_3pts_ml_to_mH_S0_av);
  free(corr_3pts_rest_mH_to_ml_V0_av);
  free(corr_3pts_rest_mH_to_ml_S0_av);
  
  free(ratio_ml_to_mH_V0);
  free(ratio_mH_to_ml_V0);
  free(ratio_Hl_sym_V0);
  free(ratio_Hl_sym_V0_av);

  free(ratio_ml_to_mH_Vi);
  free(ratio_mH_to_ml_Vi);
  free(ratio_Hl_sym_Vi);
  free(ratio_Hl_sym_Vi_av);

  free(ratio_ml_to_mH_S0);
  free(ratio_mH_to_ml_S0);
  free(ratio_Hl_sym_S0);
  free(ratio_Hl_sym_S0_av);
  free(ratio_Hl_rest_sym_V0_av);
  free(ratio_Hl_rest_sym_S0_av);
    
  /////// SIGMA Vi
  double *sigma_Vi_ml_to_mH=(double*)malloc(sizeof(double)*Nt), *sigma_Vi_mH_to_ml=(double*)malloc(sizeof(double)*Nt);
  double *sigma_Vi_Hl_sym=(double*)malloc(sizeof(double)*Nt), *sigma_Vi_Hl_sym_av=(double*)malloc(sizeof(double)*Nt);

  sigma_matrix_el( sigma_Vi_ml_to_mH, Vi_ml_to_mH, Tsep_num[ibeta], Nt, clusterfile);
  sigma_matrix_el( sigma_Vi_mH_to_ml, Vi_mH_to_ml, Tsep_num[ibeta], Nt, clusterfile);
  sigma_matrix_el( sigma_Vi_Hl_sym, Vi_Hl_sym, Tsep_num[ibeta], Nt, clusterfile);
  sigma_matrix_el( sigma_Vi_Hl_sym_av, Vi_Hl_sym_av, Tsep_num[ibeta], Nt, clusterfile);

  /////// SIGMA V0
  double *sigma_V0_ml_to_mH=(double*)malloc(sizeof(double)*Nt), *sigma_V0_mH_to_ml=(double*)malloc(sizeof(double)*Nt);
  double *sigma_V0_Hl_sym=(double*)malloc(sizeof(double)*Nt), *sigma_V0_Hl_sym_av=(double*)malloc(sizeof(double)*Nt);
  double *sigma_V0_Hl_rest_sym_av=(double*)malloc(sizeof(double)*Nt);

  sigma_matrix_el( sigma_V0_ml_to_mH, V0_ml_to_mH, Tsep_num[ibeta], Nt, clusterfile);
  sigma_matrix_el( sigma_V0_mH_to_ml, V0_mH_to_ml, Tsep_num[ibeta], Nt, clusterfile);
  sigma_matrix_el( sigma_V0_Hl_sym, V0_Hl_sym, Tsep_num[ibeta], Nt, clusterfile);
  sigma_matrix_el( sigma_V0_Hl_sym_av, V0_Hl_sym_av, Tsep_num[ibeta], Nt, clusterfile);
  sigma_matrix_el( sigma_V0_Hl_rest_sym_av, V0_Hl_rest_sym_av, Tsep_num[ibeta], Nt, clusterfile);

  /////// SIGMA V0
  double *sigma_S0_ml_to_mH=(double*)malloc(sizeof(double)*Nt), *sigma_S0_mH_to_ml=(double*)malloc(sizeof(double)*Nt);
  double *sigma_S0_Hl_sym=(double*)malloc(sizeof(double)*Nt), *sigma_S0_Hl_sym_av=(double*)malloc(sizeof(double)*Nt);
  double *sigma_S0_Hl_rest_sym_av=(double*)malloc(sizeof(double)*Nt);
  
  sigma_matrix_el( sigma_S0_ml_to_mH, S0_ml_to_mH, Tsep_num[ibeta], Nt, clusterfile);
  sigma_matrix_el( sigma_S0_mH_to_ml, S0_mH_to_ml, Tsep_num[ibeta], Nt, clusterfile);
  sigma_matrix_el( sigma_S0_Hl_sym, S0_Hl_sym, Tsep_num[ibeta], Nt, clusterfile);
  sigma_matrix_el( sigma_S0_Hl_sym_av, S0_Hl_sym_av, Tsep_num[ibeta], Nt, clusterfile);
  sigma_matrix_el( sigma_S0_Hl_rest_sym_av, S0_Hl_rest_sym_av, Tsep_num[ibeta], Nt, clusterfile);

  
  ///////// FACCIO IL FIT A COSTANTE

  //// V0 ////
  double *V0_fit_ml_to_mH=(double*)malloc(sizeof(double)*clusterfile), *V0_fit_mH_to_ml=(double*)malloc(sizeof(double)*clusterfile);
  double *V0_fit_Hl_sym=(double*)malloc(sizeof(double)*clusterfile), *V0_fit_Hl_sym_av=(double*)malloc(sizeof(double)*clusterfile);
  double *V0_fit_Hl_rest_sym_av=(double*)malloc(sizeof(double)*clusterfile);

  constant_fit_scorrelated( V0_fit_ml_to_mH, V0_ml_to_mH, sigma_V0_ml_to_mH, clusterfile, T_max[ibeta], ( Tsep_num[ibeta]/2 - fit_interval[ibeta] ), ( Tsep_num[ibeta]/2 + fit_interval[ibeta] ));
  constant_fit_scorrelated( V0_fit_mH_to_ml, V0_mH_to_ml, sigma_V0_mH_to_ml, clusterfile, T_max[ibeta], ( Tsep_num[ibeta]/2 - fit_interval[ibeta] ), ( Tsep_num[ibeta]/2 + fit_interval[ibeta] ));
  constant_fit_scorrelated( V0_fit_Hl_sym, V0_Hl_sym, sigma_V0_Hl_sym, clusterfile, T_max[ibeta], ( Tsep_num[ibeta]/2 - fit_interval[ibeta] ), ( Tsep_num[ibeta]/2 + fit_interval[ibeta] ));
  constant_fit_scorrelated( V0_fit_Hl_sym_av, V0_Hl_sym_av, sigma_V0_Hl_sym_av, clusterfile, T_max[ibeta], ( Tsep_num[ibeta]/2 - fit_interval[ibeta] ), ( Tsep_num[ibeta]/2 + fit_interval[ibeta] ));  
  constant_fit_scorrelated( V0_fit_Hl_rest_sym_av, V0_Hl_rest_sym_av, sigma_V0_Hl_rest_sym_av, clusterfile, T_max[ibeta], ( Tsep_num[ibeta]/2 - fit_interval[ibeta] ), ( Tsep_num[ibeta]/2 + fit_interval[ibeta] ));

  //// Vi ////
  double *Vi_fit_ml_to_mH=(double*)malloc(sizeof(double)*clusterfile), *Vi_fit_mH_to_ml=(double*)malloc(sizeof(double)*clusterfile);
  double *Vi_fit_Hl_sym=(double*)malloc(sizeof(double)*clusterfile), *Vi_fit_Hl_sym_av=(double*)malloc(sizeof(double)*clusterfile);

  constant_fit_scorrelated( Vi_fit_ml_to_mH, Vi_ml_to_mH, sigma_Vi_ml_to_mH, clusterfile, T_max[ibeta], ( Tsep_num[ibeta]/2 - fit_interval[ibeta] ), ( Tsep_num[ibeta]/2 + fit_interval[ibeta] ));
  constant_fit_scorrelated( Vi_fit_mH_to_ml, Vi_mH_to_ml, sigma_Vi_mH_to_ml, clusterfile, T_max[ibeta], ( Tsep_num[ibeta]/2 - fit_interval[ibeta] ), ( Tsep_num[ibeta]/2 + fit_interval[ibeta] ));
  constant_fit_scorrelated( Vi_fit_Hl_sym, Vi_Hl_sym, sigma_Vi_Hl_sym, clusterfile, T_max[ibeta], ( Tsep_num[ibeta]/2 - fit_interval[ibeta] ), ( Tsep_num[ibeta]/2 + fit_interval[ibeta] ));
  constant_fit_scorrelated( Vi_fit_Hl_sym_av, Vi_Hl_sym_av, sigma_Vi_Hl_sym_av, clusterfile, T_max[ibeta], ( Tsep_num[ibeta]/2 - fit_interval[ibeta] ), ( Tsep_num[ibeta]/2 + fit_interval[ibeta] ));  

  //// S0 ////
  double *S0_fit_ml_to_mH=(double*)malloc(sizeof(double)*clusterfile), *S0_fit_mH_to_ml=(double*)malloc(sizeof(double)*clusterfile);
  double *S0_fit_Hl_sym=(double*)malloc(sizeof(double)*clusterfile), *S0_fit_Hl_sym_av=(double*)malloc(sizeof(double)*clusterfile);
  double *S0_fit_Hl_rest_sym_av=(double*)malloc(sizeof(double)*clusterfile);

  constant_fit_scorrelated( S0_fit_ml_to_mH, S0_ml_to_mH, sigma_S0_ml_to_mH, clusterfile, T_max[ibeta], ( Tsep_num[ibeta]/2 - fit_interval[ibeta] ), ( Tsep_num[ibeta]/2 + fit_interval[ibeta] ));
  constant_fit_scorrelated( S0_fit_mH_to_ml, S0_mH_to_ml, sigma_S0_mH_to_ml, clusterfile, T_max[ibeta], ( Tsep_num[ibeta]/2 - fit_interval[ibeta] ), ( Tsep_num[ibeta]/2 + fit_interval[ibeta] ));
  constant_fit_scorrelated( S0_fit_Hl_sym, S0_Hl_sym, sigma_S0_Hl_sym, clusterfile, T_max[ibeta], ( Tsep_num[ibeta]/2 - fit_interval[ibeta] ), ( Tsep_num[ibeta]/2 + fit_interval[ibeta] ));
  constant_fit_scorrelated( S0_fit_Hl_sym_av, S0_Hl_sym_av, sigma_S0_Hl_sym_av, clusterfile, T_max[ibeta], ( Tsep_num[ibeta]/2 - fit_interval[ibeta] ), ( Tsep_num[ibeta]/2 + fit_interval[ibeta] ));  
  constant_fit_scorrelated( S0_fit_Hl_rest_sym_av, S0_Hl_rest_sym_av, sigma_S0_Hl_rest_sym_av, clusterfile, T_max[ibeta], ( Tsep_num[ibeta]/2 - fit_interval[ibeta] ), ( Tsep_num[ibeta]/2 + fit_interval[ibeta] ));
  
  ////////////////////////
  //                    //
  // SIGMA MEDIA PESATA //
  //                    //
  ////////////////////////

  double sigma_Vi_fit_ml_to_mH, sigma_Vi_fit_mH_to_ml, sigma_Vi_fit_Hl_sym, sigma_Vi_fit_Hl_sym_av;
  double sigma_V0_fit_ml_to_mH, sigma_V0_fit_mH_to_ml, sigma_V0_fit_Hl_sym, sigma_V0_fit_Hl_sym_av;
  double sigma_S0_fit_ml_to_mH, sigma_S0_fit_mH_to_ml, sigma_S0_fit_Hl_sym, sigma_S0_fit_Hl_sym_av;
  double sigma_V0_fit_Hl_rest_sym_av, sigma_S0_fit_Hl_rest_sym_av;
  
  sigma_Vi_fit_ml_to_mH = sigma_JK(Vi_fit_ml_to_mH, clusterfile);
  sigma_Vi_fit_mH_to_ml = sigma_JK(Vi_fit_mH_to_ml, clusterfile);
  sigma_Vi_fit_Hl_sym = sigma_JK(Vi_fit_Hl_sym, clusterfile);
  sigma_Vi_fit_Hl_sym_av = sigma_JK(Vi_fit_Hl_sym_av, clusterfile);
  
  sigma_V0_fit_ml_to_mH = sigma_JK(V0_fit_ml_to_mH, clusterfile);
  sigma_V0_fit_mH_to_ml = sigma_JK(V0_fit_mH_to_ml, clusterfile);
  sigma_V0_fit_Hl_sym = sigma_JK(V0_fit_Hl_sym, clusterfile);
  sigma_V0_fit_Hl_sym_av = sigma_JK(V0_fit_Hl_sym_av, clusterfile);
  sigma_V0_fit_Hl_rest_sym_av = sigma_JK(V0_fit_Hl_rest_sym_av, clusterfile);
  
  sigma_S0_fit_ml_to_mH = sigma_JK(S0_fit_ml_to_mH, clusterfile);
  sigma_S0_fit_mH_to_ml = sigma_JK(S0_fit_mH_to_ml, clusterfile);
  sigma_S0_fit_Hl_sym = sigma_JK(S0_fit_Hl_sym, clusterfile);
  sigma_S0_fit_Hl_sym_av = sigma_JK(S0_fit_Hl_sym_av, clusterfile);
  sigma_S0_fit_Hl_rest_sym_av = sigma_JK(S0_fit_Hl_rest_sym_av, clusterfile);


  
  //////////////////////
  //                  //
  //  PLATEAUX GRACE  //
  //                  //
  //////////////////////

  ////////////////
  //            //
  //     Vi     //
  //            //
  ////////////////
  
  ////////// LIGHT_TO_HEAVY
  plot_plateau( Vi_fit_ml_to_mH, sigma_Vi_fit_ml_to_mH, Vi_ml_to_mH, sigma_Vi_ml_to_mH, fit_interval, clusterfile, Nt, Tsep_num, strategy, dir_S, dir_E, beta_V, mu_sea_1, op, mu_sea_2, m_, th_, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 0, 1);

  free(Vi_ml_to_mH);
  free(sigma_Vi_ml_to_mH);
  
  ////////// HEAVY_TO_LIGHT
  plot_plateau( Vi_fit_mH_to_ml, sigma_Vi_fit_mH_to_ml, Vi_mH_to_ml, sigma_Vi_mH_to_ml, fit_interval, clusterfile, Nt, Tsep_num, strategy, dir_S, dir_E, beta_V, mu_sea_1, op, mu_sea_2, m_, th_, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 1, 1);

  free(Vi_mH_to_ml);
  free(sigma_Vi_mH_to_ml);
  
  ////////// HL_SYM
  plot_plateau( Vi_fit_Hl_sym, sigma_Vi_fit_Hl_sym, Vi_Hl_sym, sigma_Vi_Hl_sym, fit_interval, clusterfile, Nt, Tsep_num, strategy, dir_S, dir_E, beta_V, mu_sea_1, op, mu_sea_2, m_, th_, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 2, 1);

  free(Vi_Hl_sym);
  free(sigma_Vi_Hl_sym);
  
  ////////// HL_SYM_AV
  plot_plateau( Vi_fit_Hl_sym_av, sigma_Vi_fit_Hl_sym_av, Vi_Hl_sym_av, sigma_Vi_Hl_sym_av, fit_interval, clusterfile, Nt, Tsep_num, strategy, dir_S, dir_E, beta_V, mu_sea_1, op, mu_sea_2, m_, th_, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 3, 1);

  free(Vi_Hl_sym_av);
  free(sigma_Vi_Hl_sym_av);

  
  ////////////////
  //            //
  //     V0     //
  //            //
  ////////////////
  
  ////////// LIGHT_TO_HEAVY
  plot_plateau( V0_fit_ml_to_mH, sigma_V0_fit_ml_to_mH, V0_ml_to_mH, sigma_V0_ml_to_mH, fit_interval, clusterfile, Nt, Tsep_num, strategy, dir_S, dir_E, beta_V, mu_sea_1, op, mu_sea_2, m_, th_, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 0, 0);

  free(V0_ml_to_mH);
  free(sigma_V0_ml_to_mH);
  
  ////////// HEAVY_TO_LIGHT
  plot_plateau( V0_fit_mH_to_ml, sigma_V0_fit_mH_to_ml, V0_mH_to_ml, sigma_V0_mH_to_ml, fit_interval, clusterfile, Nt, Tsep_num, strategy, dir_S, dir_E, beta_V, mu_sea_1, op, mu_sea_2, m_, th_, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 1, 0);

  free(V0_mH_to_ml);
  free(sigma_V0_mH_to_ml);
  
  ////////// HL_SYM
  plot_plateau( V0_fit_Hl_sym, sigma_V0_fit_Hl_sym, V0_Hl_sym, sigma_V0_Hl_sym, fit_interval, clusterfile, Nt, Tsep_num, strategy, dir_S, dir_E, beta_V, mu_sea_1, op, mu_sea_2, m_, th_, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 2, 0);

  free(V0_Hl_sym);
  free(sigma_V0_Hl_sym);
  
  ////////// HL_SYM_AV
  plot_plateau( V0_fit_Hl_sym_av, sigma_V0_fit_Hl_sym_av, V0_Hl_sym_av, sigma_V0_Hl_sym_av, fit_interval, clusterfile, Nt, Tsep_num, strategy, dir_S, dir_E, beta_V, mu_sea_1, op, mu_sea_2, m_, th_, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 3, 0);

  free(V0_Hl_sym_av);
  free(sigma_V0_Hl_sym_av);

  
  ////////////////
  //            //
  //     S0     //
  //            //
  ////////////////
  
  ////////// LIGHT_TO_HEAVY
  plot_plateau( S0_fit_ml_to_mH, sigma_S0_fit_ml_to_mH, S0_ml_to_mH, sigma_S0_ml_to_mH, fit_interval, clusterfile, Nt, Tsep_num, strategy, dir_S, dir_E, beta_V, mu_sea_1, op, mu_sea_2, m_, th_, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 0, 2);

  free(S0_ml_to_mH);
  free(sigma_S0_ml_to_mH);
  
  ////////// HEAVY_TO_LIGHT
  plot_plateau( S0_fit_mH_to_ml, sigma_S0_fit_mH_to_ml, S0_mH_to_ml, sigma_S0_mH_to_ml, fit_interval, clusterfile, Nt, Tsep_num, strategy, dir_S, dir_E, beta_V, mu_sea_1, op, mu_sea_2, m_, th_, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 1, 2);

  free(S0_mH_to_ml);
  free(sigma_S0_mH_to_ml);
  
  ////////// HL_SYM
  plot_plateau( S0_fit_Hl_sym, sigma_S0_fit_Hl_sym, S0_Hl_sym, sigma_S0_Hl_sym, fit_interval, clusterfile, Nt, Tsep_num, strategy, dir_S, dir_E, beta_V, mu_sea_1, op, mu_sea_2, m_, th_, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 2, 2);

  free(S0_Hl_sym);
  free(sigma_S0_Hl_sym);
  
  ////////// HL_SYM_AV
  plot_plateau( S0_fit_Hl_sym_av, sigma_S0_fit_Hl_sym_av, S0_Hl_sym_av, sigma_S0_Hl_sym_av, fit_interval, clusterfile, Nt, Tsep_num, strategy, dir_S, dir_E, beta_V, mu_sea_1, op, mu_sea_2, m_, th_, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 3, 2);

  free(S0_Hl_sym_av);
  free(sigma_S0_Hl_sym_av);
  
  ///////// FINE PLATEAU
  

  ////////////////////////////
  //                        //
  //    OUTPUT Vi, V0 ED S0 //
  //                        //
  ////////////////////////////

  ////////////////
  //            //
  //     Vi     //
  //            //
  ////////////////
  
  ////////// LIGHT_TO_HEAVY  
  write_output_matrix_el_fit( Vi_fit_ml_to_mH, dir_S, dir_E, op, strategy, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 1, 0);
  
  ////////// HEAVY_TO_LIGHT
  write_output_matrix_el_fit( Vi_fit_mH_to_ml, dir_S, dir_E, op, strategy, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 1, 1);
  
  ////////// HL_SYM
  write_output_matrix_el_fit( Vi_fit_Hl_sym, dir_S, dir_E, op, strategy, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 1, 2);
  
  ////////// HL_SYM_AV
  write_output_matrix_el_fit( Vi_fit_Hl_sym_av, dir_S, dir_E, op, strategy, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 1, 3);
  
  
  ////////////////
  //            //
  //     V0     //
  //            //
  ////////////////

  ////////// LIGHT_TO_HEAVY  
  write_output_matrix_el_fit( V0_fit_ml_to_mH, dir_S, dir_E, op, strategy, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 0, 0);
  
  ////////// HEAVY_TO_LIGHT
  write_output_matrix_el_fit( V0_fit_mH_to_ml, dir_S, dir_E, op, strategy, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 0, 1);
  
  ////////// HL_SYM
  write_output_matrix_el_fit( V0_fit_Hl_sym, dir_S, dir_E, op, strategy, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 0, 2);
  
  ////////// HL_SYM_AV
  write_output_matrix_el_fit( V0_fit_Hl_sym_av, dir_S, dir_E, op, strategy, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 0, 3);
  
  
  ////////////////
  //            //
  //     S0     //
  //            //
  ////////////////
  
  ////////// LIGHT_TO_HEAVY  
  write_output_matrix_el_fit( S0_fit_ml_to_mH, dir_S, dir_E, op, strategy, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 2, 0);
  
  ////////// HEAVY_TO_LIGHT
  write_output_matrix_el_fit( S0_fit_mH_to_ml, dir_S, dir_E, op, strategy, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 2, 1);
  
  ////////// HL_SYM
  write_output_matrix_el_fit( S0_fit_Hl_sym, dir_S, dir_E, op, strategy, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 2, 2);
  
  ////////// HL_SYM_AV
  write_output_matrix_el_fit( S0_fit_Hl_sym_av, dir_S, dir_E, op, strategy, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 2, 3);

  ///////// FINE OUTPUT Vi, V0 ED S0


  ///////////////////////////////////////////////////
  //                                               //
  // CALCOLO I FATTORI DI FORMA PER FARE IL CHECK  //
  //                                               //
  ///////////////////////////////////////////////////

  double *q2=(double*)malloc(sizeof(double)*clusterfile);

  double *fplus_check=(double*)malloc(sizeof(double)*clusterfile);
  double *fzero_check=(double*)malloc(sizeof(double)*clusterfile);
  double *fminus_check=(double*)malloc(sizeof(double)*clusterfile);
  
  double sigma_fplus_check, sigma_fzero_check, sigma_fminus_check;
  
  for(int ijk = 0; ijk < clusterfile; ijk++){
    q2[ijk] = pow( energy_2pts_m2[ijk]-energy_2pts_m1[ijk], 2) - 3*pow( momentum_2-momentum_1, 2);
  }
  
  for(int ijk = 0; ijk < clusterfile; ijk++){
    
    fplus_check[ijk] = ( V0_fit_ml_to_mH[ijk]*(momentum_2-momentum_1) - Vi_fit_ml_to_mH[ijk]*(energy_2pts_m2[ijk]-energy_2pts_m1[ijk]) )/(2*(momentum_2*energy_2pts_m1[ijk]-momentum_1*energy_2pts_m2[ijk] )  ); 
    fminus_check[ijk] = ( V0_fit_ml_to_mH[ijk]*(momentum_2+momentum_1) - Vi_fit_ml_to_mH[ijk]*(energy_2pts_m2[ijk]+energy_2pts_m1[ijk]) )/(2*(-momentum_2*energy_2pts_m1[ijk]+momentum_1*energy_2pts_m2[ijk] )  );
    fzero_check[ijk] = fplus_check[ijk] + (fminus_check[ijk]*q2[ijk])/(pow(mass_2pts_mH[ijk],2)-pow(mass_2pts_ml[ijk],2));
    
    if(ith1 == 3 && ith2 == 3){
      
      fzero_check[ijk] = V0_fit_ml_to_mH[ijk]/(mass_2pts_mH[ijk] + mass_2pts_ml[ijk]);
      fplus_check[ijk] = 0;
      
    }// if
  }// ijk
  
  sigma_fplus_check = sigma_JK(fplus_check, clusterfile);
  sigma_fzero_check = sigma_JK(fzero_check, clusterfile);
  sigma_fminus_check = sigma_JK(fminus_check, clusterfile);
  //////// FINE CALCOLO I FATTORI DI FORMA PER FARE IL CHECK


  
  ///////////////////////////////////////////////////
  //                                               //
  // CALCOLO IL FATTORE DI FORMA f0 DALLO SCALARE  //
  //                                               //
  ///////////////////////////////////////////////////
     
  double *f0_S_ml_to_mH=(double*)malloc(sizeof(double)*clusterfile), *f0_S_mH_to_ml=(double*)malloc(sizeof(double)*clusterfile);
  double *f0_S_Hl_sym=(double*)malloc(sizeof(double)*clusterfile), *f0_S_Hl_sym_av=(double*)malloc(sizeof(double)*clusterfile);
  double *f0_rest_Hl_sym_av=(double*)malloc(sizeof(double)*clusterfile);

  double sigma_f0_S_ml_to_mH, sigma_f0_S_mH_to_ml, sigma_f0_S_Hl_sym, sigma_f0_S_Hl_sym_av, sigma_f0_rest_Hl_sym_av;
  
  for(int ijk = 0; ijk < clusterfile; ijk++){
    
    f0_S_ml_to_mH[ijk]     = ( (mq[ibeta][im2] - mq[ibeta][im1])/( pow(mass_2pts_mH[ijk],2) - pow(mass_2pts_ml[ijk],2) ))*S0_fit_ml_to_mH[ijk];
    f0_S_mH_to_ml[ijk]     = ( (mq[ibeta][im2] - mq[ibeta][im1])/( pow(mass_2pts_mH[ijk],2) - pow(mass_2pts_ml[ijk],2) ))*S0_fit_mH_to_ml[ijk];
    f0_S_Hl_sym[ijk]       = ( (mq[ibeta][im2] - mq[ibeta][im1])/( pow(mass_2pts_mH[ijk],2) - pow(mass_2pts_ml[ijk],2) ))*S0_fit_Hl_sym[ijk];
    f0_S_Hl_sym_av[ijk]    = ( (mq[ibeta][im2] - mq[ibeta][im1])/( pow(mass_2pts_mH[ijk],2) - pow(mass_2pts_ml[ijk],2) ))*S0_fit_Hl_sym_av[ijk];

    f0_rest_Hl_sym_av[ijk] = (S0_fit_Hl_sym_av[ijk]/S0_fit_Hl_rest_sym_av[ijk])*( V0_fit_Hl_rest_sym_av[ijk]/(mass_2pts_mH[ijk] + mass_2pts_ml[ijk] ));
  }
  free(S0_fit_ml_to_mH);
  free(S0_fit_mH_to_ml);
  free(S0_fit_Hl_sym);  
  free(S0_fit_Hl_sym_av);
  free(S0_fit_Hl_rest_sym_av);
  free(V0_fit_Hl_rest_sym_av);
  
  sigma_f0_S_ml_to_mH = sigma_JK(f0_S_ml_to_mH, clusterfile);
  sigma_f0_S_mH_to_ml = sigma_JK(f0_S_mH_to_ml, clusterfile);
  sigma_f0_S_Hl_sym = sigma_JK(f0_S_Hl_sym, clusterfile);
  sigma_f0_S_Hl_sym_av = sigma_JK(f0_S_Hl_sym_av, clusterfile);
  sigma_f0_rest_Hl_sym_av = sigma_JK(f0_rest_Hl_sym_av, clusterfile);

  ////// SCRIVO L'OUTPUT DI f0 DALLO SCALARE
  FILE *fout_fit_f0_S_ml_to_mH, *fout_fit_f0_S_mH_to_ml, *fout_fit_f0_S_Hl_sym, *fout_fit_f0_S_Hl_sym_av;
  
  char *file_out_fit_f0_S_ml_to_mH=(char*)malloc(sizeof(char)*LEN_NAME);
  char *file_out_fit_f0_S_mH_to_ml=(char*)malloc(sizeof(char)*LEN_NAME);
  char *file_out_fit_f0_S_Hl_sym=(char*)malloc(sizeof(char)*LEN_NAME);
  char *file_out_fit_f0_S_Hl_sym_av=(char*)malloc(sizeof(char)*LEN_NAME);
  
  sprintf(file_out_fit_f0_S_ml_to_mH, "OUTPUT_SMEAR/%s/%s/f0_S/light_to_Heavy/%s/%s/f0_S.%s.m1m2_%s%s.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  sprintf(file_out_fit_f0_S_mH_to_ml, "OUTPUT_SMEAR/%s/%s/f0_S/Heavy_to_light/%s/%s/f0_S.%s.m1m2_%s%s.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  sprintf(file_out_fit_f0_S_Hl_sym, "OUTPUT_SMEAR/%s/%s/f0_S/Hl_sym/%s/%s/f0_S.%s.m1m2_%s%s.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  sprintf(file_out_fit_f0_S_Hl_sym_av, "OUTPUT_SMEAR/%s/%s/f0_S/Hl_sym_av/%s/%s/f0_S.%s.m1m2_%s%s.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  
  if ((fout_fit_f0_S_ml_to_mH = fopen(file_out_fit_f0_S_ml_to_mH, "w")) == NULL ){
    printf("Error opening the output file: file_out_fit_f0_S_ml_to_mH\n");
    exit(EXIT_FAILURE);
  }
  free(file_out_fit_f0_S_ml_to_mH);

  if ((fout_fit_f0_S_mH_to_ml = fopen(file_out_fit_f0_S_mH_to_ml, "w")) == NULL ){
    printf("Error opening the output file: file_out_fit_f0_S_mH_to_ml\n");
    exit(EXIT_FAILURE);
  }
  free(file_out_fit_f0_S_mH_to_ml);
  
  if ((fout_fit_f0_S_Hl_sym = fopen(file_out_fit_f0_S_Hl_sym, "w")) == NULL ){
    printf("Error opening the output file: file_out_fit_f0_S_Hl_sym\n");
    exit(EXIT_FAILURE);
  }
  free(file_out_fit_f0_S_Hl_sym);
  
  if ((fout_fit_f0_S_Hl_sym_av = fopen(file_out_fit_f0_S_Hl_sym_av, "w")) == NULL ){
    printf("Error opening the output file: file_out_fit_f0_S_Hl_sym_av\n");
    exit(EXIT_FAILURE);
  }
  free(file_out_fit_f0_S_Hl_sym_av);
  
  for(int ijk = 0; ijk < clusterfile; ijk++){
    
    fprintf(fout_fit_f0_S_ml_to_mH, "%f\n", f0_S_ml_to_mH[ijk]);
    fprintf(fout_fit_f0_S_mH_to_ml, "%f\n", f0_S_mH_to_ml[ijk]);
    fprintf(fout_fit_f0_S_Hl_sym, "%f\n", f0_S_Hl_sym[ijk]);
    fprintf(fout_fit_f0_S_Hl_sym_av, "%f\n", f0_S_Hl_sym_av[ijk]);
  }
  fclose(fout_fit_f0_S_ml_to_mH);
  fclose(fout_fit_f0_S_mH_to_ml);
  fclose(fout_fit_f0_S_Hl_sym);
  fclose(fout_fit_f0_S_Hl_sym_av);
  ///////// FINE CALCOLO IL FATTORE DI FORMA f0 DALLO SCALARE

  ////////////////////////////////
  //                            //
  // CALCOLO I FATTORI DI FORMA //
  //                            //
  ////////////////////////////////

  ////// CASO IN CUI ALMENO UNO DEI DUE MESONI È IN MOTO: FACCIO IL FIT
  double *fzero_fit_ml_to_mH=(double*)malloc(sizeof(double)*clusterfile);
  double *fzero_fit_mH_to_ml=(double*)malloc(sizeof(double)*clusterfile);
  double *fzero_fit_Hl_sym=(double*)malloc(sizeof(double)*clusterfile);
  double *fzero_fit_Hl_sym_av=(double*)malloc(sizeof(double)*clusterfile);

  double *fplus_fit_ml_to_mH=(double*)malloc(sizeof(double)*clusterfile);
  double *fplus_fit_mH_to_ml=(double*)malloc(sizeof(double)*clusterfile);
  double *fplus_fit_Hl_sym=(double*)malloc(sizeof(double)*clusterfile);
  double *fplus_fit_Hl_sym_av=(double*)malloc(sizeof(double)*clusterfile);

  double *fminus_fit_ml_to_mH=(double*)malloc(sizeof(double)*clusterfile);
  double *fminus_fit_mH_to_ml=(double*)malloc(sizeof(double)*clusterfile);
  double *fminus_fit_Hl_sym=(double*)malloc(sizeof(double)*clusterfile);
  double *fminus_fit_Hl_sym_av=(double*)malloc(sizeof(double)*clusterfile);

  double sigma_fzero_fit_ml_to_mH, sigma_fzero_fit_mH_to_ml, sigma_fzero_fit_Hl_sym, sigma_fzero_fit_Hl_sym_av;
  double sigma_fplus_fit_ml_to_mH, sigma_fplus_fit_mH_to_ml, sigma_fplus_fit_Hl_sym, sigma_fplus_fit_Hl_sym_av;
  double sigma_fminus_fit_ml_to_mH, sigma_fminus_fit_mH_to_ml, sigma_fminus_fit_Hl_sym, sigma_fminus_fit_Hl_sym_av;
  
  double outpar[npar], err[npar];
  //double par[npar] = {0.005, 0.005};
  double par[npar]  = {0.60, 0.60};
  double step[npar] = {0.01, 0.01};
  double min[npar]  = {0.00, 0.00};
  double max[npar]  = {0.00, 0.00};
  string cpar[npar] = {"f_zero", "f_plus"};

  if(ith1 != 3 || ith2 != 3){
    
    for(int choice = 0; choice <= 3; choice++){
      for(int ijk = 0; ijk < clusterfile; ijk++){
	
        q2_glb  = q2[ijk];
        E_l_glb = energy_2pts_m1[ijk]; 
        E_H_glb = energy_2pts_m2[ijk];
        pl_glb  = momentum_1; 
        pH_glb  = momentum_2; 
        M_l_glb = mass_2pts_ml[ijk];
        M_H_glb = mass_2pts_mH[ijk];
	
        if(choice == 0){
	  
          V0_fit_glb = V0_fit_ml_to_mH[ijk];
          Vi_fit_glb = Vi_fit_ml_to_mH[ijk];
          f0_S_fit_glb = f0_S_ml_to_mH[ijk];

          sigma_V0_fit_glb = sigma_V0_fit_ml_to_mH;
          sigma_Vi_fit_glb = sigma_Vi_fit_ml_to_mH;
          sigma_f0_S_fit_glb = sigma_f0_S_ml_to_mH*sigma_S_enhancement;
        }// choice 0 ml_to_mH
	
        if(choice == 1){

          V0_fit_glb = V0_fit_mH_to_ml[ijk];
          Vi_fit_glb = Vi_fit_mH_to_ml[ijk];
          f0_S_fit_glb = f0_S_mH_to_ml[ijk];

          sigma_V0_fit_glb = sigma_V0_fit_mH_to_ml;
          sigma_Vi_fit_glb = sigma_Vi_fit_mH_to_ml;
          sigma_f0_S_fit_glb = sigma_f0_S_mH_to_ml*sigma_S_enhancement;
	      }// choice 1 mH_to_ml
	
        if(choice == 2){
	  
          V0_fit_glb = V0_fit_Hl_sym[ijk];
          Vi_fit_glb = Vi_fit_Hl_sym[ijk];
          f0_S_fit_glb = f0_S_Hl_sym[ijk];

          sigma_V0_fit_glb = sigma_V0_fit_Hl_sym;
          sigma_Vi_fit_glb = sigma_Vi_fit_Hl_sym;
          sigma_f0_S_fit_glb = sigma_f0_S_Hl_sym*sigma_S_enhancement;
        }// choice 2 Hl_sym
	
        if(choice == 3){

          V0_fit_glb = V0_fit_Hl_sym_av[ijk];
          Vi_fit_glb = Vi_fit_Hl_sym_av[ijk];
	  
          if(if0 == 0){
            f0_S_fit_glb = f0_S_Hl_sym_av[ijk];
          }else if(if0 == 1){
            f0_S_fit_glb = f0_rest_Hl_sym_av[ijk];
          }

          sigma_V0_fit_glb = sigma_V0_fit_Hl_sym_av;
          sigma_Vi_fit_glb = sigma_Vi_fit_Hl_sym_av;
	  
          if(if0 == 0){
            sigma_f0_S_fit_glb = sigma_f0_S_Hl_sym_av*sigma_S_enhancement;
          }else if(if0 == 1){
            sigma_f0_S_fit_glb = sigma_f0_rest_Hl_sym_av*sigma_S_enhancement;
          } 
        }// choice 3 Hl_sym_av

	TMinuit minuit(npar);

	minuit.SetFCN(chi2);

	minuit.SetErrorDef(1.);

	for(int j = 0; j < npar; j++){
	  minuit.DefineParameter(j, cpar[j].c_str(), par[j], step[j], min[j], max[j]);
	}

	minuit.Migrad();
	
	for(int j = 0; j < npar; j++){
	  minuit.GetParameter(j, outpar[j], err[j]); 
	}
	
	if(choice == 0){
	  fzero_fit_ml_to_mH[ijk] = outpar[0];
	  fplus_fit_ml_to_mH[ijk] = outpar[1];
	  fminus_fit_ml_to_mH[ijk] = (outpar[0]-outpar[1])*( ( pow(M_H_glb ,2)-pow(M_l_glb ,2) )/q2_glb  );		 
	}
	if(choice == 1){
	  fzero_fit_mH_to_ml[ijk] = outpar[0];
	  fplus_fit_mH_to_ml[ijk] = outpar[1];
	  fminus_fit_mH_to_ml[ijk] = (outpar[0]-outpar[1])*( ( pow(M_H_glb ,2)-pow(M_l_glb ,2) )/q2_glb  );
	}
	if(choice == 2){
	  fzero_fit_Hl_sym[ijk] = outpar[0];
	  fplus_fit_Hl_sym[ijk] = outpar[1];
	  fminus_fit_Hl_sym[ijk] = (outpar[0]-outpar[1])*( ( pow(M_H_glb ,2)-pow(M_l_glb ,2) )/q2_glb  );
	}
	if(choice == 3){
	  fzero_fit_Hl_sym_av[ijk] = outpar[0];
	  fplus_fit_Hl_sym_av[ijk] = outpar[1];
	  fminus_fit_Hl_sym_av[ijk] = (outpar[0]-outpar[1])*( ( pow(M_H_glb ,2)-pow(M_l_glb ,2) )/q2_glb  );
	}
	
	if(ijk == 15){
	  
	  chi2_f1 = pow( (V0_fit_glb - V0_function( outpar)), 2)/( pow(sigma_V0_fit_glb, 2)) ;
	  
	  chi2_f2 = pow( (Vi_fit_glb - Vi_function( outpar)), 2)/( pow(sigma_Vi_fit_glb, 2)) ;
	  
	  chi2_f3 = pow( f0_S_fit_glb -  outpar[0], 2)/( pow(sigma_f0_S_fit_glb, 2)) ;
	  
	  chi2_f = chi2_f1 + chi2_f2 + chi2_f3;

	  printf("\n\nChi2\n\n");
	  printf("Chi2 = %f   per la scelta %d\n", chi2_f, choice);
	  
	}// if(ijk == 15)
	
      } // ijk
      
      if(choice == 0){
	sigma_fzero_fit_ml_to_mH = sigma_JK(fzero_fit_ml_to_mH, clusterfile);
	sigma_fplus_fit_ml_to_mH = sigma_JK(fplus_fit_ml_to_mH, clusterfile);
	sigma_fminus_fit_ml_to_mH = sigma_JK(fminus_fit_ml_to_mH, clusterfile);
      }
      if(choice == 1){
	sigma_fzero_fit_mH_to_ml = sigma_JK(fzero_fit_mH_to_ml, clusterfile);
	sigma_fplus_fit_mH_to_ml = sigma_JK(fplus_fit_mH_to_ml, clusterfile);
	sigma_fminus_fit_mH_to_ml = sigma_JK(fminus_fit_mH_to_ml, clusterfile);
      }
      if(choice == 2){
	sigma_fzero_fit_Hl_sym = sigma_JK(fzero_fit_Hl_sym, clusterfile);
	sigma_fplus_fit_Hl_sym = sigma_JK(fplus_fit_Hl_sym, clusterfile);
	sigma_fminus_fit_Hl_sym = sigma_JK(fminus_fit_Hl_sym, clusterfile);
      }
      if(choice == 3){
	sigma_fzero_fit_Hl_sym_av = sigma_JK(fzero_fit_Hl_sym_av, clusterfile);
	sigma_fplus_fit_Hl_sym_av = sigma_JK(fplus_fit_Hl_sym_av, clusterfile);
	sigma_fminus_fit_Hl_sym_av = sigma_JK(fminus_fit_Hl_sym_av, clusterfile);
      }
      
    } // choice

  }// if(ith1 != 3 || ith2 != 3)
  free(Vi_fit_ml_to_mH);
  free(Vi_fit_mH_to_ml);
  free(Vi_fit_Hl_sym);
  free(Vi_fit_Hl_sym_av);
  
  free(energy_2pts_m1);
  free(energy_2pts_m2);

  
  ////// CASO CON ENTRAMBI I MESONI FERMI
  double *fzero_ferma_ml_to_mH=(double*)malloc(sizeof(double)*clusterfile), *fzero_ferma_mH_to_ml=(double*)malloc(sizeof(double)*clusterfile);
  double *fzero_ferma_Hl_sym=(double*)malloc(sizeof(double)*clusterfile), *fzero_ferma_Hl_sym_av=(double*)malloc(sizeof(double)*clusterfile);
  
  double sigma_fzero_ferma_ml_to_mH, sigma_fzero_ferma_mH_to_ml, sigma_fzero_ferma_Hl_sym, sigma_fzero_ferma_Hl_sym_av;

  double num1 = 0, den1 = 0, num2 = 0, den2 = 0, num3 = 0, den3 = 0, num4 = 0, den4 = 0;
  
  if(ith1 == 3 && ith2 == 3){
    
    for(int choice = 0; choice <= 3; choice++){
      for(int ijk = 0; ijk < clusterfile; ijk++){
	
	if(choice == 0){
	  fzero_ferma_ml_to_mH[ijk] = V0_fit_ml_to_mH[ijk]/(mass_2pts_mH[ijk] + mass_2pts_ml[ijk]);
	  fplus_fit_ml_to_mH[ijk] = 0;
	  fminus_fit_ml_to_mH[ijk] = 0;
	}
	if(choice == 1){
	  fzero_ferma_mH_to_ml[ijk] = V0_fit_mH_to_ml[ijk]/(mass_2pts_mH[ijk] + mass_2pts_ml[ijk]);
	  fplus_fit_mH_to_ml[ijk] = 0;
	  fminus_fit_mH_to_ml[ijk] = 0;
	}
	if(choice == 2){
	  fzero_ferma_Hl_sym[ijk] = V0_fit_Hl_sym[ijk]/(mass_2pts_mH[ijk] + mass_2pts_ml[ijk]);
	  fplus_fit_Hl_sym[ijk] = 0;
	  fminus_fit_Hl_sym[ijk] = 0;	 
	}
	if(choice == 3){
	  fzero_ferma_Hl_sym_av[ijk] = V0_fit_Hl_sym_av[ijk]/(mass_2pts_mH[ijk] + mass_2pts_ml[ijk]);
	  fplus_fit_Hl_sym_av[ijk] = 0;
	  fminus_fit_Hl_sym_av[ijk] = 0;	 
	}
      }// ijk
      
      if(choice == 0){
	sigma_fzero_ferma_ml_to_mH = sigma_JK(fzero_ferma_ml_to_mH, clusterfile);
	sigma_fplus_fit_ml_to_mH = sigma_JK(fplus_fit_ml_to_mH, clusterfile);
	sigma_fminus_fit_ml_to_mH = sigma_JK(fminus_fit_ml_to_mH, clusterfile);
      }
      if(choice == 1){
	sigma_fzero_ferma_mH_to_ml = sigma_JK(fzero_ferma_mH_to_ml, clusterfile);
	sigma_fplus_fit_mH_to_ml = sigma_JK(fplus_fit_mH_to_ml, clusterfile);
	sigma_fminus_fit_mH_to_ml = sigma_JK(fminus_fit_mH_to_ml, clusterfile);
      }
      if(choice == 2){
	sigma_fzero_ferma_Hl_sym = sigma_JK(fzero_ferma_Hl_sym, clusterfile);
	sigma_fplus_fit_Hl_sym = sigma_JK(fplus_fit_Hl_sym, clusterfile);
	sigma_fminus_fit_Hl_sym = sigma_JK(fminus_fit_Hl_sym, clusterfile);
      }
      if(choice == 3){
	sigma_fzero_ferma_Hl_sym_av = sigma_JK(fzero_ferma_Hl_sym_av, clusterfile);
	sigma_fplus_fit_Hl_sym_av = sigma_JK(fplus_fit_Hl_sym_av, clusterfile);
	sigma_fminus_fit_Hl_sym_av = sigma_JK(fminus_fit_Hl_sym_av, clusterfile);
      }
      
      //////// CALCOLO f0 TRAMITE FIT A COSTANTE
      for(int ijk = 0; ijk < clusterfile; ijk++){
	
	if(choice == 0){
	  num1 = ( fzero_ferma_ml_to_mH[ijk]*(1/pow(sigma_fzero_ferma_ml_to_mH,2)) ) + ( f0_S_ml_to_mH[ijk]*(1/pow(sigma_f0_S_ml_to_mH,2)) );
	  den1 = ( 1/pow(sigma_fzero_ferma_ml_to_mH,2)  ) + ( 1/pow(sigma_f0_S_ml_to_mH,2)  );
	  fzero_fit_ml_to_mH[ijk] = num1/den1;
	}
	if(choice == 1){
	  num2 = ( fzero_ferma_mH_to_ml[ijk]*(1/pow(sigma_fzero_ferma_mH_to_ml,2)) ) + ( f0_S_mH_to_ml[ijk]*(1/pow(sigma_f0_S_mH_to_ml,2)) );
	  den2 = ( 1/pow(sigma_fzero_ferma_mH_to_ml,2) ) + ( 1/pow(sigma_f0_S_mH_to_ml,2) );
	  fzero_fit_mH_to_ml[ijk] = num2/den2;
	}
	if(choice == 2){
	  num3 = ( fzero_ferma_Hl_sym[ijk]*(1/pow(sigma_fzero_ferma_Hl_sym,2)) ) + ( f0_S_Hl_sym[ijk]*(1/pow(sigma_f0_S_Hl_sym,2)) );
	  den3 = ( 1/pow(sigma_fzero_ferma_Hl_sym,2) ) + ( 1/pow(sigma_f0_S_Hl_sym,2) );
	  fzero_fit_Hl_sym[ijk] = num3/den3;
	}
	if(choice == 3){
	  num4 = ( fzero_ferma_Hl_sym_av[ijk]*(1/pow(sigma_fzero_ferma_Hl_sym_av,2)) ) + ( f0_S_Hl_sym_av[ijk]*(1/pow(sigma_f0_S_Hl_sym_av,2)) );
	  den4 = ( 1/pow(sigma_fzero_ferma_Hl_sym_av,2) ) + ( 1/pow(sigma_f0_S_Hl_sym_av,2) );
	  fzero_fit_Hl_sym_av[ijk] = num4/den4;
	}
	
      }// ijk
      
      if(choice == 0){
	sigma_fzero_fit_ml_to_mH = sigma_JK(fzero_fit_ml_to_mH, clusterfile);
      }
      if(choice == 1){
	sigma_fzero_fit_mH_to_ml = sigma_JK(fzero_fit_mH_to_ml, clusterfile);
      }
      if(choice == 2){
	sigma_fzero_fit_Hl_sym = sigma_JK(fzero_fit_Hl_sym, clusterfile);
      }
      if(choice == 3){
	sigma_fzero_fit_Hl_sym_av = sigma_JK(fzero_fit_Hl_sym_av, clusterfile);
      }
      
    }// choice    
    free(fzero_ferma_ml_to_mH);
    free(fzero_ferma_mH_to_ml);
    free(fzero_ferma_Hl_sym);
    free(fzero_ferma_Hl_sym_av);
    
  }// if(ith1 == 3 && ith2 == 3)

  free(V0_fit_ml_to_mH);
  free(V0_fit_mH_to_ml);
  free(V0_fit_Hl_sym);
  free(V0_fit_Hl_sym_av);
  
  free(f0_S_mH_to_ml);
  free(f0_S_Hl_sym);
  free(f0_S_Hl_sym_av);
  
  free(mass_2pts_ml);
  free(mass_2pts_mH);

  
  /////////////////////////
  //                     //
  // PRINT DI CONTROLLO  //
  //                     //
  /////////////////////////

  printf("Results for: %s %s %s %s %s %s\n", beta_V[ibeta].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  
  printf("f0SlH  %f %f\n", f0_S_ml_to_mH[15], sigma_f0_S_ml_to_mH);
  printf("f0res  %f %f\n", f0_rest_Hl_sym_av[15], sigma_f0_rest_Hl_sym_av);
  printf("check  %f %f %f %f %f %f\n", fzero_check[15], sigma_fzero_check, fplus_check[15], sigma_fplus_check, fminus_check[15], sigma_fminus_check);
  printf("fltoH  %f %f %f %f\n", fzero_fit_ml_to_mH[15], sigma_fzero_fit_ml_to_mH, fplus_fit_ml_to_mH[15], sigma_fplus_fit_ml_to_mH);
  printf("fHtol  %f %f %f %f\n", fzero_fit_mH_to_ml[15], sigma_fzero_fit_mH_to_ml, fplus_fit_mH_to_ml[15], sigma_fplus_fit_mH_to_ml);
  printf("fHlsy  %f %f %f %f\n", fzero_fit_Hl_sym[15], sigma_fzero_fit_Hl_sym, fplus_fit_Hl_sym[15], sigma_fplus_fit_Hl_sym);
  printf("fHlav  %f %f %f %f\n", fzero_fit_Hl_sym_av[15], sigma_fzero_fit_Hl_sym_av, fplus_fit_Hl_sym_av[15], sigma_fplus_fit_Hl_sym_av);
  printf("q2=%f, chi2_Hlav=%f\n", q2[15], chi2_f);
  
  free(fzero_check);
  free(fplus_check);
  free(fminus_check);
  
  free(f0_S_ml_to_mH);
  free(f0_rest_Hl_sym_av);
  

  ////////////
  //        //
  // OUTPUT //
  //        //
  ////////////
  
  ////// FZERO
  write_output_form_factors_fit( fzero_fit_ml_to_mH, q2, dir_S, dir_E, form_factors, strategy, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 0, 0);
  write_output_form_factors_fit( fzero_fit_mH_to_ml, q2, dir_S, dir_E, form_factors, strategy, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 0, 1);
  write_output_form_factors_fit( fzero_fit_Hl_sym, q2, dir_S, dir_E, form_factors, strategy, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 0, 2);
  write_output_form_factors_fit( fzero_fit_Hl_sym_av, q2, dir_S, dir_E, form_factors, strategy, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 0, 3);
  
  if(ith1 !=3 || ith2 !=3 ){
    
    ////// FPLUS
    write_output_form_factors_fit( fplus_fit_ml_to_mH, q2, dir_S, dir_E, form_factors, strategy, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 1, 0);
    write_output_form_factors_fit( fplus_fit_mH_to_ml, q2, dir_S, dir_E, form_factors, strategy, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 1, 1);
    write_output_form_factors_fit( fplus_fit_Hl_sym, q2, dir_S, dir_E, form_factors, strategy, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 1, 2);
    write_output_form_factors_fit( fplus_fit_Hl_sym_av, q2, dir_S, dir_E, form_factors, strategy, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 1, 3);
    
    ////// FMINUS
    write_output_form_factors_fit( fminus_fit_ml_to_mH, q2, dir_S, dir_E, form_factors, strategy, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 2, 0);
    write_output_form_factors_fit( fminus_fit_mH_to_ml, q2, dir_S, dir_E, form_factors, strategy, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 2, 1);
    write_output_form_factors_fit( fminus_fit_Hl_sym, q2, dir_S, dir_E, form_factors, strategy, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 2, 2);
    write_output_form_factors_fit( fminus_fit_Hl_sym_av, q2, dir_S, dir_E, form_factors, strategy, beta_V, mu_sea_1, mu_sea_2, m_, th_, clusterfile, iS, iE, ibeta, imusea, im1, im2, ith1, ith2, 2, 3);
    
  } // if(ith1 !=3 || ith2 !=3) 

  free(q2);
  
  return 0;
  
}


double V0_function( double *par){

  // par[0] = f_zero, par[1] = f_plus;
  return  par[1]*(E_H_glb + E_l_glb) + (par[0]-par[1])*( ( pow(M_H_glb ,2)-pow(M_l_glb ,2) )/q2_glb  )*(E_H_glb - E_l_glb) ;
}

double Vi_function( double *par){

  // par[0] = f_zero, par[1] = f_plus;
  return  par[1]*(pH_glb + pl_glb) + (par[0]-par[1])*( ( pow(M_H_glb ,2)-pow(M_l_glb ,2) )/q2_glb  )*(pH_glb - pl_glb) ;
}


void chi2( int &npar, double *deriv, double &f, double *par, int iflag){
  
  double f1, f2, f3;
  f = 0;

  f1 = pow( (V0_fit_glb - V0_function( par)), 2)/( pow(sigma_V0_fit_glb, 2)) ;

  f2 = pow( (Vi_fit_glb - Vi_function( par)), 2)/( pow(sigma_Vi_fit_glb, 2)) ;

  f3 = pow( f0_S_fit_glb -  par[0], 2)/( pow(sigma_f0_S_fit_glb, 2)) ;

  f = f1 + f2 + f3;
  
}

void sigma_matrix_el( double *sigma_matix_el, double *matix_el, int Tmax, int Nt, int clusterfile){
  
  double *temp=(double*)malloc(sizeof(double)*(clusterfile));
  
  for(int t = 0; t <= Tmax; t++){
    for(int ijk = 0; ijk < clusterfile; ijk++){
      
      temp[ijk] = matix_el[ijk*Nt+t];
      
    }// ijk

    sigma_matix_el[t] = sigma_JK(temp, clusterfile);
    
  }// t

  free(temp);
  
}// sigma_matrix_el

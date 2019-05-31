#include <TMinuit.h>
#include <TF1.h>
#include <fstream>
#include "lettura_ult_inp.h"
#include "stat_analysis_func.h"
#include "definitions.h"

using namespace std;

#define LEN_NAME 1024
#define PI 3.141592653589793

const int Nanalysis = 8;
int analysis_in = 1, analysis_fin = 8;
int start_ev = (analysis_in - 1)*Nev_an + 1;  // Da questi cambi il numero di eventi su cui fare il fit
int end_ev = analysis_fin*Nev_an;             // Da questi cambi il numero di eventi su cui fare il fit

const int npar = 45;

double LambdaQCD = 0.35;
double alpha_eff = 1;

double inv_MDsStar2 = 0.224188;

#if defined(STD) || defined(CUT) || defined(NO_HYP)
int block_par = 1; // polo_f+ fixed
#endif
#if defined(A4) || defined(A4_CUT)
int block_par = 5; // polo_f+ fixed & Prior sui parametri a4
#endif

#if defined(STD)
const int nfix = 30;
double fix_par[nfix] = {4, 7, 10, 13, 15, 16, 17, 18, 19, 21, 23, 24, 25, 26, 27, 28, 30, 31, 32, 33, 34, 35, 36, 37, 39, 40, 41, 42, 43, 44}; // 30
#endif
#if defined(CUT)
const int nfix = 31;
double fix_par[nfix] = {3, 4, 7, 10, 13, 15, 16, 17, 18, 19, 21, 23, 24, 25, 26, 27, 28, 30, 31, 32, 33, 34, 35, 36, 37, 39, 40, 41, 42, 43, 44}; // 31
#endif
#if defined(A4)
const int nfix = 26;
double fix_par[nfix] = {4, 7, 10, 13, 15, 16, 17, 18, 19, 21, 23, 24, 25, 26, 27, 28, 30, 32, 33, 34, 35, 36, 37, 39, 40, 41}; // 26
#endif
#if defined(A4_CUT)
const int nfix = 27;
double fix_par[nfix] = {3, 4, 7, 10, 13, 15, 16, 17, 18, 19, 21, 23, 24, 25, 26, 27, 28, 30, 32, 33, 34, 35, 36, 37, 39, 40, 41}; // 27
#endif
#if defined(NO_HYP)
const int nfix = 35;
double fix_par[nfix] = {4, 7, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 39, 40, 41, 42, 43, 44}; // 35
#endif

int out = 1;
int dati_sint = 1;

/*
                          {"f_0 -> [0]", 
			  "f_ml -> [1]", 
                          "f_a2 -> [2]",  
                          "f_q2 -> [3]",
                          "polo_fzero_0 -> [4]", 
                          "polo_fplus_0 -> [5]",
			  "A_fzero_0 -> [6]",
                          "A_fzero_ml -> [7]",
                          "A_fzero_a2 -> [8]",
			  "A_fplus_0 -> [9]",
                          "A_fplus_ml -> [10]",
                          "A_fplus_a2 -> [11]",  
                          "a_H1_0 -> [12]", 
                          "a_H1_ml -> [13]", 
                          "b_H1_0 -> [14]", 
                          "b_H1_ml -> [15]", 
                          "a_H2_0 -> [16]", 
                          "a_H2_ml -> [17]", 
                          "b_H2_0 -> [18]", 
                          "b_H2_ml -> [19]",
                          "a_H3_0 -> [20]", 
                          "a_H3_ml -> [21]", 
                          "b_H3_0 -> [22]", 
                          "b_H3_ml -> [23]",
                          "a_H4_0 -> [24]", 
                          "a_H4_ml -> [25]", 
                          "b_H4_0 -> [26]", 
                          "b_H4_ml -> [27]"
                          "WI_vi_0 -> [28]",
                          "WI_vi_ml -> [29]",
                          "WI_vi_q2 -> [30]",
                          "f0_S_a2_dis -> [31]",
                          "Az_f0_S_a2_dis -> [32]",
                          "FSE_fplus -> [33]",
                          "FSE_fzero -> [34]",
                          "polo_fzero_ml -> [35]",
                          "polo_fzero_a2 -> [36]", 
                          "polo_fplus_ml -> [37]", 
                          "polo_fplus_a2 -> [38]",
                          "WI_2_vi_0 -> [39]",
                          "WI_2_vi_ml -> [40]",
                          "WI_2_vi_q2 -> [41]",
                          "A_z_term_2_0 -> [42]",
                          "A_z_term_2_ml -> [43]",
                          "A_z_term_2_a2 -> [44]"};
 */

//////////////////////////////

int energia_da_sinhDR = 0;  
int energia_da_stdDR = 1;
int energia_dal_fit = 0;

//////////////////////////////

//////////////////////////////

int f0_with_S = 1;
int f0_without_S = 0;

//////////////////////////////


///////////////////  QUESTE SONO LE GRANDEZZE CHE VANNO DICHIARATE GLOBALI
double PH_fit[Nbeta][Nmusea][Nth1][Nth2], PLi_fit[Nbeta][Nmusea][Nth1][Nth2], a_fit[Nbeta];
double V0_fit[Nbeta][Nmusea][Nth1][Nth2], sigma_V0_fit[Nbeta][Nmusea][Nth1][Nth2];
double Vi_fit[Nbeta][Nmusea][Nth1][Nth2], sigma_Vi_fit[Nbeta][Nmusea][Nth1][Nth2];
double f0_S_fit[Nbeta][Nmusea][Nth1][Nth2], sigma_f0_S_fit[Nbeta][Nmusea][Nth1][Nth2];
double MLi_fit[Nbeta][Nmusea], MH_fit[Nbeta][Nmusea], Mpi_fit[Nbeta][Nmusea];
double ml_fit[Nbeta][Nmusea], Xi_l_fit[Nbeta][Nmusea];

double a_glb, q0, qi, q2, q4, P0, Pi, P4, t0, tp, z, z0, z_term;
double Xil, ml, MH, MLi, EH, ELi, PLi, PH, Mpion;
double V0_glb, Vi_glb, f0_S_glb, sigma_V0_glb, sigma_Vi_glb, sigma_f0_S_glb;
///////////////////

double chi2_num( double *);
void chi2( int &, double *, double &, double *, int);
double fplus_no_hyp_num( double *, double, double, double, double, double, double, int);
double fzero_no_hyp_num( double *, double, double, double, double, double, double, int);

int num_par = 0, num_points = 0;

int main(){
    
  //////////////////////////////

  string dir_E[3] = {"E_from_sinhDR", "E_from_stdDR", "E_from_fit"}, strategy[4] = { "light_to_Heavy", "Heavy_to_light", "Hl_sym", "Hl_sym_av"};
  int iE, istrategy = 3; // questo fissa la cartella Hl_sym_av
  
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


  ////////  LETTURA Ultimate Input
  
  double mlight[Nev+1], mstrange[Nev+1], mcharm[Nev+1], a[Nbeta][Nev+1], ainv[Nbeta][Nev+1], r0[Nev+1], Zev[Nbeta][Nev+1], ZTev[Nbeta][Nev+1], f0[Nev+1], B0[Nev+1], fkfpi[Nev+1];
  int  iboot[Nbeta][Nmusea][Nev+1];
  
  lettura_ultimate_input( mlight, mstrange, mcharm, a, ainv, r0, Zev, iboot, f0, B0, fkfpi, ZTev);
  
  ////////  FINE LETTURA Ultimate Input


  ///////////  Xil_phys
  
  double Xil_phys[Nev+1];
  
  for(int iev = 0; iev <= Nev; iev++){
    
    Xil_phys[iev] = (B0[iev]*mlight[iev])/pow( 4*PI*f0[iev] ,2);
    
  } // iev
  
  ///////////  FINE Xil_phys
  
  
  /////////// CHIAMO V0, Vi, f0_S
  
  FILE *fr_V0, *fr_Vi, *fr_f0_S;
  char open_V0[LEN_NAME], open_Vi[LEN_NAME], open_f0_S[LEN_NAME];
  
  double *****V0=(double*****)malloc(sizeof(double)*(Nbeta));
  for(int ibeta = 0; ibeta < Nbeta; ibeta++) V0[ibeta]=(double****)malloc(sizeof(double)*Nmusea);
  for(int ibeta = 0; ibeta < Nbeta; ibeta++) for(int imusea = 0; imusea < Nmusea; imusea++) V0[ibeta][imusea]=(double***)malloc(sizeof(double)*Nth1);
  for(int ibeta = 0; ibeta < Nbeta; ibeta++) for(int imusea = 0; imusea < Nmusea; imusea++) for(int ith1 = 0; ith1 < Nth1; ith1++) V0[ibeta][imusea][ith1]=(double**)malloc(sizeof(double)*Nth2);
  for(int ibeta = 0; ibeta < Nbeta; ibeta++) for(int imusea = 0; imusea < Nmusea; imusea++) for(int ith1 = 0; ith1 < Nth1; ith1++) for(int ith2 = 0; ith2 < Nth2; ith2++) V0[ibeta][imusea][ith1][ith2]=(double*)malloc(sizeof(double)*(Nev+1));
  
  double *****sigma_V0=(double*****)malloc(sizeof(double)*(Nbeta));
  for(int ibeta = 0; ibeta < Nbeta; ibeta++) sigma_V0[ibeta]=(double****)malloc(sizeof(double)*Nmusea);
  for(int ibeta = 0; ibeta < Nbeta; ibeta++) for(int imusea = 0; imusea < Nmusea; imusea++) sigma_V0[ibeta][imusea]=(double***)malloc(sizeof(double)*Nth1);
  for(int ibeta = 0; ibeta < Nbeta; ibeta++) for(int imusea = 0; imusea < Nmusea; imusea++) for(int ith1 = 0; ith1 < Nth1; ith1++) sigma_V0[ibeta][imusea][ith1]=(double**)malloc(sizeof(double)*Nth2);
  for(int ibeta = 0; ibeta < Nbeta; ibeta++) for(int imusea = 0; imusea < Nmusea; imusea++) for(int ith1 = 0; ith1 < Nth1; ith1++) for(int ith2 = 0; ith2 < Nth2; ith2++) sigma_V0[ibeta][imusea][ith1][ith2]=(double*)malloc(sizeof(double)*(Nanalysis));
  
  double *temp_V0=(double*)malloc(sizeof(double)*(Nev+1));
  
  for(int ibeta = 0; ibeta < Nbeta; ibeta++){
    for(int imusea = 0; imusea < NKL[ibeta]; imusea++){
      
      for(int ith1 = 0; ith1 < Nth1; ith1++){
	for(int ith2 = 0; ith2 < Nth2; ith2++){
	  
	  sprintf(open_V0, "OUTPUT_SMEAR/%s/%s/V0/%s/%s/%s/V0.%s.m1m2_sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
	  
	  if ((fr_V0 = fopen(open_V0, "r")) == NULL ){
	    printf("Error opening the input file: %s\n",open_V0);
	    exit(EXIT_FAILURE);
	  }	  
	  
	  for(int iev = 0; iev <= Nev; iev++){
	    
	    fscanf(fr_V0,"%lf\n", &V0[ibeta][imusea][ith1][ith2][iev]);
	    V0[ibeta][imusea][ith1][ith2][iev] = V0[ibeta][imusea][ith1][ith2][iev]/a[ibeta][iev];
	    temp_V0[iev] = V0[ibeta][imusea][ith1][ith2][iev];	    
	  }
	  
	  fclose(fr_V0);
	  
	  for(int ianalysis = 0; ianalysis < Nanalysis; ianalysis++){
	    
	    sigma_V0[ibeta][imusea][ith1][ith2][ianalysis] = sigma_JK_modified_2(temp_V0, ianalysis*Nev_an + 1, ianalysis*Nev_an + Nev_an, clusterfile);
	    
	  }// ianalysis
	  
	  
	} // ith1
      } // ith2
      
    } // imusea
  } // ibeta
  
  free(V0);
  free(temp_V0);
  
  
  double *****Vi=(double*****)malloc(sizeof(double)*(Nbeta));
  for(int ibeta = 0; ibeta < Nbeta; ibeta++) Vi[ibeta]=(double****)malloc(sizeof(double)*Nmusea);
  for(int ibeta = 0; ibeta < Nbeta; ibeta++) for(int imusea = 0; imusea < Nmusea; imusea++) Vi[ibeta][imusea]=(double***)malloc(sizeof(double)*Nth1);
  for(int ibeta = 0; ibeta < Nbeta; ibeta++) for(int imusea = 0; imusea < Nmusea; imusea++) for(int ith1 = 0; ith1 < Nth1; ith1++) Vi[ibeta][imusea][ith1]=(double**)malloc(sizeof(double)*Nth2);
  for(int ibeta = 0; ibeta < Nbeta; ibeta++) for(int imusea = 0; imusea < Nmusea; imusea++) for(int ith1 = 0; ith1 < Nth1; ith1++) for(int ith2 = 0; ith2 < Nth2; ith2++) Vi[ibeta][imusea][ith1][ith2]=(double*)malloc(sizeof(double)*(Nev+1));
  
  double *****sigma_Vi=(double*****)malloc(sizeof(double)*(Nbeta));
  for(int ibeta = 0; ibeta < Nbeta; ibeta++) sigma_Vi[ibeta]=(double****)malloc(sizeof(double)*Nmusea);
  for(int ibeta = 0; ibeta < Nbeta; ibeta++) for(int imusea = 0; imusea < Nmusea; imusea++) sigma_Vi[ibeta][imusea]=(double***)malloc(sizeof(double)*Nth1);
  for(int ibeta = 0; ibeta < Nbeta; ibeta++) for(int imusea = 0; imusea < Nmusea; imusea++) for(int ith1 = 0; ith1 < Nth1; ith1++) sigma_Vi[ibeta][imusea][ith1]=(double**)malloc(sizeof(double)*Nth2);
  for(int ibeta = 0; ibeta < Nbeta; ibeta++) for(int imusea = 0; imusea < Nmusea; imusea++) for(int ith1 = 0; ith1 < Nth1; ith1++) for(int ith2 = 0; ith2 < Nth2; ith2++) sigma_Vi[ibeta][imusea][ith1][ith2]=(double*)malloc(sizeof(double)*(Nanalysis));
  
  double *temp_Vi=(double*)malloc(sizeof(double)*(Nev+1));
  
  
  for(int ibeta = 0; ibeta < Nbeta; ibeta++){
    for(int imusea = 0; imusea < NKL[ibeta]; imusea++){
      
      for(int ith1 = 0; ith1 < Nth1; ith1++){
	for(int ith2 = 0; ith2 < Nth2; ith2++){
	  
	  sprintf(open_Vi, "OUTPUT_SMEAR/%s/%s/Vi/%s/%s/%s/Vi.%s.m1m2_sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
	  
	  if ((fr_Vi = fopen(open_Vi, "r")) == NULL ){
	    printf("Error opening the input file: %s\n",open_Vi);
	    exit(EXIT_FAILURE);
	  }	  
	  
	  for(int iev = 0; iev <= Nev; iev++){
	    
	    fscanf(fr_Vi,"%lf\n", &Vi[ibeta][imusea][ith1][ith2][iev]);
	    Vi[ibeta][imusea][ith1][ith2][iev] = Vi[ibeta][imusea][ith1][ith2][iev]/a[ibeta][iev];
	    temp_Vi[iev] = Vi[ibeta][imusea][ith1][ith2][iev];	    
	  }
	  
	  fclose(fr_Vi);
	  
	  for(int ianalysis = 0; ianalysis < Nanalysis; ianalysis++){
	    
	    sigma_Vi[ibeta][imusea][ith1][ith2][ianalysis] = sigma_JK_modified_2(temp_Vi, ianalysis*Nev_an + 1, ianalysis*Nev_an + Nev_an, clusterfile);
	    
	  }// ianalysis
	  
	  
	} // ith1
      } // ith2
      
    } // imusea
  } // ibeta
  
  free(Vi);
  free(temp_Vi);
  
  double *****f0_S=(double*****)malloc(sizeof(double)*(Nbeta));
  for(int ibeta = 0; ibeta < Nbeta; ibeta++) f0_S[ibeta]=(double****)malloc(sizeof(double)*Nmusea);
  for(int ibeta = 0; ibeta < Nbeta; ibeta++) for(int imusea = 0; imusea < Nmusea; imusea++) f0_S[ibeta][imusea]=(double***)malloc(sizeof(double)*Nth1);
  for(int ibeta = 0; ibeta < Nbeta; ibeta++) for(int imusea = 0; imusea < Nmusea; imusea++) for(int ith1 = 0; ith1 < Nth1; ith1++) f0_S[ibeta][imusea][ith1]=(double**)malloc(sizeof(double)*Nth2);
  for(int ibeta = 0; ibeta < Nbeta; ibeta++) for(int imusea = 0; imusea < Nmusea; imusea++) for(int ith1 = 0; ith1 < Nth1; ith1++) for(int ith2 = 0; ith2 < Nth2; ith2++) f0_S[ibeta][imusea][ith1][ith2]=(double*)malloc(sizeof(double)*(Nev+1));
  
  double *****sigma_f0_S=(double*****)malloc(sizeof(double)*(Nbeta));
  for(int ibeta = 0; ibeta < Nbeta; ibeta++) sigma_f0_S[ibeta]=(double****)malloc(sizeof(double)*Nmusea);
  for(int ibeta = 0; ibeta < Nbeta; ibeta++) for(int imusea = 0; imusea < Nmusea; imusea++) sigma_f0_S[ibeta][imusea]=(double***)malloc(sizeof(double)*Nth1);
  for(int ibeta = 0; ibeta < Nbeta; ibeta++) for(int imusea = 0; imusea < Nmusea; imusea++) for(int ith1 = 0; ith1 < Nth1; ith1++) sigma_f0_S[ibeta][imusea][ith1]=(double**)malloc(sizeof(double)*Nth2);
  for(int ibeta = 0; ibeta < Nbeta; ibeta++) for(int imusea = 0; imusea < Nmusea; imusea++) for(int ith1 = 0; ith1 < Nth1; ith1++) for(int ith2 = 0; ith2 < Nth2; ith2++) sigma_f0_S[ibeta][imusea][ith1][ith2]=(double*)malloc(sizeof(double)*(Nanalysis));
  
  double *temp_f0_S=(double*)malloc(sizeof(double)*(Nev+1));
  
  for(int ibeta = 0; ibeta < Nbeta; ibeta++){
    for(int imusea = 0; imusea < NKL[ibeta]; imusea++){
      
      for(int ith1 = 0; ith1 < Nth1; ith1++){
	for(int ith2 = 0; ith2 < Nth2; ith2++){
	  
	  sprintf(open_f0_S, "OUTPUT_SMEAR/%s/%s/f0_S/%s/%s/%s/f0_S.%s.m1m2_sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
	  
	  if ((fr_f0_S = fopen(open_f0_S, "r")) == NULL ){
	    printf("Error opening the input file: %s\n",open_f0_S);
	    exit(EXIT_FAILURE);
	  }	  
	  
	  for(int iev = 0; iev <= Nev; iev++){
	    
	    fscanf(fr_f0_S,"%lf\n", &f0_S[ibeta][imusea][ith1][ith2][iev]);
	    temp_f0_S[iev] = f0_S[ibeta][imusea][ith1][ith2][iev];	    
	  }
	  
	  fclose(fr_f0_S);

	  for(int ianalysis = 0; ianalysis < Nanalysis; ianalysis++){
	    
	    sigma_f0_S[ibeta][imusea][ith1][ith2][ianalysis] = sigma_JK_modified_2(temp_f0_S, ianalysis*Nev_an + 1, ianalysis*Nev_an + Nev_an, clusterfile);
	    
	  }// ianalysis
	  
	  
	} // ith1
      } // ith2
      
    } // imusea
  } // ibeta

  free(f0_S);
  free(temp_f0_S);
  
  ////////////////// SIGMA MEDIO V0, Vi, f0_S
 
  double sigma_V0_av[Nbeta][Nmusea][Nth1][Nth2], sigma_Vi_av[Nbeta][Nmusea][Nth1][Nth2], sigma_f0_S_av[Nbeta][Nmusea][Nth1][Nth2];

  for(int ibeta = 0; ibeta < Nbeta; ibeta++){
    for(int imusea = 0; imusea < NKL[ibeta]; imusea++){
      for(int ith1 = 0; ith1 < Nth1; ith1++){
	for(int ith2 = 0; ith2 < Nth2; ith2++){
	  
	  sigma_V0_av[ibeta][imusea][ith1][ith2] = 0; 
	  sigma_Vi_av[ibeta][imusea][ith1][ith2] = 0; 
	  sigma_f0_S_av[ibeta][imusea][ith1][ith2] = 0;
	  
	}
      }
    }
  }

  for(int ibeta = 0; ibeta < Nbeta; ibeta++){
    for(int imusea = 0; imusea < NKL[ibeta]; imusea++){
      for(int ith1 = 0; ith1 < Nth1; ith1++){
	for(int ith2 = 0; ith2 < Nth2; ith2++){
	  for(int ianalysis = 0; ianalysis < Nanalysis; ianalysis++){
	    
	    sigma_V0_av[ibeta][imusea][ith1][ith2] = sigma_V0_av[ibeta][imusea][ith1][ith2] + sigma_V0[ibeta][imusea][ith1][ith2][ianalysis]/Nanalysis;
	    sigma_Vi_av[ibeta][imusea][ith1][ith2] = sigma_Vi_av[ibeta][imusea][ith1][ith2] + sigma_Vi[ibeta][imusea][ith1][ith2][ianalysis]/Nanalysis;
	    sigma_f0_S_av[ibeta][imusea][ith1][ith2] = sigma_f0_S_av[ibeta][imusea][ith1][ith2] + sigma_f0_S[ibeta][imusea][ith1][ith2][ianalysis]/Nanalysis;
	    
	  }
	}
      }
    }
  }
  
  //////////////////// FINE SIGMA MEDIO V0, Vi, f0_S
     
  /////////// FINE CHIAMO V0, Vi, f0_S  



  
  ////////////////////////
  //                    //
  //       fit          //
  //                    //
  ////////////////////////
    
  double parameter[npar][Nev+1];

  double epsilon = 0.000001;
  double chi2_par[npar], chi2_f;

  double outpar[npar], err[npar];
  
  double step[npar] = {0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
  
  double par[npar] = {0.60, 0.00, 0.00, 0.00, 0.03, 0.20, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};
  
  double min[npar] = {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};
  
  double max[npar] = {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};
  
  string cpar[npar] = {"f_0", "f_ml", "f_a2", "f_q2",
		       "polo_fzero_0",
		       "polo_fplus_0",
		       "A_fzero_0", "A_fzero_ml", "A_fzero_a2",
		       "A_fplus_0", "A_fplus_ml", "A_fplus_a2",
		       "a_H1_0", "a_H1_ml", "b_H1_0", "b_H1_ml",
		       "a_H2_0", "a_H2_ml", "b_H2_0", "b_H2_ml",
		       "a_H3_0", "a_H3_ml", "b_H3_0", "b_H3_ml",
		       "a_H4_0", "a_H4_ml", "b_H4_0", "b_H4_ml",
		       "WI_vi_0", "WI_vi_ml", "WI_vi_q2",
		       "f0_S_a2_dis", "Az_f0_S_a2_dis",
		       "FSE_fplus", "FSE_fzero",
		       "polo_fzero_ml", "polo_fzero_a2",
		       "polo_fplus_ml", "polo_fplus_a2",
		       "WI_2_vi_0", "WI_2_vi_ml", "WI_2_vi_q2",
		       "A_z_term_2_0", "A_z_term_2_ml", "A_z_term_2_a2"};
  
  int ian;

  ///////////// QUANTITÀ FIT
  FILE *fr_M_from_2pts_1, *fr_M_from_2pts_2, *fr_Mpi_from_2pts;
  char open_M_from_2pts_fit_1[LEN_NAME], open_M_from_2pts_fit_2[LEN_NAME], open_Mpi_from_2pts[LEN_NAME];
  int Nchar_row = 10;
  ///////////// FINE QUANTITÀ FIT  

  
  for(int iev = 0; iev <= Nev; iev++){
    
    if( (iev == 0) || ( iev >= start_ev && iev <= end_ev ) ){
      
      if(iev == 0){
	ian = 0;
      } else if(iev != 0){
	ian = ((iev-1)/Nev_an);
      }
      
      // TUTTE LE QUANTITÀ USATE NEL FIT SONO IN GEV
      
      for(int ibeta = 0; ibeta < Nbeta; ibeta++){
	for(int imusea = 0; imusea < NKL[ibeta]; imusea++){
	  
	  /////////// CHIAMO LE MASSE
	  sprintf(open_Mpi_from_2pts, "OUTPUT_SMEAR/%s/%s/Mpi/%s/%s/Mpi.%s.m1m2_00.th_3.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str());       
	  sprintf(open_M_from_2pts_fit_1, "OUTPUT_SMEAR/%s/%s/MK/%s/%s/MK.%s.m1m2_0s.th_3.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str());
	  sprintf(open_M_from_2pts_fit_2, "OUTPUT_SMEAR/%s/%s/MD/%s/%s/MD.%s.m1m2_0c.th_3.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str());
	  
	  if ((fr_Mpi_from_2pts = fopen(open_Mpi_from_2pts, "r")) == NULL ){
	    printf("Error opening the file to read: %s\n",fr_Mpi_from_2pts);
	    exit(EXIT_FAILURE);
	  }
	  if ((fr_M_from_2pts_1 = fopen(open_M_from_2pts_fit_1, "r")) == NULL ){
	    printf("Error opening the file to read: %s\n",open_M_from_2pts_fit_1);
	    exit(EXIT_FAILURE);
	  }
	  if ((fr_M_from_2pts_2 = fopen(open_M_from_2pts_fit_2, "r")) == NULL ){
	    printf("Error opening the file to read: %s\n",open_M_from_2pts_fit_2);
	    exit(EXIT_FAILURE);
	  }
	  
	  fseek(fr_M_from_2pts_1, (SEEK_SET + Nchar_row*iev), SEEK_SET);
	  fseek(fr_M_from_2pts_2, (SEEK_SET + Nchar_row*iev), SEEK_SET);
	  fseek(fr_Mpi_from_2pts, (SEEK_SET + Nchar_row*iev), SEEK_SET);
	  
	  fscanf(fr_M_from_2pts_1,"%lf\n", &MLi_fit[ibeta][imusea]);
	  fscanf(fr_M_from_2pts_2,"%lf\n", &MH_fit[ibeta][imusea]);
	  fscanf(fr_Mpi_from_2pts,"%lf\n", &Mpi_fit[ibeta][imusea]);
	  
	  MLi_fit[ibeta][imusea] = MLi_fit[ibeta][imusea]/a[ibeta][iev];
	  MH_fit[ibeta][imusea] = MH_fit[ibeta][imusea]/a[ibeta][iev];
	  Mpi_fit[ibeta][imusea] = Mpi_fit[ibeta][imusea]/a[ibeta][iev];
	  /////////// FINE CHIAMO LE MASSE
	  
	  //if( (ibeta == 0) && (imusea == 0)){
	  //  printf("iev=%d\tMpi=%f\tMpi=%f\tMD=%f\n", iev, Mpi_fit[ibeta][imusea], MLi_fit[ibeta][imusea], MH_fit[ibeta][imusea]);
	  //}
	  
	  for(int ith1 = 0; ith1 < Nth1; ith1++){
	    for(int ith2 = 0; ith2 < Nth2; ith2++){
	      
	      /////////// QUANTITÀ CINEMATICHE
	      PLi_fit[ibeta][imusea][ith1][ith2] = ((theta_value[ibeta][ith1]*PI)/L[ibeta])/a[ibeta][iev];
	      PH_fit[ibeta][imusea][ith1][ith2] = ((theta_value[ibeta][ith2]*PI)/L[ibeta])/a[ibeta][iev];
	      
	      
	      ///////////  MASSA LEGGERO
	      ml_fit[ibeta][imusea] = mq_l[ibeta][imusea]/(Zev[ibeta][iev]*a[ibeta][iev]);
	      
	      
	      ///////////  Xi_l
	      Xi_l_fit[ibeta][imusea] = (B0[iev]*ml_fit[ibeta][imusea])/pow( 4*PI*f0[iev] ,2);	      
	      
	      
	      /////////// LATTICE SPACING
	      a_fit[ibeta] = a[ibeta][iev];

	      
              ////////// MATRIX ELEMENTS
	      sprintf(open_V0, "OUTPUT_SMEAR/%s/%s/V0/%s/%s/%s/V0.%s.m1m2_sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
	      sprintf(open_Vi, "OUTPUT_SMEAR/%s/%s/Vi/%s/%s/%s/Vi.%s.m1m2_sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
	      sprintf(open_f0_S, "OUTPUT_SMEAR/%s/%s/f0_S/%s/%s/%s/f0_S.%s.m1m2_sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
	      
	      if ((fr_V0 = fopen(open_V0, "r")) == NULL ){
		printf("Error opening the input file: %s\n",open_V0);
		exit(EXIT_FAILURE);
	      }
	      if ((fr_Vi = fopen(open_Vi, "r")) == NULL ){
		printf("Error opening the input file: %s\n",open_Vi);
		exit(EXIT_FAILURE);
	      }
	      if ((fr_f0_S = fopen(open_f0_S, "r")) == NULL ){
		printf("Error opening the input file: %s\n",open_f0_S);
		exit(EXIT_FAILURE);
	      }
	      
	      fseek(fr_V0, (SEEK_SET + Nchar_row*iev), SEEK_SET);
	      fseek(fr_Vi, (SEEK_SET + Nchar_row*iev), SEEK_SET);
	      fseek(fr_f0_S, (SEEK_SET + Nchar_row*iev), SEEK_SET);
	      
	      fscanf(fr_V0,"%lf\n", &V0_fit[ibeta][imusea][ith1][ith2]);
	      fscanf(fr_Vi,"%lf\n", &Vi_fit[ibeta][imusea][ith1][ith2]);
	      fscanf(fr_f0_S,"%lf\n", &f0_S_fit[ibeta][imusea][ith1][ith2]);

	      //if( (ibeta == 0) && (imusea == 0) && (ith1 == 0) && (ith2 == 0) ){
	      //  printf("VCVC\tiev=%d\tV0=%f\tVi=%f\tf0_S=%f\n", iev, V0_fit[ibeta][imusea][ith1][ith2], Vi_fit[ibeta][imusea][ith1][ith2], f0_S_fit[ibeta][imusea][ith1][ith2]);
	      //}
	      
	      V0_fit[ibeta][imusea][ith1][ith2] = V0_fit[ibeta][imusea][ith1][ith2]/a[ibeta][iev];
	      Vi_fit[ibeta][imusea][ith1][ith2] = Vi_fit[ibeta][imusea][ith1][ith2]/a[ibeta][iev];
	      
	      sigma_V0_fit[ibeta][imusea][ith1][ith2] = sigma_V0[ibeta][imusea][ith1][ith2][ian];
	      sigma_Vi_fit[ibeta][imusea][ith1][ith2] = sigma_Vi[ibeta][imusea][ith1][ith2][ian];
	      sigma_f0_S_fit[ibeta][imusea][ith1][ith2] = sigma_f0_S[ibeta][imusea][ith1][ith2][ian];
	      
	      if( ith1==3 && ith2 ==3 ){
		Vi_fit[ibeta][imusea][ith1][ith2] = 0;
		sigma_Vi_fit[ibeta][imusea][ith1][ith2] = 100000;
	      }
	      
	      if(iev == 0){
		
		sigma_V0_fit[ibeta][imusea][ith1][ith2] = sigma_V0_av[ibeta][imusea][ith1][ith2];
		sigma_Vi_fit[ibeta][imusea][ith1][ith2] = sigma_Vi_av[ibeta][imusea][ith1][ith2];
		sigma_f0_S_fit[ibeta][imusea][ith1][ith2] = sigma_f0_S_av[ibeta][imusea][ith1][ith2];

		if( ith1==3 && ith2 ==3 ){
		  sigma_Vi_fit[ibeta][imusea][ith1][ith2] = 100000;
		}
	      }
	      
	      fclose(fr_V0);
	      fclose(fr_Vi);
	      fclose(fr_f0_S);
	      ////////// FINE MATRIX ELEMENTS
	      
	    } // ith2
	  } // ith1
	  
	  ////// CHIUSURA FILE MASSE
	  fclose(fr_Mpi_from_2pts);
	  fclose(fr_M_from_2pts_1);
	  fclose(fr_M_from_2pts_2);

	  
	} //imusea
      } // ibeta
      
      TMinuit minuit(npar);
	
      minuit.SetFCN(chi2);
      
      minuit.SetErrorDef(1.);
      
      for(int j = 0; j < npar; j++){
	minuit.DefineParameter(j, cpar[j].c_str(), par[j], step[j], min[j], max[j]);
      }
      
      for(int ipar = 0; ipar < npar; ipar++){
	for(int j = 0; j < nfix; j++){
	  
	  if( ipar == fix_par[j] ){
	    
	    minuit.DefineParameter(ipar, cpar[ipar].c_str(), 0, step[j], min[j], max[j]);
	    
	    minuit.FixParameter(ipar);
	    
	  }
	}
      }
      
      minuit.Migrad();
      
      for(int j = 0; j < npar; j++){
	minuit.GetParameter(j, outpar[j], err[j]);
      }
      
      for(int ipar = 0; ipar < npar; ipar++){
	
	parameter[ipar][iev] = outpar[ipar];
	
	printf("LLL%dL%d\t iev=%d\t %s=%f\n", ipar, ian, iev, cpar[ipar].c_str(), parameter[ipar][iev]);
	
      }
      
      for(int ipar = 0; ipar < npar; ipar++){
	
	if( fabs(outpar[ipar]) > epsilon){
	  num_par = num_par + 1;
	}
	
	chi2_par[ipar] = parameter[ipar][iev];
      }
      
      chi2_f = chi2_num(chi2_par);
      
      printf("CCCC chi2_ev[%d] = %f\t num_par-block_par = %d\t num_points = %d\n", iev, chi2_f, num_par-block_par, num_points);
      
      num_par = 0;
      num_points = 0;
      
    } // if( (iev == 0) || ( iev >= start_ev && iev <= end_ev ) )
    
  } // iev

  
  ////////// STUDIO I PARAMETRI DAL FIT  
  
  double parameter_av[npar][Nanalysis];
  double sigma_parameter[npar];
  
  for(int ipar = 0; ipar < npar; ipar++){
    for( ian = 0; ian < Nanalysis; ian++){
      
      parameter_av[ipar][ian] = 0;
      
    } // ian
  } // ipar
  
  
  for(int ipar = 0; ipar < npar; ipar++){
    for(int iev = 0; iev <= Nev; iev++){
      
      if( iev >= start_ev && iev <= end_ev ){
	
	if(iev == 0){
	  ian = 0;
	} else if(iev != 0){
	  ian = ((iev-1)/Nev_an);
	}
	
	parameter_av[ipar][ian] = parameter_av[ipar][ian] + parameter[ipar][iev]/Nev_an;
	
      } // if( iev >= start_ev && iev <= end_ev )
      
    } // iev
    
    for(ian = 0; ian < Nanalysis; ian++){
      printf("BBBB ian=%d %s_av=%f\n", ian, cpar[ipar].c_str(), parameter_av[ipar][ian]);
    }
    
  } // ipar
  
  for(int ipar = 0; ipar < npar; ipar++){
    sigma_parameter[ipar] = sigma_bootstrap( parameter[ipar], analysis_in, analysis_fin, Nev_an, clusterfile); 
  }
  
  for(int ipar = 0; ipar < npar; ipar++){
    printf("PPPP  %s\t=\t %f\t +- %f\n", cpar[ipar].c_str(), parameter[ipar][0], sigma_parameter[ipar]);
  }

  ////////// STUDIO I PARAMETRI DAL FIT
  
  ////////////////// FINE FIT


  ////////////////// TEST FATTORI DI FORMA SPERIMENTALI
  
  double q2_start = 0, q2_step;
  
  double par_plot[npar];
  
  const int Nstep = 50;
  
  double fzero_plot[Nstep][Nev+1], fplus_plot[Nstep][Nev+1];

  double fzero_plot_sup[Nstep], fplus_plot_sup[Nstep], fzero_plot_inf[Nstep], fplus_plot_inf[Nstep];
  
  double MD_phys = 1.867, MK_phys = 0.4942, Mpi_phys = 0.135;

  double q2_max = pow((MD_phys - MK_phys), 2);
  
  printf("TTTT\n");
  
  for(int iq = 0; iq < Nstep; iq++){
    
    q2_step = q2_start + (q2_max/( (double)  (Nstep - 1)))*iq;
    
    
    for(int iev = 0; iev <= Nev; iev++){
      
      if( (iev == 0) || ( iev >= start_ev && iev <= end_ev ) ){
	
	for(int ipar = 0; ipar < npar; ipar++){
	  
	  par_plot[ipar] = parameter[ipar][iev];
	  
	}
	
	fplus_plot[iq][iev] = fplus_no_hyp_num( par_plot, q2_step, MD_phys, MK_phys, Mpi_phys, Xil_phys[iev], 0.0, 1000);
	
	fzero_plot[iq][iev] = fzero_no_hyp_num( par_plot, q2_step, MD_phys, MK_phys, Mpi_phys, Xil_phys[iev], 0.0, 1000);
	
      } // if
      
    } // iev
    
    
    /////////////////// SIGMA DEI FATTORI DI FORMA PER IL PLOT
    
    double sigma_fzero_plot[Nstep], sigma_fplus_plot[Nstep];
    
    sigma_fzero_plot[iq] = sigma_bootstrap( fzero_plot[iq], analysis_in, analysis_fin, Nev_an, clusterfile);
    
    sigma_fplus_plot[iq] = sigma_bootstrap( fplus_plot[iq], analysis_in, analysis_fin, Nev_an, clusterfile);
    
    /////////////////////////// FINE SIGMA DEI FATTORI DI FORMA PER IL PLOT

    
    fplus_plot_sup[iq] = fplus_plot[iq][0] + sigma_fplus_plot[iq];
    
    fplus_plot_inf[iq] = fplus_plot[iq][0] - sigma_fplus_plot[iq];
    
    fzero_plot_sup[iq] = fzero_plot[iq][0] + sigma_fzero_plot[iq];
    
    fzero_plot_inf[iq] = fzero_plot[iq][0] - sigma_fzero_plot[iq];
    
  } // iq
  
  
  printf("fplus\n");
  
  printf("@type xy\n");
  
  for(int iq = 0; iq < Nstep; iq++){
    
    q2_step = q2_start + (q2_max/( (double)  (Nstep - 1)))*iq;
    
    printf("%f %f\n", q2_step, fplus_plot_sup[iq]);
       
  }
  
  for(int iq = Nstep-1; iq >= 0; iq--){
    
    q2_step = q2_start + (q2_max/( (double)  (Nstep - 1)))*iq;
    
    printf("%f %f\n", q2_step, fplus_plot_inf[iq]);
    
  }
  
  printf("fzero\n");
  
  printf("@type xy\n");
  
  for(int iq = 0; iq < Nstep; iq++){
    
    q2_step = q2_start + (q2_max/( (double)  (Nstep - 1)))*iq;
    
    printf("%f %f\n", q2_step, fzero_plot_sup[iq]);
    
  }
  
  for(int iq = Nstep-1; iq >= 0; iq--){
    
    q2_step = q2_start + (q2_max/( (double)  (Nstep - 1)))*iq;
    
    printf("%f %f\n", q2_step, fzero_plot_inf[iq]);
    
  }
  
  ////////////////// FINE TEST FATTORI DI FORMA SPERIMENTALI


  
  
  ////////////////// CREAZIONE DATI SINTETICI
  
  FILE *fout_sint_f0_q2, *fout_sint_fp_q2;
  
  char file_out_sint_f0_q2[LEN_NAME], file_out_sint_fp_q2[LEN_NAME];
  
  const int Nstep_sint = 8;     
  
  double fzero_sint[Nstep_sint][Nev+1], fplus_sint[Nstep_sint][Nev+1];
  
  if(dati_sint == 1){
    
    for(int iq = 0; iq < Nstep_sint; iq++){
      
      q2_step = q2_start + (q2_max/( (double)  (Nstep_sint-1) ))*iq;
      
      sprintf(file_out_sint_f0_q2, "OUTPUT_SMEAR/%s/%s/dati_sintetici/fzero/fzero_sint_q2_%d.out", dir_S[iS].c_str(), dir_E[iE].c_str(), iq);
      sprintf(file_out_sint_fp_q2, "OUTPUT_SMEAR/%s/%s/dati_sintetici/fplus/fplus_sint_q2_%d.out", dir_S[iS].c_str(), dir_E[iE].c_str(), iq);
      
      if ((fout_sint_f0_q2 = fopen( file_out_sint_f0_q2, "a")) == NULL ){
	printf("Error opening the output file: %s\n",file_out_sint_f0_q2);
	exit(EXIT_FAILURE);
      }
      if ((fout_sint_fp_q2 = fopen( file_out_sint_fp_q2, "a")) == NULL ){
	printf("Error opening the output file: %s\n",file_out_sint_fp_q2);
	exit(EXIT_FAILURE);
      }
      
      for(int iev = 0; iev <= Nev; iev++){
	
	if( (iev == 0) || ( iev >= start_ev && iev <= end_ev ) ){
	  
	  for(int ipar = 0; ipar < npar; ipar++){
	    
	    par_plot[ipar] = parameter[ipar][iev];
	    
	  }
	  
	  fplus_sint[iq][iev] = fplus_no_hyp_num( par_plot, q2_step, MD_phys, MK_phys, Mpi_phys, Xil_phys[iev], 0.0, 1000);
	  fzero_sint[iq][iev] = fzero_no_hyp_num( par_plot, q2_step, MD_phys, MK_phys, Mpi_phys, Xil_phys[iev], 0.0, 1000);
	  
	  if( iev > 0 ){
	    
	    fprintf(fout_sint_f0_q2, "%f\n", fzero_sint[iq][iev]);
	    fprintf(fout_sint_fp_q2, "%f\n", fplus_sint[iq][iev]);
	    
	  }
	  
	} // if
      } // iev
      
      fclose(fout_sint_f0_q2);
      fclose(fout_sint_fp_q2);
      
    } // iq
    
  } // if dati_sint
  
  ////////////////// FINE CREAZIONE DATI SINTETICI

  
  
  ////////////////// OUTPUT DELLE CORREZIONI V0_hyper, Vi_hyper, f0_S_hyper

  FILE *fout_V0_hyp, *fout_Vi_hyp, *fout_f0_S_hyp;
  char file_out_V0_hyp[LEN_NAME], file_out_Vi_hyp[LEN_NAME], file_out_f0_S_hyp[LEN_NAME];

  double res, A_fplus, A_fplus_2, A_fzero, polo_fplus, polo_fzero, fse_fzero, fse_fplus;
  
  double a2, a4, V0_hyper, Vi_hyper, f0_S_hyper;
  double H1, H2, H3, H4, H5, H6, WI_violation;
  
  double temp, diff_z2, g = 0.61;
  
  double delta_fp, fzero_ph, fplus_ph, fminus_ph;
  

  if( out == 1 ){    
    
    for(int ibeta = 0; ibeta < Nbeta; ibeta++){
      for(int imusea = 0; imusea < NKL[ibeta]; imusea++){	
	for(int ith1 = 0; ith1 < Nth1; ith1++){
	  for(int ith2 = 0; ith2 < Nth2; ith2++){
	    
	    /////////// FILE MASSE
	    sprintf(open_Mpi_from_2pts, "OUTPUT_SMEAR/%s/%s/Mpi/%s/%s/Mpi.%s.m1m2_00.th_3.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str());       
	    sprintf(open_M_from_2pts_fit_1, "OUTPUT_SMEAR/%s/%s/MK/%s/%s/MK.%s.m1m2_0s.th_3.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str());
	    sprintf(open_M_from_2pts_fit_2, "OUTPUT_SMEAR/%s/%s/MD/%s/%s/MD.%s.m1m2_0c.th_3.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str());
	    
	    if ((fr_Mpi_from_2pts = fopen(open_Mpi_from_2pts, "r")) == NULL ){
	      printf("Error opening the file to read: Mpi_from_2pts\n");
	      exit(EXIT_FAILURE);
	    }
	    if ((fr_M_from_2pts_1 = fopen(open_M_from_2pts_fit_1, "r")) == NULL ){
	      printf("Error opening the file to read: Mass_from_2pts_fit_1\n");
	      exit(EXIT_FAILURE);
	    }
	    if ((fr_M_from_2pts_2 = fopen(open_M_from_2pts_fit_2, "r")) == NULL ){
	      printf("Error opening the file to read: Mass_from_2pts_fit_2\n");
	      exit(EXIT_FAILURE);
	    }
	    

	    /////////// FILE CORREZIONI IPERCUBICHE
	    sprintf(file_out_V0_hyp, "OUTPUT_SMEAR/%s/%s/hyp_corrections/V0/%s/%s/%s/V0_hyp_corr.%s.m1m2_sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
	    sprintf(file_out_Vi_hyp, "OUTPUT_SMEAR/%s/%s/hyp_corrections/Vi/%s/%s/%s/Vi_hyp_corr.%s.m1m2_sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
	    sprintf(file_out_f0_S_hyp, "OUTPUT_SMEAR/%s/%s/hyp_corrections/f0_S/%s/%s/%s/f0_S_hyp_corr.%s.m1m2_sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
	    
	    if ((fout_V0_hyp = fopen(file_out_V0_hyp, "w")) == NULL ){
	      printf("Error opening the input file: %s\n",file_out_V0_hyp);
	      exit(EXIT_FAILURE);
	    }
	    if ((fout_Vi_hyp = fopen(file_out_Vi_hyp, "w")) == NULL ){
	      printf("Error opening the input file: %s\n",file_out_Vi_hyp);
	      exit(EXIT_FAILURE);
	    }
	    if ((fout_f0_S_hyp = fopen(file_out_f0_S_hyp, "w")) == NULL ){
	      printf("Error opening the input file: %s\n",file_out_f0_S_hyp);
	      exit(EXIT_FAILURE);
	    }
	    
	    
	    for(int iev = 0; iev <= Nev; iev++){
	      
	      if( (iev == 0) || ( iev >= start_ev && iev <= end_ev ) ){
		
		/////////// CHIAMO LE MASSE
		fscanf(fr_M_from_2pts_1,"%lf\n", &MLi);
		fscanf(fr_M_from_2pts_2,"%lf\n", &MH);
		fscanf(fr_Mpi_from_2pts,"%lf\n", &Mpion);
		
		MLi = MLi/a[ibeta][iev];
		MH = MH/a[ibeta][iev];
		Mpion = Mpion/a[ibeta][iev];

		/////////// QUANTITÀ CINEMATICHE
		PLi = ((theta_value[ibeta][ith1]*PI)/L[ibeta])/a[ibeta][iev];
		PH = ((theta_value[ibeta][ith2]*PI)/L[ibeta])/a[ibeta][iev];
		
		///////////  MASSA LEGGERO
		ml = mq_l[ibeta][imusea]/(Zev[ibeta][iev]*a[ibeta][iev]);
		
		///////////  Xi_l
		Xil = (B0[iev]*ml)/pow( 4*PI*f0[iev] ,2);	      
		
		/////////// LATTICE SPACING
		a_glb = a[ibeta][iev];

		///// QUANTITÀ COSTRUITE
		a2 = pow( a_glb, 2);
		a4 = pow( a_glb, 4);

		EH = sqrt(pow( MH, 2) + 3*pow( PH, 2));
		ELi = sqrt(pow( MLi, 2) + 3*pow( PLi, 2));
		
		q0 = EH - ELi;
		qi = PH - PLi;
		q2 = pow(q0,2) - 3*pow(qi,2);
		P0 = EH + ELi;
		Pi = PH + PLi;
		q4 = pow(q0,4) + 3*pow(qi,4);
		P4 = pow(P0,4) + 3*pow(Pi,4);
		
		t0 = (MH + MLi)*pow(sqrt(MH) - sqrt(MLi), 2);
		tp = pow(MH + MLi, 2);
		z = ( sqrt(tp - q2) - sqrt(tp - t0) )/( sqrt(tp - q2) + sqrt(tp - t0) );
		z0 = ( sqrt(tp) - sqrt(tp - t0) )/( sqrt(tp) + sqrt(tp - t0) );
		z_term = (z - z0)*( 1 + (z + z0)/2. );
		
		
		res = parameter[0][iev]*( 1 + (1./2.)*Xil*log(Xil) + parameter[1][iev]*Xil + parameter[2][iev]*a2 + parameter[3][iev]*pow(Xil, 2) + parameter[42][iev]*a4*pow(LambdaQCD, 4) );
		
		polo_fzero = parameter[4][iev] + parameter[35][iev]*Xil + parameter[36][iev]*a2;
		polo_fplus = parameter[5][iev] + parameter[37][iev]*Xil + parameter[38][iev]*a2 + parameter[31][iev]*inv_MDsStar2*a4*pow(LambdaQCD, 4);
		
		fse_fplus = parameter[33][iev]*Xil*( exp(-Mpion*a_glb*L[ibeta])/pow( Mpion*a_glb*L[ibeta],  alpha_eff) );
		fse_fzero = parameter[34][iev]*Xil*( exp(-Mpion*a_glb*L[ibeta])/pow( Mpion*a_glb*L[ibeta],  alpha_eff) );
	     	
		A_fzero = parameter[6][iev] + parameter[7][iev]*Xil + parameter[8][iev]*a2 + parameter[43][iev]*a4*pow(LambdaQCD, 4);
		A_fplus = parameter[9][iev] + parameter[10][iev]*Xil + parameter[11][iev]*a2 + parameter[44][iev]*a4*pow(LambdaQCD, 4);
		//A_fplus_2 = parameter[42][iev];
		
		
		//diff_z2 = pow( z - z0 , 2);
		diff_z2 = 1;
		
		temp = Xil;
		
		H1 = parameter[12][iev]*diff_z2 + parameter[13][iev]*temp + (parameter[14][iev] + parameter[15][iev]*Xil)*(z-z0);
		
		H2 = parameter[16][iev]*diff_z2 + parameter[17][iev]*temp + (parameter[18][iev] + parameter[19][iev]*Xil)*(z-z0);
		
		H3 = parameter[20][iev]*diff_z2 + parameter[21][iev]*temp + (parameter[22][iev] + parameter[23][iev]*Xil)*(z-z0);
		
		H4 = parameter[24][iev]*diff_z2 + parameter[25][iev]*temp + (parameter[26][iev] + parameter[27][iev]*Xil)*(z-z0);
		
		H5 = parameter[28][iev]*diff_z2 + parameter[29][iev]*temp + parameter[30][iev]*(z-z0);
		
		H6 = parameter[39][iev]*diff_z2 + parameter[40][iev]*temp + parameter[41][iev]*(z-z0);
		
		
		printf("HHHH ev=%d beta=%d musea=%d th1th2=%d%d H1=%f H2=%f H3=%f H4=%f\n", iev, ibeta, imusea, ith1, ith2, H1, H2, H3, H4);
		
		printf("H1H1 ev=%d beta=%d musea=%d th1th2=%d%d H1/H3=%f H2/H4=%f\n", iev, ibeta, imusea, ith1, ith2, H1/H3, H2/H4);	     
		
		
		WI_violation = ( H5*a2*q4 + H6*a2*P4 )/( pow(MH,2) - pow(MLi,2) );
		
		fzero_ph = ( res + A_fzero*z_term )/( 1 - (1 + fse_fzero)*polo_fzero*q2  );
		fplus_ph = ( res + A_fplus*z_term )/( 1 - (1 + fse_fplus)*polo_fplus*q2  );
		fminus_ph = (fzero_ph - fplus_ph)*( (pow(MH,2) - pow(MLi,2))/q2 );
		
		V0_hyper = +a2*(pow( q0, 3)*H1 + pow( q0, 2)*P0*H2 + pow( P0, 2)*q0*H3+ pow( P0, 3)*H4);
		Vi_hyper = -a2*(pow( qi, 3)*H1 + pow( qi, 2)*Pi*H2 + pow( Pi, 2)*qi*H3+ pow( Pi, 3)*H4);
		f0_S_hyper = ( 1/(pow(MH,2) - pow(MLi,2)) )*( q0*V0_hyper - 3*qi*Vi_hyper ) + WI_violation;
		
		delta_fp = -( ( PLi - PH )*V0_hyper - ( ELi - EH )*Vi_hyper )/( 2*( EH*PLi - ELi*PH ) ); // NOTA CHE C'È IL CAMBIAMENTO DI SEGNO
	        
		
		fprintf(fout_V0_hyp, "%f\n", V0_hyper);
		fprintf(fout_Vi_hyp, "%f\n", Vi_hyper);
		fprintf(fout_f0_S_hyp, "%f\n", f0_S_hyper);
		
		
		if(iev == 0){		  
		  printf("FFFF ibeta=%d\timusea=%d\tith1=%d\tith2=%d\tdelta_fplus=%f\n", ibeta, imusea, ith1, ith2, delta_fp);
		}
		
	      }// if
	      
	      //if( (ibeta == 1) && (imusea == 2) && (ith1 == 2) && (ith2 == 5)){
	      //printf("CMCM iev=%d\tMpi=%f\tMpi=%f\tMD=%f\n", iev, Mpion*a_glb, MLi*a_glb, MH*a_glb);
	      //}
	      
	    }// iev

	    ////// CHIUSURA FILE CORREZIONI IPERCUBICHE	    
	    fclose(fout_V0_hyp);
	    fclose(fout_Vi_hyp);
	    fclose(fout_f0_S_hyp);

	    ////// CHIUSURA FILE MASSE
	    fclose(fr_Mpi_from_2pts);
	    fclose(fr_M_from_2pts_1);
	    fclose(fr_M_from_2pts_2);
	    
	  }//ith2
	}// ith1	
      }// imusea
    }// ibeta
    
  }// if( out == 1 )
  
  ////////////////// FINE OUTPUT DELLE CORREZIONI V0_hyper, Vi_hyper, f0_S_hyper
  
  ///////////////////////////// WI CHECK per l'evento 0
  
  char file_out_WIv_vs_q2[LEN_NAME], file_out_WIv_vs_q4[LEN_NAME];
  FILE *fout_WIv_vs_q2, *fout_WIv_vs_q4;
  
  double WI[Nev+1], qV, sigma_WI, temp_q2, temp_q4;
  
  for(int ibeta = 0; ibeta < Nbeta; ibeta++){
    for(int imusea = 0; imusea < NKL[ibeta]; imusea++){
      
      sprintf(file_out_WIv_vs_q2, "OUTPUT_SMEAR/%s/%s/hyp_WI/WIv_vs_q2_%s.%s.m1m2_sc.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V_2[ibeta].c_str(), mu_sea_2[imusea].c_str());
      sprintf(file_out_WIv_vs_q4, "OUTPUT_SMEAR/%s/%s/hyp_WI/WIv_vs_q4_%s.%s.m1m2_sc.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V_2[ibeta].c_str(), mu_sea_2[imusea].c_str());
      
      if ((fout_WIv_vs_q2 = fopen(file_out_WIv_vs_q2, "w")) == NULL ){
	printf("Error opening the output file: %s\n",file_out_WIv_vs_q2);
	exit(EXIT_FAILURE);
      }
      if ((fout_WIv_vs_q4 = fopen(file_out_WIv_vs_q4, "w")) == NULL ){
	printf("Error opening the output file: %s\n",file_out_WIv_vs_q4);
	exit(EXIT_FAILURE);
      }
      
      fprintf(fout_WIv_vs_q2, "@type xydy\n");
      fprintf(fout_WIv_vs_q4, "@type xydy\n");
      
      for(int ith1 = 0; ith1 < Nth1; ith1++){
	for(int ith2 = 0; ith2 < Nth2; ith2++){

	  /////////// CHIAMO LE MASSE       
	  sprintf(open_M_from_2pts_fit_1, "OUTPUT_SMEAR/%s/%s/MK/%s/%s/MK.%s.m1m2_0s.th_3.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str());
	  sprintf(open_M_from_2pts_fit_2, "OUTPUT_SMEAR/%s/%s/MD/%s/%s/MD.%s.m1m2_0c.th_3.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str());
	  
	  if ((fr_M_from_2pts_1 = fopen(open_M_from_2pts_fit_1, "r")) == NULL ){
	    printf("Error opening the file to read: %s\n",open_M_from_2pts_fit_1);
	    exit(EXIT_FAILURE);
	  }
	  if ((fr_M_from_2pts_2 = fopen(open_M_from_2pts_fit_2, "r")) == NULL ){
	    printf("Error opening the file to read: %s\n",open_M_from_2pts_fit_2);
	    exit(EXIT_FAILURE);
	  }

	  ////////// MATRIX ELEMENTS
	  sprintf(open_V0, "OUTPUT_SMEAR/%s/%s/V0/%s/%s/%s/V0.%s.m1m2_sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
	  sprintf(open_Vi, "OUTPUT_SMEAR/%s/%s/Vi/%s/%s/%s/Vi.%s.m1m2_sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
	  sprintf(open_f0_S, "OUTPUT_SMEAR/%s/%s/f0_S/%s/%s/%s/f0_S.%s.m1m2_sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
	  
	  if ((fr_V0 = fopen(open_V0, "r")) == NULL ){
	    printf("Error opening the input file: %s\n",open_V0);
	    exit(EXIT_FAILURE);
	  }
	  if ((fr_Vi = fopen(open_Vi, "r")) == NULL ){
	    printf("Error opening the input file: %s\n",open_Vi);
	    exit(EXIT_FAILURE);
	  }
	  if ((fr_f0_S = fopen(open_f0_S, "r")) == NULL ){
	    printf("Error opening the input file: %s\n",open_f0_S);
	    exit(EXIT_FAILURE);
	  }
	  	  
	  for(int iev = 0; iev <= Nev; iev++){

	    /////////// CHIAMO LE MASSE
	    fscanf(fr_M_from_2pts_1,"%lf\n", &MLi);
	    fscanf(fr_M_from_2pts_2,"%lf\n", &MH);
	    
	    MLi = MLi/a[ibeta][iev];
	    MH = MH/a[ibeta][iev];

	    /////////// CHIAMO GLI ELEMENTI DI MATRICE
	    fscanf(fr_V0,"%lf\n", &V0_glb);
	    fscanf(fr_Vi,"%lf\n", &Vi_glb);
	    fscanf(fr_f0_S,"%lf\n", &f0_S_glb);
	    
	    V0_glb = V0_glb/a[ibeta][iev];
	    Vi_glb = Vi_glb/a[ibeta][iev];

	    /////////// QUANTITÀ CINEMATICHE
	    PLi = ((theta_value[ibeta][ith1]*PI)/L[ibeta])/a[ibeta][iev];
	    PH = ((theta_value[ibeta][ith2]*PI)/L[ibeta])/a[ibeta][iev];
	    
	    EH = sqrt(pow( MH, 2) + 3*pow( PH, 2));
	    ELi = sqrt(pow( MLi, 2) + 3*pow( PLi, 2));
	    
	    q0 = EH - ELi;
	    qi = PH - PLi;
	    q2 = pow(q0,2) - 3*pow(qi,2);
	    q4 = pow(q0,4) + 3*pow(qi,4);

	    if(iev == 0){
	      temp_q2 = q2;
	      temp_q4 = q4;
	    }
	    
	    qV = q0*V0_glb - 3*qi*Vi_glb;
	    
	    WI[iev] = qV - ( pow( MH, 2) - pow( MLi, 2) )*f0_S_glb;
	    
	  }// iev

	  sigma_WI = sigma_bootstrap( WI, analysis_in, analysis_fin, Nev_an, clusterfile);

	  //////// questo vedilo alla fine
	  fprintf(fout_WIv_vs_q2, "%f %f %f\n", temp_q2, WI[0], sigma_WI);
	  fprintf(fout_WIv_vs_q4, "%f %f %f\n", temp_q4, WI[0], sigma_WI);
	  //////// questo vedilo alla fine
	  
	  ////// CHIUSURA FILE MASSE
	  fclose(fr_M_from_2pts_1);
	  fclose(fr_M_from_2pts_2);

	  ////// CHIUSURA FILE MATRIX ELEMENTS	  
	  fclose(fr_V0);
	  fclose(fr_Vi);
	  fclose(fr_f0_S);
	  
	  
	}// ith2
	fprintf(fout_WIv_vs_q2, "&\n");
	fprintf(fout_WIv_vs_q4, "&\n");	
      }// ith1
      fclose(fout_WIv_vs_q2);
      fclose(fout_WIv_vs_q4);
    }// imusea
  }// ibeta
  
  ///////////////////////////// FINE WI CHECK per l'evento 0
  
  return 0;
  
}


void chi2( int &npar, double *deriv, double &f, double *par, int iflag){
  
  int NKL[5] = { 3, 4, 4, 1, 3}, L[5] = { 32, 24, 32, 24, 48};
  double g = 0.61;
  
  double f1, f2, f3, f4;
  double V0_function, V0_mink, V0_hyper, Vi_function, Vi_mink, Vi_hyper, f0_S_function, f0_S_mink, f0_S_hyper;
  
  double fzero, fplus, fminus;
  
  double a2, a4;
  
  double res, polo_fzero, polo_fplus, A_fzero, A_fplus, A_fplus_2, fse_fplus, fse_fzero;
  
  double H1, H2, H3, H4, H5, H6, WI_violation;
  
  double temp, diff_z2;
  
  f = 0;
  
  for(int ibeta = 0; ibeta < Nbeta; ibeta++){
    for(int imusea = 0; imusea < NKL[ibeta]; imusea++){
      for(int ith1 = 0; ith1 < Nth1; ith1++){
	for(int ith2 = 0; ith2 < Nth2; ith2++){
 
#if defined(CUT) || defined(A4_CUT)
	  if( !(ibeta == 2 && imusea == 3) && !(ibeta == 1 && imusea == 2) &&  !(ibeta == 1 && imusea == 3) && !(ibeta == 3 && imusea == 0) ){  // CUT
#endif
	    
	  //if( ith2 == 3 ){  // Heavy meson at rest
	    
	    if( !(ith1 == 3 && ith2 > 3) ){

	      ///// QUANTITÀ LETTE DA FUORI
	      MLi = MLi_fit[ibeta][imusea];
	      MH = MH_fit[ibeta][imusea];
	      Mpion = Mpi_fit[ibeta][imusea];
	      
	      PH = PH_fit[ibeta][imusea][ith1][ith2];
	      PLi = PLi_fit[ibeta][imusea][ith1][ith2];

	      ml = ml_fit[ibeta][imusea];
	      Xil = Xi_l_fit[ibeta][imusea];
	      a_glb = a_fit[ibeta];

	      V0_glb = V0_fit[ibeta][imusea][ith1][ith2];
	      Vi_glb = Vi_fit[ibeta][imusea][ith1][ith2];
	      f0_S_glb = f0_S_fit[ibeta][imusea][ith1][ith2];
	      
	      sigma_V0_glb = sigma_V0_fit[ibeta][imusea][ith1][ith2];
	      sigma_Vi_glb = sigma_Vi_fit[ibeta][imusea][ith1][ith2];
	      sigma_f0_S_glb = sigma_f0_S_fit[ibeta][imusea][ith1][ith2];

	      ///// QUANTITÀ COSTRUITE
	      a2 = pow( a_glb, 2);
	      a4 = pow( a_glb, 4);

	      EH = sqrt(pow( MH, 2) + 3*pow( PH, 2));
	      ELi = sqrt(pow( MLi, 2) + 3*pow( PLi, 2));

	      q0 = EH - ELi;
	      qi = PH - PLi;
	      q2 = pow(q0,2) - 3*pow(qi,2);
	      P0 = EH + ELi;
	      Pi = PH + PLi;
	      q4 = pow(q0,4) + 3*pow(qi,4);
	      P4 = pow(P0,4) + 3*pow(Pi,4);
	      
	      t0 = (MH + MLi)*pow(sqrt(MH) - sqrt(MLi), 2);
	      tp = pow(MH + MLi, 2);
	      z = ( sqrt(tp - q2) - sqrt(tp - t0) )/( sqrt(tp - q2) + sqrt(tp - t0) );
	      z0 = ( sqrt(tp) - sqrt(tp - t0) )/( sqrt(tp) + sqrt(tp - t0) );
	      z_term = (z - z0)*( 1 + (z + z0)/2. );
	      //////////////////////////////////////

	      
	      res = par[0]*( 1 + (1./2.)*Xil*log(Xil) + par[1]*Xil + par[2]*a2 + par[3]*pow(Xil, 2) + par[42]*a4*pow(LambdaQCD, 4) );
	      
	      polo_fzero = par[4] + par[35]*Xil + par[36]*a2;
	      polo_fplus = par[5] + par[37]*Xil + par[38]*a2 + par[31]*inv_MDsStar2*a4*pow(LambdaQCD, 4);
	      
	      fse_fplus = par[33]*Xil*( exp(-Mpion*a_glb*L[ibeta])/pow( Mpion*a_glb*L[ibeta],  alpha_eff) );
	      fse_fzero = par[34]*Xil*( exp(-Mpion*a_glb*L[ibeta])/pow( Mpion*a_glb*L[ibeta],  alpha_eff) );
	      
	      A_fzero = par[6] + par[7]*Xil + par[8]*a2 + par[43]*a4*pow(LambdaQCD, 4);
	      A_fplus = par[9] + par[10]*Xil + par[11]*a2 + par[44]*a4*pow(LambdaQCD, 4);
	      //A_fplus_2 = par[42];
	      
	      fzero = ( res + A_fzero*z_term )/( 1 - (1 + fse_fzero)*polo_fzero*q2  );
	      fplus = ( res + A_fplus*z_term )/( 1 - polo_fplus*q2  );	
	      fminus = (fzero - fplus)*( (pow(MH,2) - pow(MLi,2))/q2 );
	      
	      temp = Xil;
	      
	      //diff_z2 = pow( z - z0 , 2);
	      diff_z2 = 1;
	      
	      
	      H5 = par[28]*diff_z2 + par[29]*temp + par[30]*(z-z0);
	      
	      H6 = par[39]*diff_z2 + par[40]*temp + par[41]*(z-z0);
	      
	      WI_violation = ( H5*a2*q4 )/( pow(MH,2) - pow(MLi,2) );
	      
	      
	      H1 = par[12]*diff_z2 + par[13]*temp + (par[14] + par[15]*Xil)*(z-z0) ;
	      
	      H2 = par[16]*diff_z2 + par[17]*temp + (par[18] + par[19]*Xil)*(z-z0) ;
	      
	      H3 = par[20]*diff_z2 + par[21]*temp + (par[22] + par[23]*Xil)*(z-z0) ;
	      
	      H4 = par[24]*diff_z2 + par[25]*temp + (par[26] + par[27]*Xil)*(z-z0) ;
	      
	      V0_hyper = a2*(pow( q0, 3)*H1 + pow( q0, 2)*P0*H2 + pow( P0, 2)*q0*H3+ pow( P0, 3)*H4);
	      Vi_hyper = -a2*(pow( qi, 3)*H1 + pow( qi, 2)*Pi*H2 + pow( Pi, 2)*qi*H3+ pow( Pi, 3)*H4);
	      f0_S_hyper = ( 1/(pow(MH,2) - pow(MLi,2)) )*( q0*V0_hyper -3*qi*Vi_hyper ) + WI_violation;
	      
	      V0_mink = P0*fplus + q0*fminus;
	      Vi_mink = Pi*fplus + qi*fminus;
	      f0_S_mink = fzero;
	      
	      V0_function = V0_mink + V0_hyper;
	      Vi_function = Vi_mink + Vi_hyper;
	      f0_S_function = f0_S_mink + f0_S_hyper;
	      
	      
	      f1 = pow( (V0_glb - V0_function), 2)/( pow(sigma_V0_glb, 2)) ;
	      
	      f2 = pow( (Vi_glb - Vi_function), 2)/( pow(sigma_Vi_glb, 2)) ;
	      
	      if(ith1 == 3 && ith2 == 3){
		
		f2 = 0 ;
		
	      }
	      
	      f3 = pow( f0_S_glb -  f0_S_function, 2)/( pow(sigma_f0_S_glb, 2)) ;
	      
	      
	      f = f + f1 + f2 + f3;
	      
	      
	    } // !(ith1 == 3 && ith2 >3)
	    
	    //} // if( ith2 == 3 ) Heavy meson at rest

#if defined(CUT) || defined(A4_CUT)
	    } // if( !(ibeta == 2 && imusea == 3) && !(ibeta == 1 && imusea == 2) &&  !(ibeta == 1 && imusea == 3) && !(ibeta == 3 && imusea == 0) ) CUT
#endif  
	    
	  
	} // ith2
      } // ith1
    } // imusea
  } // ibeta
  
  //f4 = pow( par[4] - 0.0390, 2)/pow( 0.000001, 2); // polo_f0 fixed 4H

#if defined(STD) || defined(CUT) || defined(NO_HYP)
  f4 = pow( par[5] - inv_MDsStar2, 2)/pow(0.00001, 2); // polo_f+ fixed
#endif
#if defined(A4) || defined(A4_CUT)
  f4 = pow( par[5] - inv_MDsStar2, 2)/pow(0.00001, 2) + pow( par[31], 2)/9 + pow( par[42], 2)/9 + pow( par[43], 2)/9 + pow( par[44], 2)/9; // polo_f+ fixed & Prior sui parametri a4
#endif
  
  f = f + f4;
  
}


double chi2_num( double *parameters){

  int NKL[5] = { 3, 4, 4, 1, 3}, L[5] = { 32, 24, 32, 24, 48}, Npts_Li_fermo = 44;
  double f = 0, temp_fens = 0, temp_f_Li_fermo = 0, f_Li_fermo_tot = 0;
  double f1, f2, f3, f4, fens[Nbeta][Nmusea], f_Li_fermo[Nbeta][Nmusea], f1_tot = 0, f2_tot = 0, f3_tot = 0;

  double g = 0.61;

  double V0_function, V0_mink, V0_hyper, Vi_function, Vi_mink, Vi_hyper, f0_S_function, f0_S_mink, f0_S_hyper;
  
  double fzero, fplus, fminus;
  
  double a2, a4;
  
  double res, polo_fzero, polo_fplus, A_fzero, A_fplus, A_fplus_2, fse_fplus, fse_fzero;
  
  double H1, H2, H3, H4, H5, H6, WI_violation;
  
  double temp, diff_z2;
  
  for(int ibeta = 0; ibeta < Nbeta; ibeta++){
    for(int imusea = 0; imusea < NKL[ibeta]; imusea++){
      for(int ith1 = 0; ith1 < Nth1; ith1++){
	for(int ith2 = 0; ith2 < Nth2; ith2++){
	  
#if defined(CUT) || defined(A4_CUT)
	  if( !(ibeta == 2 && imusea == 3) && !(ibeta == 1 && imusea == 2) &&  !(ibeta == 1 && imusea == 3) && !(ibeta == 3 && imusea == 0) ){  // CUT
#endif
	    
	  //if( ith2 == 3 ){  // Heavy meson at rest
	    
	    if( !(ith1 == 3 && ith2 > 3) ){
	      
	      ///// QUANTITÀ LETTE DA FUORI
	      MLi = MLi_fit[ibeta][imusea];
	      MH = MH_fit[ibeta][imusea];
	      Mpion = Mpi_fit[ibeta][imusea];
	      
	      PH = PH_fit[ibeta][imusea][ith1][ith2];
	      PLi = PLi_fit[ibeta][imusea][ith1][ith2];

	      ml = ml_fit[ibeta][imusea];
	      Xil = Xi_l_fit[ibeta][imusea];
	      a_glb = a_fit[ibeta];

	      V0_glb = V0_fit[ibeta][imusea][ith1][ith2];
	      Vi_glb = Vi_fit[ibeta][imusea][ith1][ith2];
	      f0_S_glb = f0_S_fit[ibeta][imusea][ith1][ith2];
	      
	      sigma_V0_glb = sigma_V0_fit[ibeta][imusea][ith1][ith2];
	      sigma_Vi_glb = sigma_Vi_fit[ibeta][imusea][ith1][ith2];
	      sigma_f0_S_glb = sigma_f0_S_fit[ibeta][imusea][ith1][ith2];

	      ///// QUANTITÀ COSTRUITE
	      a2 = pow( a_glb, 2);
	      a4 = pow( a_glb, 4);

	      EH = sqrt(pow( MH, 2) + 3*pow( PH, 2));
	      ELi = sqrt(pow( MLi, 2) + 3*pow( PLi, 2));

	      q0 = EH - ELi;
	      qi = PH - PLi;
	      q2 = pow(q0,2) - 3*pow(qi,2);
	      P0 = EH + ELi;
	      Pi = PH + PLi;
	      q4 = pow(q0,4) + 3*pow(qi,4);
	      P4 = pow(P0,4) + 3*pow(Pi,4);
	      
	      t0 = (MH + MLi)*pow(sqrt(MH) - sqrt(MLi), 2);
	      tp = pow(MH + MLi, 2);
	      z = ( sqrt(tp - q2) - sqrt(tp - t0) )/( sqrt(tp - q2) + sqrt(tp - t0) );
	      z0 = ( sqrt(tp) - sqrt(tp - t0) )/( sqrt(tp) + sqrt(tp - t0) );
	      z_term = (z - z0)*( 1 + (z + z0)/2. );
	      //////////////////////////////////////
	      
	      
	      res = parameters[0]*( 1 + (1./2.)*Xil*log(Xil) + parameters[1]*Xil  + parameters[2]*a2 + parameters[3]*pow(Xil, 2) + parameters[42]*a4*pow(LambdaQCD, 4)  );
	      
	      polo_fzero = parameters[4] + parameters[35]*Xil + parameters[36]*a2;
	      polo_fplus = parameters[5] + parameters[37]*Xil + parameters[38]*a2 + parameters[31]*inv_MDsStar2*a4*pow(LambdaQCD, 4);
	      
	      fse_fplus = parameters[33]*Xil*( exp(-Mpion*a_glb*L[ibeta])/pow( Mpion*a_glb*L[ibeta],  alpha_eff) );
	      fse_fzero = parameters[34]*Xil*( exp(-Mpion*a_glb*L[ibeta])/pow( Mpion*a_glb*L[ibeta],  alpha_eff) );
	      
	      A_fzero = parameters[6] + parameters[7]*Xil + parameters[8]*a2 + parameters[43]*a4*pow(LambdaQCD, 4);
	      A_fplus = parameters[9] + parameters[10]*Xil + parameters[11]*a2 + parameters[44]*a4*pow(LambdaQCD, 4);
	      //A_fplus_2 = parameters[42];
	      
	      fzero = ( res + A_fzero*z_term )/( 1 - (1 + fse_fzero)*polo_fzero*q2  );
	      fplus = ( res + A_fplus*z_term )/( 1 - polo_fplus*q2  );
	      fminus = (fzero - fplus)*( (pow(MH,2) - pow(MLi,2))/q2 );

	      
	      temp = Xil;
	      
	      //diff_z2 = pow( z - z0 , 2);
	      diff_z2 = 1;
	      
	      
	      H5 = parameters[28]*diff_z2 + parameters[29]*temp + parameters[30]*(z-z0);
	      
	      H6 = parameters[39]*diff_z2 + parameters[40]*temp + parameters[41]*(z-z0);
	      
	      WI_violation = ( H5*a2*q4 )/( pow(MH,2) - pow(MLi,2) );
	      
	      
	      
	      H1 = parameters[12]*diff_z2 + parameters[13]*temp + (parameters[14] + parameters[15]*Xil)*(z-z0) ;
	      
	      H2 = parameters[16]*diff_z2 + parameters[17]*temp + (parameters[18] + parameters[19]*Xil)*(z-z0) ;
	      
	      H3 = parameters[20]*diff_z2 + parameters[21]*temp + (parameters[22] + parameters[23]*Xil)*(z-z0) ;
	      
	      H4 = parameters[24]*diff_z2 + parameters[25]*temp + (parameters[26] + parameters[27]*Xil)*(z-z0) ;	
	      
	      V0_hyper = a2*(pow( q0, 3)*H1 + pow( q0, 2)*P0*H2 + pow( P0, 2)*q0*H3+ pow( P0, 3)*H4);
	      Vi_hyper = -a2*(pow( qi, 3)*H1 + pow( qi, 2)*Pi*H2 + pow( Pi, 2)*qi*H3+ pow( Pi, 3)*H4);
	      f0_S_hyper = ( 1/(pow(MH,2) - pow(MLi,2)) )*( q0*V0_hyper -3*qi*Vi_hyper ) + WI_violation;
	      
	      V0_mink = P0*fplus + q0*fminus;
	      Vi_mink = Pi*fplus + qi*fminus;
	      f0_S_mink = fzero;        
	      
	      V0_function = V0_mink + V0_hyper;
	      Vi_function = Vi_mink + Vi_hyper;
	      f0_S_function = f0_S_mink + f0_S_hyper;
	      
	      
	      f1 = pow( (V0_glb - V0_function), 2)/( pow(sigma_V0_glb, 2)) ;
	      
	      f2 = pow( (Vi_glb - Vi_function), 2)/( pow(sigma_Vi_glb, 2)) ;
	      
	      if(ith1 == 3 && ith2 == 3){
		
		f2 = 0 ;
		
		num_points = num_points - 1;
		
	      }
	      
	      f3 = pow( f0_S_glb -  f0_S_function, 2)/( pow(sigma_f0_S_glb, 2)) ;
	      
	      
	      f = f + f1 + f2 + f3;
	      
	      
	      temp_fens = temp_fens + f1 + f2 + f3;
	      
	      
	      num_points = num_points + 3; // DENTRO if
	      
	      f1_tot = f1_tot + f1;
	      f2_tot = f2_tot + f2;
	      f3_tot = f3_tot + f3;
	      
	      if( (ibeta == 2 && imusea == 0) || (ibeta == 3 && imusea == 0) || (ibeta == 4 && imusea == 1) || (ibeta == 4 && imusea == 2) ){
		if(ith1 == 3){
		  
		  temp_f_Li_fermo = temp_f_Li_fermo + f1 + f2 + f3;
		  
		  f_Li_fermo_tot = f_Li_fermo_tot + f1 + f2 + f3;
		  
		  
		} // if(ith1 == 3)
	      } // if( (ibeta == 2 && imusea == 0) || (ibeta == 3 && imusea == 0) || (ibeta == 4 && imusea == 1) )
	      
	      
	      
	    } // !(ith1 == 3 && ith2 >3)
	    
	    //} // if( ith2 == 3 ) Heavy meson at rest
	  
#if defined(CUT) || defined(A4_CUT)
	  } // if( !(ibeta == 2 && imusea == 3) && !(ibeta == 1 && imusea == 2) &&  !(ibeta == 1 && imusea == 3) && !(ibeta == 3 && imusea == 0) ) CUT
#endif  
	  
	} // ith2
      } // ith1
      
      fens[ibeta][imusea] = temp_fens;
      
      f_Li_fermo[ibeta][imusea] = temp_f_Li_fermo;
      
      temp_fens = 0;
      
      temp_f_Li_fermo = 0;
      
    } // imusea
  } // ibeta

  //f4 = pow( parameters[4] - 0.0390, 2)/pow( 0.000001, 2); // polo_f0 fixed 4H

#if defined(STD) || defined(CUT) || defined(NO_HYP)
  f4 = pow( parameters[5] - inv_MDsStar2, 2)/pow(0.00001, 2); // polo_f+ fixed
#endif
#if defined(A4) || defined(A4_CUT)
  f4 = pow( parameters[5] - inv_MDsStar2, 2)/pow(0.00001, 2) + pow( parameters[31], 2)/9 + pow( parameters[42], 2)/9 + pow( parameters[43], 2)/9 + pow( parameters[44], 2)/9; // polo_f+ fixed & Prior sui parametri a4
#endif
  
  f = f + f4;
  
  ///////////////////////////////
  
  for(int ibeta = 0; ibeta < Nbeta; ibeta++){
    for(int imusea = 0; imusea < NKL[ibeta]; imusea++){
      
      printf("KKKK ibeta=%d\timusea=%d\tchi2=%f\n", ibeta, imusea, fens[ibeta][imusea]/(num_points - (num_par - block_par)));
      
      if( (ibeta == 2 && imusea == 0) || (ibeta == 3 && imusea == 0) || (ibeta == 4 && imusea == 1) || (ibeta == 4 && imusea == 2) ){
	
	printf("RRRR ibeta=%d\timusea=%d\tchi2=%f\n", ibeta, imusea, f_Li_fermo[ibeta][imusea]/(Npts_Li_fermo - (num_par - block_par)));
	
      }
      
    }
  }
  
  printf("SSSS chi2_Li_fermo_tot=%f\n", f_Li_fermo_tot/(Npts_Li_fermo - (num_par - block_par)));
  
  printf("chi2_V0=%f\tchi2_Vi=%f\tchi2_f0_S=%f\n", f1_tot/(num_points - (num_par - block_par)), f2_tot/(num_points - (num_par - block_par)), f3_tot/(num_points - (num_par - block_par)));

  
  return f/(num_points - (num_par - block_par));
   
}


double fplus_no_hyp_num( double *parameters, double q2, double MH, double MLi, double Mpion, double Xil, double a, int L){

  double a2, a4, LambdaQCD = 0.35, g = 0.61;
  double tp, t0, z, z0, z_term;
  double res, polo_fplus, fse_fplus, A_fplus, A_fplus_2, fplus;

  a2 = pow( a, 2);
  a4 = pow( a, 4);

  t0 = (MH + MLi)*pow(sqrt(MH) - sqrt(MLi), 2);
  tp = pow(MH + MLi, 2);
  z = ( sqrt(tp - q2) - sqrt(tp - t0) )/( sqrt(tp - q2) + sqrt(tp - t0) );
  z0 = ( sqrt(tp) - sqrt(tp - t0) )/( sqrt(tp) + sqrt(tp - t0) );
  z_term = (z - z0)*( 1 + (z + z0)/2. );
  
  res = parameters[0]*( 1 + (1./2.)*Xil*log(Xil) + parameters[1]*Xil  + parameters[2]*a2 + parameters[3]*pow(Xil, 2) + parameters[42]*a4*pow(LambdaQCD, 4)  );
  
  polo_fplus = parameters[5] + parameters[37]*Xil + parameters[38]*a2 + parameters[31]*inv_MDsStar2*a4*pow(LambdaQCD, 4);
  
  fse_fplus = parameters[33]*Xil*( exp(-Mpion*L)/pow( Mpion*L,  alpha_eff) ); // GLI PASSO SOLO "L" E NON "a*L" PERCHÈ NON È DEFINITO AL CONTINUO
  
  A_fplus = parameters[9] + parameters[10]*Xil + parameters[11]*a2 + parameters[44]*a4*pow(LambdaQCD, 4);
  
  //A_fplus_2 = parameters[42];
  

  fplus = ( res + A_fplus*z_term )/( 1 - (1 + fse_fplus)*polo_fplus*q2  );
  
  return fplus;
  
}


double fzero_no_hyp_num( double *parameters, double q2, double MH, double MLi, double Mpion, double Xil, double a, int L){
  
  double a2, a4, LambdaQCD = 0.35, g = 0.61;
  double tp, t0, z, z0, z_term;
  double res, polo_fzero, fse_fzero, A_fzero, fzero;

  a2 = pow( a, 2);
  a4 = pow( a, 4);

  t0 = (MH + MLi)*pow(sqrt(MH) - sqrt(MLi), 2);
  tp = pow(MH + MLi, 2);
  z = ( sqrt(tp - q2) - sqrt(tp - t0) )/( sqrt(tp - q2) + sqrt(tp - t0) );
  z0 = ( sqrt(tp) - sqrt(tp - t0) )/( sqrt(tp) + sqrt(tp - t0) );
  z_term = (z - z0)*( 1 + (z + z0)/2. );
  
  res = parameters[0]*( 1 + (1./2.)*Xil*log(Xil) + parameters[1]*Xil  + parameters[2]*a2 + parameters[3]*pow(Xil, 2) + parameters[42]*a4*pow(LambdaQCD, 4)  );
  
  polo_fzero = parameters[4] + parameters[35]*Xil + parameters[36]*a2;
  
  fse_fzero = parameters[34]*Xil*( exp(-Mpion*L)/pow( Mpion*L,  alpha_eff) ); // GLI PASSO SOLO "L" E NON "a*L" PERCHÈ NON È DEFINITO AL CONTINUO

  A_fzero = parameters[6] + parameters[7]*Xil + parameters[8]*a2 + parameters[43]*a4*pow(LambdaQCD, 4);

  
  fzero = ( res + A_fzero*z_term )/( 1 - (1 + fse_fzero)*polo_fzero*q2  );
  
  return fzero;
  
}

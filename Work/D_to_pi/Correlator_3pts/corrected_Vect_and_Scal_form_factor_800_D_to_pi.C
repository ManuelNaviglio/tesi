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
const int npar = 2;

//////////////////////////////

int f0_with_S = 1;
int f0_without_S = 0;

//////////////////////////////

//////////////////////////////
   
int energia_da_sinhDR = 0;
int energia_da_stdDR = 1;
int energia_dal_fit = 0;

//////////////////////////////


string beta_V[5] =  {"1.90/32", "1.90/24", "1.95/32", "1.95/24", "2.10/48"};
string mu_sea_2[4] = {"musea_0", "musea_1", "musea_2", "musea_3"};
string m[15] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "A", "B", "C", "D", "E"};
string th_[7] = {"0", "1", "2", "3", "4", "5", "6"};

int ibeta, imusea, ith1, ith2, no_correction;

/*
double theta_value[5][7]={{-1.733, -0.933, -0.400, 0.0, 0.400, 0.933, 1.733}, {-1.300, -0.700, -0.300, 0.0, 0.300, 0.700, 1.300}, 
                          {-1.588, -0.854, -0.366, 0.0, 0.366, 0.854, 1.588}, {-1.191, -0.641, -0.275, 0.0, 0.275, 0.641, 1.191}, 
                          {-1.832, -0.986, -0.424, 0.0, 0.424, 0.986, 1.832}};
*/

double q2_glb, E_H_glb, E_l_glb, pl_glb, pH_glb, M_l_glb, M_H_glb;
double V0_glb, Vi_glb, f0_S_glb;
double sigma_V0_glb, sigma_Vi_glb, sigma_f0_S_glb;

void chi2( int &, double *, double &, double *, int);
double chi2_num( double *);

int main(){
  
  ////////  LETTURA Ultimate Input
  
  double mlight[Nev+1], mstrange[Nev+1], mcharm[Nev+1], a[Nbeta][Nev+1], ainv[Nbeta][Nev+1], r0[Nev+1], Zev[Nbeta][Nev+1], ZTev[Nbeta][Nev+1], f0[Nev+1], B0[Nev+1], fkfpi[Nev+1];
  int  iboot[Nbeta][Nmusea][Nev+1];
  
  lettura_ultimate_input( mlight, mstrange, mcharm, a, ainv, r0, Zev, iboot, f0, B0, fkfpi, ZTev);
  
  ////////  FINE LETTURA Ultimate Input
  
  
  FILE *fin;
  
  if ((fin = fopen("Input_corr_3pts/file_input_corr_3pts.out", "r")) == NULL ){
    printf("Error opening the input file!!\n");
    exit(EXIT_FAILURE);
  }
  fscanf(fin, "%d %d %d %d %d", &ibeta, &imusea, &ith1, &ith2, &no_correction); 
  fclose(fin);
  
  //////////////////////////////
  
  string strategy[4] = { "light_to_Heavy", "Heavy_to_light", "Hl_sym", "Hl_sym_av"};
  int iE, istrategy = 3; // questo fissa la cartella Hl_sym_av
  
  string dir_E[3] = {"E_from_sinhDR", "E_from_stdDR", "E_from_fit"};
  
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
  
  ////////////// LEGGO GLI EL. DI MATRICE, f0_S E LE LORO CORREZIONI
  
  FILE *fr_V0, *fr_Vi, *fr_f0_S, *fr_V0_correction, *fr_Vi_correction, *fr_f0_S_correction;
  char open_V0[LEN_NAME], open_Vi[LEN_NAME], open_f0_S[LEN_NAME];
  char open_V0_correction[LEN_NAME], open_Vi_correction[LEN_NAME], open_f0_S_correction[LEN_NAME];

  double V0[Nev+1], V0_correction[Nev+1], V0_corrected[Nev+1], sigma_V0_corrected[Nanalysis];
  double Vi[Nev+1], Vi_correction[Nev+1], Vi_corrected[Nev+1], sigma_Vi_corrected[Nanalysis];
  double f0_S[Nev+1], f0_S_correction[Nev+1], f0_S_corrected[Nev+1], sigma_f0_S_corrected[Nanalysis];
  
  sprintf(open_V0, "OUTPUT_SMEAR/%s/%s/V0/%s/%s/%s/V0.%s.m1m2_0c.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  sprintf(open_Vi, "OUTPUT_SMEAR/%s/%s/Vi/%s/%s/%s/Vi.%s.m1m2_0c.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  sprintf(open_f0_S, "OUTPUT_SMEAR/%s/%s/f0_S/%s/%s/%s/f0_S.%s.m1m2_0c.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  
  sprintf(open_Vi_correction, "OUTPUT_SMEAR/%s/%s/hyp_corrections/Vi/%s/%s/%s/Vi_hyp_corr.%s.m1m2_0c.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  sprintf(open_V0_correction, "OUTPUT_SMEAR/%s/%s/hyp_corrections/V0/%s/%s/%s/V0_hyp_corr.%s.m1m2_0c.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  sprintf(open_f0_S_correction, "OUTPUT_SMEAR/%s/%s/hyp_corrections/f0_S/%s/%s/%s/f0_S_hyp_corr.%s.m1m2_0c.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), th_[ith1].c_str(), th_[ith2].c_str());

  if ((fr_Vi = fopen(open_Vi, "r")) == NULL ){
    printf("Error opening the input file: open_Vi\n");
    exit(EXIT_FAILURE);
  }
  if ((fr_V0 = fopen(open_V0, "r")) == NULL ){
    printf("Error opening the input file: open_V0\n");
    exit(EXIT_FAILURE);
  }
  if ((fr_f0_S = fopen(open_f0_S, "r")) == NULL ){
    printf("Error opening the input file: open_f0_S\n");
    exit(EXIT_FAILURE);
  }
  if ((fr_Vi_correction = fopen(open_Vi_correction, "r")) == NULL ){
    printf("Error opening the input file: open_Vi\n");
    exit(EXIT_FAILURE);
  }
  if ((fr_V0_correction = fopen(open_V0_correction, "r")) == NULL ){
    printf("Error opening the input file: open_V0\n");
    exit(EXIT_FAILURE);
  }
  if ((fr_f0_S_correction = fopen(open_f0_S_correction, "r")) == NULL ){
    printf("Error opening the input file: open_f0_S\n");
    exit(EXIT_FAILURE);
  }

  for(int iev = 0; iev <= Nev; iev++){

    fscanf(fr_V0,"%lf\n", &V0[iev]);
    fscanf(fr_Vi,"%lf\n", &Vi[iev]);
    fscanf(fr_f0_S,"%lf\n", &f0_S[iev]);
    
    V0[iev] = V0[iev]/a[ibeta][iev];  // TRASFORMATO IN GeV
    Vi[iev] = Vi[iev]/a[ibeta][iev];  // TRASFORMATO IN GeV
    
    fscanf(fr_V0_correction,"%lf\n", &V0_correction[iev]);      // QUESTE SONO GIÀ IN GeV
    fscanf(fr_Vi_correction,"%lf\n", &Vi_correction[iev]);      // QUESTE SONO GIÀ IN GeV
    fscanf(fr_f0_S_correction,"%lf\n", &f0_S_correction[iev]);
    
    if(no_correction == 1){
      
      V0_correction[iev] = 0;
      Vi_correction[iev] = 0;
      f0_S_correction[iev] = 0;
    }
    
    V0_corrected[iev] = V0[iev] - V0_correction[iev];
    Vi_corrected[iev] = Vi[iev] - Vi_correction[iev];
    f0_S_corrected[iev] = f0_S[iev] - f0_S_correction[iev];
    
  }// iev

  fclose(fr_V0);
  fclose(fr_Vi);
  fclose(fr_f0_S);
  fclose(fr_V0_correction);
  fclose(fr_Vi_correction);
  fclose(fr_f0_S_correction);

  for(int ianalysis = 0; ianalysis < Nanalysis; ianalysis++){
    
    sigma_V0_corrected[ianalysis] = sigma_JK_modified_2(V0_corrected, ianalysis*Nev_an + 1, ianalysis*Nev_an + Nev_an, clusterfile);
    sigma_Vi_corrected[ianalysis] = sigma_JK_modified_2(Vi_corrected, ianalysis*Nev_an + 1, ianalysis*Nev_an + Nev_an, clusterfile);
    sigma_f0_S_corrected[ianalysis] = sigma_JK_modified_2(f0_S_corrected, ianalysis*Nev_an + 1, ianalysis*Nev_an + Nev_an, clusterfile);
    
  }// ianalysis

  ////////////// FINE LEGGO GLI EL. DI MATRICE E LE LORO CORREZIONI



  /////////// CHIAMO LE MASSE

  FILE *fr_M_from_2pts_1, *fr_M_from_2pts_2;
  char open_M_from_2pts_fit_1[LEN_NAME], open_M_from_2pts_fit_2[LEN_NAME];

  double mass_2pts_m1[Nev+1], mass_2pts_m2[Nev+1];

  sprintf(open_M_from_2pts_fit_1, "OUTPUT_SMEAR/%s/%s/Mpi/%s/%s/Mpi.%s.m1m2_00.th_3.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str());
  sprintf(open_M_from_2pts_fit_2, "OUTPUT_SMEAR/%s/%s/MD/%s/%s/MD.%s.m1m2_%sc.th_3.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m[0].c_str());

  if ((fr_M_from_2pts_1 = fopen(open_M_from_2pts_fit_1, "r")) == NULL ){
    printf("Error opening the file to read: Mass_from_2pts_fit_1\n");
    exit(EXIT_FAILURE);
  }
  if ((fr_M_from_2pts_2 = fopen(open_M_from_2pts_fit_2, "r")) == NULL ){
    printf("Error opening the file to read: Mass_from_2pts_fit_2\n");
    exit(EXIT_FAILURE);
  }

  for(int iev = 0; iev <= Nev; iev++){

    fscanf(fr_M_from_2pts_1,"%lf\n", &mass_2pts_m1[iev]);  // IN lattice unit
    fscanf(fr_M_from_2pts_2,"%lf\n", &mass_2pts_m2[iev]);  // IN lattice unit
    
    mass_2pts_m1[iev] = mass_2pts_m1[iev]/a[ibeta][iev];  // IN GeV
    mass_2pts_m2[iev] = mass_2pts_m2[iev]/a[ibeta][iev];  // IN GeV
    
  }// iev
  
  fclose(fr_M_from_2pts_1);
  fclose(fr_M_from_2pts_2);


  /////////// FINE CHIAMO LE MASSE

  /////////// COSTRUIAMO LE QUANTITÀ CINEMATICHE

  double th1_val, th2_val, momentum_1[Nev+1], momentum_2[Nev+1];

  double energy_2pts_m1[Nev+1], energy_2pts_m2[Nev+1];

  double q2[Nev+1], q0[Nev+1], qi[Nev+1];

  th1_val = theta_value[ibeta][ith1];
  th2_val = theta_value[ibeta][ith2];

  for(int iev = 0; iev <= Nev; iev++){

    momentum_1[iev] = (th1_val*PI)/(L[ibeta]*a[ibeta][iev]);                                      // In GeV
    momentum_2[iev] = (th2_val*PI)/(L[ibeta]*a[ibeta][iev]);                                      // In GeV 

    energy_2pts_m1[iev] = sqrt( pow(mass_2pts_m1[iev] ,2) + 3*pow(momentum_1[iev] ,2) );          // In GeV
    energy_2pts_m2[iev] = sqrt( pow(mass_2pts_m2[iev] ,2) + 3*pow(momentum_2[iev] ,2) );          // In GeV
     
    q0[iev] = (energy_2pts_m2[iev] - energy_2pts_m1[iev]);                                        // In GeV
    qi[iev] = (momentum_2[iev] - momentum_1[iev]);                                                // In GeV

    q2[iev] = pow(q0[iev] ,2) - 3*pow(qi[iev] ,2);                                                // In GeV

  }// iev

  /////////// FINE QUANTITÀ CINEMATICHE

  
  /////////// INIZIO FIT DI f+ ED F0

  double fzero[Nev+1], fplus[Nev+1], fminus[Nev+1];

  double chi2_f;

  double outpar[npar], err[npar];
  double par[npar] = {0.00, 0.00};
  double step[npar] = {0.01, 0.01};
  double min[npar] = {0.00, 0.00};
  double max[npar] = {0.00, 0.00};
  string cpar[npar] = {"f_zero - f_plus", "f_plus"};  // par[0] è dato come differenza per motivi di convergenza del fit

  int ian;
  
  if(ith1 != 3 || ith2 != 3){
    
    for(int iev = 0; iev <= Nev; iev++){
      
      if(iev == 0){
	ian = 0;
      } else if(iev != 0){
	ian = ((iev-1)/Nev_an);
      }
       
      q2_glb = q2[iev];
      E_H_glb = energy_2pts_m2[iev];
      E_l_glb = energy_2pts_m1[iev];
      pH_glb = momentum_2[iev];
      pl_glb = momentum_1[iev];
      M_H_glb = mass_2pts_m2[iev];
      M_l_glb = mass_2pts_m1[iev];
      
      V0_glb = V0_corrected[iev];
      Vi_glb = Vi_corrected[iev];
      f0_S_glb = f0_S_corrected[iev];
      
      sigma_V0_glb = sigma_V0_corrected[ian];
      sigma_Vi_glb = sigma_Vi_corrected[ian];
      sigma_f0_S_glb = sigma_f0_S_corrected[ian];
      
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
      
      fzero[iev] = outpar[0]+outpar[1];
      fplus[iev] = outpar[1];
      fminus[iev] = (fzero[iev]-fplus[iev])*( ( pow(M_H_glb ,2)-pow(M_l_glb ,2) )/q2_glb  );
      
      chi2_f = chi2_num(outpar);
      
      printf("CCCC chi2_ev[%d] = %f\n", iev, chi2_f);
      
    } // iev
    
  } // if(ith1 != 3 || ith2 != 3){
  
  /////////// FINE FIT DI f+ ED F0

  /////////// NEL CASO ith1 == 3 && ith2 == 3
  
  double fzero_ferma[Nev+1], sigma_fzero_ferma[Nanalysis];
  double num = 0, den = 0;
  
  if(ith1 == 3 && ith2 == 3){
    
    for(int iev = 0; iev <= Nev; iev++){
      
      fzero_ferma[iev] = V0_corrected[iev]/(mass_2pts_m2[iev] + mass_2pts_m1[iev]);
      fplus[iev] = 0;
      fminus[iev] = 0;
    }// iev
    
    for(int ianalysis = 0; ianalysis < Nanalysis; ianalysis++){
      sigma_fzero_ferma[ianalysis] = sigma_JK_modified_2(fzero_ferma, ianalysis*Nev_an + 1, ianalysis*Nev_an + Nev_an, clusterfile);
    }
    
    for(int iev = 0; iev <= Nev; iev++){
      
      if(iev == 0){
	ian = 0;
      } else if(iev != 0){
	ian = ((iev-1)/Nev_an);
      }
      
      num = ( fzero_ferma[iev]*(1/pow(sigma_fzero_ferma[ian] ,2)) ) + ( f0_S_corrected[iev]*(1/pow(sigma_f0_S_corrected[ian] ,2)) );
      den = ( 1/pow(sigma_fzero_ferma[ian] ,2) ) + ( 1/pow(sigma_f0_S_corrected[ian] ,2) );
      
      fzero[iev] = num/den;
      
    }// iev
    
  } // if(ith1 == 3 && ith2 == 3)
  
  /////////// FINE CASO ith1 == 3 && ith2 == 3
    
  /////////// INIZIO OUTPUT
  
  char file_out_fzero[LEN_NAME], file_out_fplus[LEN_NAME], file_out_fminus[LEN_NAME];
  FILE *fout_fzero, *fout_fplus, *fout_fminus;
  
  if(no_correction == 0){
    
    sprintf(file_out_fzero, "OUTPUT_SMEAR/%s/%s/fzero_corrected/%s/%s/%s/fzero_corrected.%s.m1m2_%sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m[0].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
    sprintf(file_out_fplus, "OUTPUT_SMEAR/%s/%s/fplus_corrected/%s/%s/%s/fplus_corrected.%s.m1m2_%sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m[0].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
    sprintf(file_out_fminus, "OUTPUT_SMEAR/%s/%s/fminus_corrected/%s/%s/%s/fminus_corrected.%s.m1m2_%sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m[0].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
    
  }else if(no_correction == 1){
    
    sprintf(file_out_fzero, "OUTPUT_SMEAR/%s/%s/fzero_corrected/%s/%s/%s/fzero_no_corrected.%s.m1m2_%sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m[0].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
    sprintf(file_out_fplus, "OUTPUT_SMEAR/%s/%s/fplus_corrected/%s/%s/%s/fplus_no_corrected.%s.m1m2_%sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m[0].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
    sprintf(file_out_fminus, "OUTPUT_SMEAR/%s/%s/fminus_corrected/%s/%s/%s/fminus_no_corrected.%s.m1m2_%sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m[0].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
    
  }
  
  if ((fout_fzero = fopen(file_out_fzero, "w")) == NULL ){
    printf("Error opening the output file: file_out_fzero\n");
    exit(EXIT_FAILURE);
  }
  for(int iev = 0; iev <= Nev; iev++){
    fprintf(fout_fzero,"%f\n", fzero[iev]);
  }
  fclose(fout_fzero);
  
  
  if(ith1 !=3 || ith2 !=3 ){
    
    if ((fout_fplus = fopen(file_out_fplus, "w")) == NULL ){
      printf("Error opening the output file: file_out_fplus\n");
      exit(EXIT_FAILURE);
    }
    if ((fout_fminus = fopen(file_out_fminus, "w")) == NULL ){
      printf("Error opening the output file: file_out_fminus\n");
      exit(EXIT_FAILURE);
    }
       
    for(int iev = 0; iev <= Nev; iev++){
      fprintf(fout_fplus,"%f\n", fplus[iev]);
      fprintf(fout_fminus,"%f\n", fminus[iev]);
    }
    
    fclose(fout_fplus);
    fclose(fout_fminus);       
    
  } // if(ith1 !=3 || ith2 !=3) 
  
  /////////// FINE OUTPUT

  return 0;
  
}

void chi2( int &npar, double *deriv, double &f, double *par, int iflag){
  
  double f1, f2, f3;
  double V0_function, Vi_function;

  f = 0;


  V0_function = par[1]*(E_H_glb + E_l_glb) + (par[0])*( ( pow(M_H_glb ,2)-pow(M_l_glb ,2) )/q2_glb  )*(E_H_glb - E_l_glb);

  Vi_function = par[1]*(pH_glb + pl_glb) + (par[0])*( ( pow(M_H_glb ,2)-pow(M_l_glb ,2) )/q2_glb  )*(pH_glb - pl_glb);


  f1 = pow( (V0_glb - V0_function), 2)/( pow(sigma_V0_glb, 2)) ;

  f2 = pow( (Vi_glb - Vi_function), 2)/( pow(sigma_Vi_glb, 2)) ;

  f3 = pow( f0_S_glb -  (par[0]+par[1]), 2)/( pow(sigma_f0_S_glb, 2)) ;

  f = f1 + f2 + f3;
  
}


double chi2_num( double *parameters){
  
  double f1, f2, f3;
  double V0_function, Vi_function;

  double f;


  V0_function = parameters[1]*(E_H_glb + E_l_glb) + (parameters[0])*( ( pow(M_H_glb ,2)-pow(M_l_glb ,2) )/q2_glb  )*(E_H_glb - E_l_glb);

  Vi_function = parameters[1]*(pH_glb + pl_glb) + (parameters[0])*( ( pow(M_H_glb ,2)-pow(M_l_glb ,2) )/q2_glb  )*(pH_glb - pl_glb);


  f1 = pow( (V0_glb - V0_function), 2)/( pow(sigma_V0_glb, 2)) ;

  f2 = pow( (Vi_glb - Vi_function), 2)/( pow(sigma_Vi_glb, 2)) ;

  f3 = pow( f0_S_glb -  (parameters[0]+parameters[1]), 2)/( pow(sigma_f0_S_glb, 2)) ;

  f = f1 + f2 + f3;

  return f;
  
}

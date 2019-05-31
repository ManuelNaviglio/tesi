#include <TF1.h>
#include <fstream>
#include "lettura_ult_inp.h"
#include "stat_analysis_func.h"
#include "spline_func.h"
#include "definitions.h"

using namespace std;

#define LEN_NAME 1024
#define PI 3.141592653589793

int i_mc0 = 4;
int i_ms0 = 1;

double phys_mc;
double phys_ms;

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

  ////////  LETTURA DEL FILE DI INPUT
  int ibeta, imusea, ith1, ith2, im;
  FILE *f_in;

  if ((f_in = fopen("Input_corr_3pts/file_input_corr_3pts.out", "r")) == NULL ){
    printf("Error opening the input file  !!\n");
    exit(EXIT_FAILURE);
  }
  fscanf(f_in, "%d %d %d %d %d", &ibeta, &imusea, &im, &ith1, &ith2);
  fclose(f_in);
  ////////  FINE LETTURA DEL FILE DI INPUT

  /////////////// COSTRUISCO L'ARRAY mq
  
  double mq[Nbeta][Nmasses] = {{0.0030, 0.01800, 0.02200, 0.02600, 0.21256, 0.25000, 0.29404, 0.34583, 0.40675, 0.47840, 0.56267, 0.66178, 0.77836, 0.91546, 1.07672},
			       {0.0040, 0.01800, 0.02200, 0.02600, 0.21256, 0.25000, 0.29404, 0.34583, 0.40675, 0.47840, 0.56267, 0.66178, 0.77836, 0.91546, 1.07672},
			       {0.0025, 0.01550, 0.01900, 0.02250, 0.18705, 0.22000, 0.25875, 0.30433, 0.35794, 0.42099, 0.49515, 0.58237, 0.68495, 0.80561, 0.94752},
			       {0.0085, 0.01550, 0.01900, 0.02250, 0.18705, 0.22000, 0.25875, 0.30433, 0.35794, 0.42099, 0.49515, 0.58237, 0.68495, 0.80561, 0.94752},
			       {0.0015, 0.01230, 0.01500, 0.01770, 0.14454, 0.17000, 0.19995, 0.23517, 0.27659, 0.32531, 0.38262, 0.45001, 0.52928, 0.62252, 0.73217}};
    
  mq[ibeta][0] = mq_l[ibeta][imusea]; // questo è il motivo per cui non posso inserirlo in "definitions.h"
  
  /////////////// COSTRUISCO L'ARRAY mq
  
  ////////  LETTURA Ultimate Input
    
  double mlight[Nev+1], mstrange[Nev+1], mcharm[Nev+1], a[Nbeta][Nev+1], ainv[Nbeta][Nev+1], r0[Nev+1], Zev[Nbeta][Nev+1], ZTev[Nbeta][Nev+1], f0[Nev+1], B0[Nev+1], fkfpi[Nev+1];
  int  iboot[Nbeta][Nmusea][Nev+1];

  lettura_ultimate_input( mlight, mstrange, mcharm, a, ainv, r0, Zev, iboot, f0, B0, fkfpi, ZTev);
  
  ////////  FINE LETTURA Ultimate Input

  
  ///////// BOOTSTRAP DELLA MASSA Mpi

#ifdef SPLINE_MPI
  
  double Mpi[Nev+1];
  char file_open_Mpi[LEN_NAME], file_out_Mpi[LEN_NAME];
  FILE *fout_Mpi;
  
  sprintf(file_open_Mpi, "../Correlator_2pts/OUTPUT_SMEAR/M/fit_data/%s/Fit_%s/fit_2pts_M.%s.m1m2_%s%s.th_3.sme_30_30.r1r2_Av_00_11.P5P5.out", beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[0].c_str(), m_[0].c_str());

  sprintf(file_out_Mpi, "OUTPUT_SMEAR/%s/%s/Mpi/%s/%s/Mpi.%s.m1m2_00.th_3.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str());
    
  bootstrap_sampling( Mpi, file_open_Mpi, ibeta, imusea, clusterfile);

  if ((fout_Mpi = fopen(file_out_Mpi, "w")) == NULL ){
    printf("Error opening the output file Mpi!!\n");
    exit(EXIT_FAILURE);
  }
  
  for(int iev = 0; iev <= Nev; iev++ ){
    fprintf(fout_Mpi, "%+f\n", Mpi[iev]);
  }// iev
  
  fclose(fout_Mpi);

#endif
  
  ///////// FINE BOOTSTRAP DELLA MASSA Mpi





  double *temp_x;

  ///////// SPLINE MD

#ifdef SPLINE_MD
  
  ///////// CREO I BOOTSTRAP
  
  double MD[Nmass_charm][Nev+1];
  char file_open_MD[LEN_NAME];
  
  for(int imc = i_mc0; imc < i_mc0 + Nmass_charm; imc++){
    
    sprintf(file_open_MD, "../Correlator_2pts/OUTPUT_SMEAR/M/fit_data/%s/Fit_%s/fit_2pts_M.%s.m1m2_%s%s.th_3.sme_30_30.r1r2_Av_00_11.P5P5.out", beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[0].c_str(), m_[imc].c_str());
    
    bootstrap_sampling( MD[imc-i_mc0], file_open_MD, ibeta, imusea, clusterfile);
    
  }// imc
  
  ///////// FINE CREO I BOOTSTRAP
  
  ///////// FACCIO LA SPLINE

  temp_x = (double*)malloc(sizeof(double)*(Nmass_charm));
  double *temp_y_MD = (double*)malloc(sizeof(double)*(Nmass_charm));
  double *B_MD = (double*)malloc(sizeof(double)*(Nmass_charm)), *C_MD = (double*)malloc(sizeof(double)*(Nmass_charm));
  double *yspline_MD = (double*)malloc(sizeof(double)*(Nev+1));
  
  for(int iev = 0; iev <= Nev; iev++ ){    
    for(int imc = i_mc0; imc < i_mc0 + Nmass_charm; imc++){
      
      temp_x[imc-i_mc0] = (mq[ibeta][imc]*ainv[ibeta][0])/Zev[ibeta][0];  // prendo le masse in GeV
      //temp_x[imc-i_mc0] = mq[ibeta][imc];  // prendo le masse in unità di reticolo
      
      temp_y_MD[imc-i_mc0] = MD[imc-i_mc0][iev];
      
    }// imc
    
    phys_mc = mcharm[iev];  // prendo le masse in GeV
    //phys_mc = (mcharm[iev]/ainv[ibeta][iev])*Zev[ibeta][iev];  // prendo le masse in unità di reticolo
    
    spline2(Nmass_charm, temp_x, temp_y_MD, B_MD, C_MD);
    yspline_MD[iev] = yspline2(Nmass_charm, temp_x, temp_y_MD, B_MD, C_MD, phys_mc);
    
  }// iev
  
  free(temp_x);
  free(temp_y_MD);
  free(B_MD);
  free(C_MD);
  
  ///////// FINE FACCIO LA SPLINE
  
  /////////// CHECK GRAFICO

  //////// check per l'evento iev_check_MD = 2
  int iev_check_MD = 2;
  
  FILE *fout_graph_check_MD;
  char file_out_graph_check_MD[LEN_NAME];

  sprintf(file_out_graph_check_MD, "OUTPUT_SMEAR/%s/%s/check_spline/MD/MD_spline_check.%s.%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V_2[ibeta].c_str(), mu_sea_2[imusea].c_str());
  
  if ((fout_graph_check_MD = fopen(file_out_graph_check_MD, "w")) == NULL ){
    printf("Error opening the output file file_out_graph_check_MD\n");
    exit(EXIT_FAILURE);
  }
  
  fprintf(fout_graph_check_MD, "@type xydy\n");
  
  for(int imc = i_mc0; imc < i_mc0 + Nmass_charm; imc++){
    
    //  prendendo le masse in GeV
    fprintf(fout_graph_check_MD, "%f %f %f\n", (mq[ibeta][imc]*ainv[ibeta][iev_check_MD])/Zev[ibeta][iev_check_MD], MD[imc-i_mc0][iev_check_MD], sigma_JK_modified_2(MD[imc-i_mc0], 1, 100, clusterfile));
    
    //  prendendo le masse in unità di reticolo
    //fprintf(fout_graph_check_MD, "%f %f %f\n", mq[ibeta][imc], MD[imc-i_mc0][iev_check_MD], sigma_JK_modified_2(MD[imc-i_mc0], 1, 100, clusterfile));
    
  }
  
  fprintf(fout_graph_check_MD, "&\n");
  
  //  prendendo le masse in GeV
  fprintf(fout_graph_check_MD, "%f %f %f\n", mcharm[iev_check_MD], yspline_MD[iev_check_MD], sigma_JK_modified_2(yspline_MD, 1, 100, clusterfile));
  
  //  prendendo le masse in unità di reticolo
  //fprintf(fout_graph_check_MD, "%f %f %f\n", (mcharm[iev_check_MD]/ainv[ibeta][iev_check_MD])*Zev[ibeta][iev_check_MD], yspline_MD[iev_check_MD], sigma_JK_modified_2(yspline_MD, 1, 100, clusterfile));
  
  fclose(fout_graph_check_MD);
  
  /////////// FINE CHECK GRAFICO

  ///////// OUTPUT

  char file_out_MD[LEN_NAME];
  FILE *fout_MD;
  
  sprintf(file_out_MD, "OUTPUT_SMEAR/%s/%s/MD/%s/%s/MD.%s.m1m2_%sc.th_3.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[0].c_str());

  if ((fout_MD = fopen(file_out_MD, "w")) == NULL ){
    printf("Error opening the output file MD!!\n");
    exit(EXIT_FAILURE);
  }
  
  for(int iev = 0; iev <= Nev; iev++ ){
    fprintf(fout_MD, "%+f\n", yspline_MD[iev]);
  }// iev
  
  fclose(fout_MD);

  free(yspline_MD);
  
  ///////// FINE OUTPUT

#endif
  
  ///////// FINE SPLINE MD
  



  ///////// SPLINE MK

#ifdef SPLINE_MK
  
  ///////// CREO I BOOTSTRAP
  
  double MK[Nmass_strange][Nev+1];
  char file_open_MK[LEN_NAME];
  
  for(int ims = i_ms0; ims < i_ms0 + Nmass_strange; ims++ ){
    
    sprintf(file_open_MK, "../Correlator_2pts/OUTPUT_SMEAR/M/fit_data/%s/Fit_%s/fit_2pts_M.%s.m1m2_%s%s.th_3.sme_30_30.r1r2_Av_00_11.P5P5.out", beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[0].c_str(), m_[ims].c_str());
    
    bootstrap_sampling( MK[ims-i_ms0], file_open_MK, ibeta, imusea, clusterfile);
    
  }// ims
  
  ///////// FINE CREO I BOOTSTRAP

  ///////// FACCIO LA SPLINE
  
  temp_x = (double*)malloc(sizeof(double)*(Nmass_strange));
  double *temp_y_MK = (double*)malloc(sizeof(double)*(Nmass_strange));
  double *B_MK = (double*)malloc(sizeof(double)*(Nmass_strange)), *C_MK = (double*)malloc(sizeof(double)*(Nmass_strange));
  double *yspline_MK = (double*)malloc(sizeof(double)*(Nev+1));
  
  for(int iev = 0; iev <= Nev; iev++ ){    
    for(int ims = i_ms0; ims < i_ms0 + Nmass_strange; ims++){
      
      temp_x[ims-i_ms0] = (mq[ibeta][ims]*ainv[ibeta][0])/Zev[ibeta][0];  // prendo le masse in GeV
      //temp_x[ims-i_ms0] = mq[ibeta][ims];  // prendo le masse in unità di reticolo
      
      temp_y_MK[ims-i_ms0] = MK[ims-i_ms0][iev];
      
    }// ims
    
    phys_ms = mstrange[iev];  // prendo le masse in GeV
    //phys_ms = (mstrange[iev]/ainv[ibeta][iev])*Zev[ibeta][iev];  // prendo le masse in unità di reticolo
    
    spline2(Nmass_strange, temp_x, temp_y_MK, B_MK, C_MK);
    yspline_MK[iev] = yspline2(Nmass_strange, temp_x, temp_y_MK, B_MK, C_MK, phys_ms);
    
  }// iev
  
  free(temp_x);
  free(temp_y_MK);
  free(B_MK);
  free(C_MK);
  
  ///////// FINE FACCIO LA SPLINE


  /////////// CHECK GRAFICO

  //////// check per l'evento iev_check_MK = 2
  int iev_check_MK = 2;
  
  FILE *fout_graph_check_MK;
  char file_out_graph_check_MK[LEN_NAME];
  
  sprintf(file_out_graph_check_MK, "OUTPUT_SMEAR/%s/%s/check_spline/MK/MK_spline_check.%s.%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V_2[ibeta].c_str(), mu_sea_2[imusea].c_str());
  
  if ((fout_graph_check_MK = fopen(file_out_graph_check_MK, "w")) == NULL ){
    printf("Error opening the output file file_out_graph_check_MK\n");
    exit(EXIT_FAILURE);
  }
  
  fprintf(fout_graph_check_MK, "@type xydy\n");
  
  for(int ims = i_ms0; ims < i_ms0 + Nmass_strange; ims++){
    
    //  prendendo le masse in GeV
    fprintf(fout_graph_check_MK, "%f %f %f\n", (mq[ibeta][ims]*ainv[ibeta][iev_check_MK])/Zev[ibeta][iev_check_MK], MK[ims-i_ms0][iev_check_MK], sigma_JK_modified_2(MK[ims-i_ms0], 1, 100, clusterfile));
    
    //  prendendo le masse in unità di reticolo
    //fprintf(fout_graph_check_MK, "%f %f %f\n", mq[ibeta][ims], MK[ims-i_ms0][iev_check_MK], sigma_JK_modified_2(MK[ims-i_ms0], 1, 100, clusterfile));
    
  }
  
  fprintf(fout_graph_check_MK, "&\n");
  
  //  prendendo le masse in GeV
  fprintf(fout_graph_check_MK, "%f %f %f\n", mstrange[iev_check_MK], yspline_MK[iev_check_MK], sigma_JK_modified_2(yspline_MK, 1, 100, clusterfile));
  
  //  prendendo le masse in unità di reticolo
  //fprintf(fout_graph_check_MK, "%f %f %f\n", (mstrange[iev_check_MK]/ainv[ibeta][iev_check_MK])*Zev[ibeta][iev_check_MK], yspline_MK[iev_check_MK], sigma_JK_modified_2(yspline_MK, 1, 100, clusterfile));
  
  fclose(fout_graph_check_MK);
  
  /////////// FINE CHECK GRAFICO
  

  ///////// OUTPUT

  char file_out_MK[LEN_NAME];
  FILE *fout_MK;
   
  sprintf(file_out_MK, "OUTPUT_SMEAR/%s/%s/MK/%s/%s/MK.%s.m1m2_%ss.th_3.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[0].c_str());
  
  if ((fout_MK = fopen(file_out_MK, "w")) == NULL ){
    printf("Error opening the output file MK!!\n");
    exit(EXIT_FAILURE);
  }
  
  for(int iev = 0; iev <= Nev; iev++ ){
    fprintf(fout_MK, "%+f\n", yspline_MK[iev]);
  }// iev
  
  fclose(fout_MK);

  free(yspline_MK);
  
  ///////// FINE OUTPUT

#endif
  
  ///////// FINE SPLINE MK


  ///////// SPLINE V0, Vi, S0, f0_S per il caso D_to_pi

#ifdef SPLINE_V_S_DPI
  
  ///////// CREO I BOOTSTRAP
  
  double V0_Dpi[Nmass_charm][Nev+1], Vi_Dpi[Nmass_charm][Nev+1], S0_Dpi[Nmass_charm][Nev+1], f0_S_Dpi[Nmass_charm][Nev+1];
  char file_open_V0_Dpi[LEN_NAME], file_open_Vi_Dpi[LEN_NAME], file_open_S0_Dpi[LEN_NAME], file_open_f0_S_Dpi[LEN_NAME];
  
  for(int imc = i_mc0; imc < i_mc0 + Nmass_charm; imc++){

    sprintf(file_open_V0_Dpi, "OUTPUT_SMEAR/%s/%s/V0/Hl_sym_av/%s/%s/V0_fit.%s.m1m2_%s%s.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[0].c_str(), m_[imc].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
    sprintf(file_open_Vi_Dpi, "OUTPUT_SMEAR/%s/%s/Vi/Hl_sym_av/%s/%s/Vi_fit.%s.m1m2_%s%s.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[0].c_str(), m_[imc].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
    sprintf(file_open_S0_Dpi, "OUTPUT_SMEAR/%s/%s/S0/Hl_sym_av/%s/%s/S0_fit.%s.m1m2_%s%s.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[0].c_str(), m_[imc].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
    sprintf(file_open_f0_S_Dpi, "OUTPUT_SMEAR/%s/%s/f0_S/Hl_sym_av/%s/%s/f0_S.%s.m1m2_%s%s.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[0].c_str(), m_[imc].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
    
    bootstrap_sampling( V0_Dpi[imc-i_mc0], file_open_V0_Dpi, ibeta, imusea, clusterfile);
    bootstrap_sampling( Vi_Dpi[imc-i_mc0], file_open_Vi_Dpi, ibeta, imusea, clusterfile);
    bootstrap_sampling( S0_Dpi[imc-i_mc0], file_open_S0_Dpi, ibeta, imusea, clusterfile);
    bootstrap_sampling( f0_S_Dpi[imc-i_mc0], file_open_f0_S_Dpi, ibeta, imusea, clusterfile);
    
  }// imc  
  
  ///////// FINE CREO I BOOTSTRAP

  ///////// FACCIO LA SPLINE

  temp_x = (double*)malloc(sizeof(double)*(Nmass_charm));
  double *temp_y_V0_Dpi = (double*)malloc(sizeof(double)*(Nmass_charm)), *temp_y_Vi_Dpi = (double*)malloc(sizeof(double)*(Nmass_charm));
  double *temp_y_S0_Dpi = (double*)malloc(sizeof(double)*(Nmass_charm)), *temp_y_f0_S_Dpi = (double*)malloc(sizeof(double)*(Nmass_charm));
  double *B_V0_Dpi = (double*)malloc(sizeof(double)*(Nmass_charm)), *B_Vi_Dpi = (double*)malloc(sizeof(double)*(Nmass_charm));
  double *B_S0_Dpi = (double*)malloc(sizeof(double)*(Nmass_charm)), *B_f0_S_Dpi = (double*)malloc(sizeof(double)*(Nmass_charm));
  double *C_V0_Dpi = (double*)malloc(sizeof(double)*(Nmass_charm)), *C_Vi_Dpi = (double*)malloc(sizeof(double)*(Nmass_charm));
  double *C_S0_Dpi = (double*)malloc(sizeof(double)*(Nmass_charm)), *C_f0_S_Dpi = (double*)malloc(sizeof(double)*(Nmass_charm));
  double *yspline_V0_Dpi = (double*)malloc(sizeof(double)*(Nev+1)), *yspline_Vi_Dpi = (double*)malloc(sizeof(double)*(Nev+1));
  double *yspline_S0_Dpi = (double*)malloc(sizeof(double)*(Nev+1)), *yspline_f0_S_Dpi = (double*)malloc(sizeof(double)*(Nev+1));
  
  for(int iev = 0; iev <= Nev; iev++){
    for(int imc = i_mc0; imc < i_mc0 + Nmass_charm; imc++){
      
      temp_x[imc-i_mc0] = (mq[ibeta][imc]*ainv[ibeta][0])/Zev[ibeta][0];  // prendo le masse in GeV
      //temp_x[imc-i_mc0] = mq[ibeta][imc];  // prendo le masse in unità di reticolo
      
      temp_y_V0_Dpi[imc-i_mc0] = V0_Dpi[imc-i_mc0][iev];
      temp_y_Vi_Dpi[imc-i_mc0] = Vi_Dpi[imc-i_mc0][iev];
      temp_y_S0_Dpi[imc-i_mc0] = S0_Dpi[imc-i_mc0][iev];
      temp_y_f0_S_Dpi[imc-i_mc0] = f0_S_Dpi[imc-i_mc0][iev];
      
    }// imc
    
    phys_mc = mcharm[iev];  // prendo le masse in GeV
    //phys_mc = (mcharm[iev]/ainv[ibeta][iev])*Zev[ibeta][iev];  // prendo le masse in unità di reticolo
    
    spline2(Nmass_charm, temp_x, temp_y_V0_Dpi, B_V0_Dpi, C_V0_Dpi);
    spline2(Nmass_charm, temp_x, temp_y_Vi_Dpi, B_Vi_Dpi, C_Vi_Dpi);
    spline2(Nmass_charm, temp_x, temp_y_S0_Dpi, B_S0_Dpi, C_S0_Dpi);
    spline2(Nmass_charm, temp_x, temp_y_f0_S_Dpi, B_f0_S_Dpi, C_f0_S_Dpi);
    
    yspline_V0_Dpi[iev] = yspline2(Nmass_charm, temp_x, temp_y_V0_Dpi, B_V0_Dpi, C_V0_Dpi, phys_mc);
    yspline_Vi_Dpi[iev] = yspline2(Nmass_charm, temp_x, temp_y_Vi_Dpi, B_Vi_Dpi, C_Vi_Dpi, phys_mc);
    yspline_S0_Dpi[iev] = yspline2(Nmass_charm, temp_x, temp_y_S0_Dpi, B_S0_Dpi, C_S0_Dpi, phys_mc);
    yspline_f0_S_Dpi[iev] = yspline2(Nmass_charm, temp_x, temp_y_f0_S_Dpi, B_f0_S_Dpi, C_f0_S_Dpi, phys_mc);
    
  }// iev
  
  free(temp_x);
  free(temp_y_V0_Dpi);
  free(temp_y_Vi_Dpi);
  free(temp_y_S0_Dpi);
  free(temp_y_f0_S_Dpi);
  free(B_V0_Dpi);
  free(B_Vi_Dpi);
  free(B_S0_Dpi);
  free(B_f0_S_Dpi);
  free(C_V0_Dpi);
  free(C_Vi_Dpi);
  free(C_S0_Dpi);
  free(C_f0_S_Dpi);
  
  ///////// FINE FACCIO LA SPLINE
  
  /////////// CHECK GRAFICO

  //////// check per l'evento iev_check_V_S_Dpi = 2
  int iev_check_V_S_Dpi = 2;
  
  FILE *fout_graph_check_V0_Dpi, *fout_graph_check_Vi_Dpi, *fout_graph_check_S0_Dpi, *fout_graph_check_f0_S_Dpi;
  char file_out_graph_check_V0_Dpi[LEN_NAME], file_out_graph_check_Vi_Dpi[LEN_NAME], file_out_graph_check_S0_Dpi[LEN_NAME], file_out_graph_check_f0_S_Dpi[LEN_NAME];
  
  sprintf(file_out_graph_check_V0_Dpi, "OUTPUT_SMEAR/%s/%s/check_spline/V0/V0_spline_check.%s.%s.m1m2_%sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V_2[ibeta].c_str(), mu_sea_2[imusea].c_str(), m_[0].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  sprintf(file_out_graph_check_Vi_Dpi, "OUTPUT_SMEAR/%s/%s/check_spline/Vi/Vi_spline_check.%s.%s.m1m2_%sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V_2[ibeta].c_str(), mu_sea_2[imusea].c_str(), m_[0].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  sprintf(file_out_graph_check_S0_Dpi, "OUTPUT_SMEAR/%s/%s/check_spline/S0/S0_spline_check.%s.%s.m1m2_%sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V_2[ibeta].c_str(), mu_sea_2[imusea].c_str(), m_[0].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  sprintf(file_out_graph_check_f0_S_Dpi, "OUTPUT_SMEAR/%s/%s/check_spline/f0_S/f0_S_spline_check.%s.%s.m1m2_%sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V_2[ibeta].c_str(), mu_sea_2[imusea].c_str(), m_[0].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  
  if ((fout_graph_check_V0_Dpi = fopen(file_out_graph_check_V0_Dpi, "w")) == NULL ){
    printf("Error opening the output file file_out_graph_check_V0_Dpi\n");
    exit(EXIT_FAILURE);
  }
  if ((fout_graph_check_Vi_Dpi = fopen(file_out_graph_check_Vi_Dpi, "w")) == NULL ){
    printf("Error opening the output file file_out_graph_check_Vi_Dpi\n");
    exit(EXIT_FAILURE);
  }
  if ((fout_graph_check_S0_Dpi = fopen(file_out_graph_check_S0_Dpi, "w")) == NULL ){
    printf("Error opening the output file file_out_graph_check_S0_Dpi\n");
    exit(EXIT_FAILURE);
  }
  if ((fout_graph_check_f0_S_Dpi = fopen(file_out_graph_check_f0_S_Dpi, "w")) == NULL ){
    printf("Error opening the output file file_out_graph_check_f0_S_Dpi\n");
    exit(EXIT_FAILURE);
  }
  
  fprintf(fout_graph_check_V0_Dpi, "@type xydy\n");
  fprintf(fout_graph_check_Vi_Dpi, "@type xydy\n");
  fprintf(fout_graph_check_S0_Dpi, "@type xydy\n");
  fprintf(fout_graph_check_f0_S_Dpi, "@type xydy\n");
  
  for(int imc = i_mc0; imc < i_mc0 + Nmass_charm; imc++){
    
    //  prendendo le masse in GeV
    fprintf(fout_graph_check_V0_Dpi, "%f %f %f\n", (mq[ibeta][imc]*ainv[ibeta][iev_check_V_S_Dpi])/Zev[ibeta][iev_check_V_S_Dpi], V0_Dpi[imc-i_mc0][iev_check_V_S_Dpi], sigma_JK_modified_2(V0_Dpi[imc-i_mc0], 1, 100, clusterfile));
    fprintf(fout_graph_check_Vi_Dpi, "%f %f %f\n", (mq[ibeta][imc]*ainv[ibeta][iev_check_V_S_Dpi])/Zev[ibeta][iev_check_V_S_Dpi], Vi_Dpi[imc-i_mc0][iev_check_V_S_Dpi], sigma_JK_modified_2(Vi_Dpi[imc-i_mc0], 1, 100, clusterfile));
    fprintf(fout_graph_check_S0_Dpi, "%f %f %f\n", (mq[ibeta][imc]*ainv[ibeta][iev_check_V_S_Dpi])/Zev[ibeta][iev_check_V_S_Dpi], S0_Dpi[imc-i_mc0][iev_check_V_S_Dpi], sigma_JK_modified_2(S0_Dpi[imc-i_mc0], 1, 100, clusterfile));
    fprintf(fout_graph_check_f0_S_Dpi, "%f %f %f\n", (mq[ibeta][imc]*ainv[ibeta][iev_check_V_S_Dpi])/Zev[ibeta][iev_check_V_S_Dpi], f0_S_Dpi[imc-i_mc0][iev_check_V_S_Dpi], sigma_JK_modified_2(f0_S_Dpi[imc-i_mc0], 1, 100, clusterfile));
    
    //  prendendo le masse in unità di reticolo
    //fprintf(fout_graph_check_V0_Dpi, "%f %f %f\n", mq[ibeta][imc], V0_Dpi[imc-i_mc0][iev_check_V_S_Dpi], sigma_JK_modified_2(V0_Dpi[imc-i_mc0], 1, 100, clusterfile));
    //fprintf(fout_graph_check_Vi_Dpi, "%f %f %f\n", mq[ibeta][imc], Vi_Dpi[imc-i_mc0][iev_check_V_S_Dpi], sigma_JK_modified_2(Vi_Dpi[imc-i_mc0], 1, 100, clusterfile));
    //fprintf(fout_graph_check_S0_Dpi, "%f %f %f\n", mq[ibeta][imc], S0_Dpi[imc-i_mc0][iev_check_V_S_Dpi], sigma_JK_modified_2(S0_Dpi[imc-i_mc0], 1, 100, clusterfile));
    //fprintf(fout_graph_check_f0_S_Dpi, "%f %f %f\n", mq[ibeta][imc], f0_S_Dpi[imc-i_mc0][iev_check_V_S_Dpi], sigma_JK_modified_2(f0_S_Dpi[imc-i_mc0], 1, 100, clusterfile));
    
  }// imc
  
  fprintf(fout_graph_check_V0_Dpi, "&\n");
  fprintf(fout_graph_check_Vi_Dpi, "&\n");
  fprintf(fout_graph_check_S0_Dpi, "&\n");
  fprintf(fout_graph_check_f0_S_Dpi, "&\n");
  
  
  //  prendendo le masse in GeV
  fprintf(fout_graph_check_V0_Dpi, "%f %f %f\n", mcharm[iev_check_V_S_Dpi], yspline_V0_Dpi[iev_check_V_S_Dpi], sigma_JK_modified_2(yspline_V0_Dpi, 1, 100, clusterfile));
  fprintf(fout_graph_check_Vi_Dpi, "%f %f %f\n", mcharm[iev_check_V_S_Dpi], yspline_Vi_Dpi[iev_check_V_S_Dpi], sigma_JK_modified_2(yspline_Vi_Dpi, 1, 100, clusterfile));
  fprintf(fout_graph_check_S0_Dpi, "%f %f %f\n", mcharm[iev_check_V_S_Dpi], yspline_S0_Dpi[iev_check_V_S_Dpi], sigma_JK_modified_2(yspline_S0_Dpi, 1, 100, clusterfile));
  fprintf(fout_graph_check_f0_S_Dpi, "%f %f %f\n", mcharm[iev_check_V_S_Dpi], yspline_f0_S_Dpi[iev_check_V_S_Dpi], sigma_JK_modified_2(yspline_f0_S_Dpi, 1, 100, clusterfile));
  
  //  prendendo le masse in unità di reticolo
  //fprintf(fout_graph_check_V0_Dpi, "%f %f %f\n", (mcharm[iev_check_V_S_Dpi]/ainv[ibeta][iev_check_V_S_Dpi])*Zev[ibeta][iev_check_V_S_Dpi], yspline_V0_Dpi[iev_check_V_S_Dpi], sigma_JK_modified_2(yspline_V0_Dpi, 1, 100, clusterfile));
  //fprintf(fout_graph_check_Vi_Dpi, "%f %f %f\n", (mcharm[iev_check_V_S_Dpi]/ainv[ibeta][iev_check_V_S_Dpi])*Zev[ibeta][iev_check_V_S_Dpi], yspline_Vi_Dpi[iev_check_V_S_Dpi], sigma_JK_modified_2(yspline_Vi_Dpi, 1, 100, clusterfile));
  //fprintf(fout_graph_check_S0_Dpi, "%f %f %f\n", (mcharm[iev_check_V_S_Dpi]/ainv[ibeta][iev_check_V_S_Dpi])*Zev[ibeta][iev_check_V_S_Dpi], yspline_S0_Dpi[iev_check_V_S_Dpi], sigma_JK_modified_2(yspline_S0_Dpi, 1, 100, clusterfile));
  //fprintf(fout_graph_check_f0_S_Dpi, "%f %f %f\n", (mcharm[iev_check_V_S_Dpi]/ainv[ibeta][iev_check_V_S_Dpi])*Zev[ibeta][iev_check_V_S_Dpi], yspline_f0_S_Dpi[iev_check_V_S_Dpi], sigma_JK_modified_2(yspline_f0_S_Dpi, 1, 100, clusterfile));
  
  fclose(fout_graph_check_V0_Dpi);
  fclose(fout_graph_check_Vi_Dpi);
  fclose(fout_graph_check_S0_Dpi);
  fclose(fout_graph_check_f0_S_Dpi);
  
  /////////// FINE CHECK GRAFICO  
  

   /////////// OUTPUT

  char file_out_Vi_Dpi[LEN_NAME], file_out_V0_Dpi[LEN_NAME], file_out_S0_Dpi[LEN_NAME], file_out_f0_S_Dpi[LEN_NAME];
  FILE *fout_Vi_Dpi, *fout_V0_Dpi, *fout_S0_Dpi, *fout_f0_S_Dpi;

  sprintf(file_out_Vi_Dpi, "OUTPUT_SMEAR/%s/%s/Vi/Hl_sym_av/%s/%s/Vi.%s.m1m2_%sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[0].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  sprintf(file_out_V0_Dpi, "OUTPUT_SMEAR/%s/%s/V0/Hl_sym_av/%s/%s/V0.%s.m1m2_%sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[0].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  sprintf(file_out_S0_Dpi, "OUTPUT_SMEAR/%s/%s/S0/Hl_sym_av/%s/%s/S0.%s.m1m2_%sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[0].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  sprintf(file_out_f0_S_Dpi, "OUTPUT_SMEAR/%s/%s/f0_S/Hl_sym_av/%s/%s/f0_S.%s.m1m2_%sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[0].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  
  if ((fout_Vi_Dpi = fopen(file_out_Vi_Dpi, "w")) == NULL ){
    printf("Error opening the output file Vi_Dpi!!\n");
    exit(EXIT_FAILURE);
  }
  if ((fout_V0_Dpi = fopen(file_out_V0_Dpi, "w")) == NULL ){
    printf("Error opening the output file V0_Dpi!!\n");
    exit(EXIT_FAILURE);
  }
  if ((fout_S0_Dpi = fopen(file_out_S0_Dpi, "w")) == NULL ){
    printf("Error opening the output file S0_Dpi!!\n");
    exit(EXIT_FAILURE);
  }
  if ((fout_f0_S_Dpi = fopen(file_out_f0_S_Dpi, "w")) == NULL ){
    printf("Error opening the output file f0_S_Dpi!!\n");
    exit(EXIT_FAILURE);
  }
  
  
  for(int iev = 0; iev <= Nev; iev++ ){
    
    fprintf(fout_V0_Dpi, "%+f\n", yspline_V0_Dpi[iev]);
    fprintf(fout_S0_Dpi, "%+f\n", yspline_S0_Dpi[iev]);
    fprintf(fout_f0_S_Dpi, "%+f\n", yspline_f0_S_Dpi[iev]);
    fprintf(fout_Vi_Dpi, "%+f\n", yspline_Vi_Dpi[iev]);
    
  }// iev

  free(yspline_V0_Dpi);
  free(yspline_Vi_Dpi);
  free(yspline_S0_Dpi);
  free(yspline_f0_S_Dpi);
  
  fclose(fout_V0_Dpi);
  fclose(fout_S0_Dpi);
  fclose(fout_f0_S_Dpi);
  fclose(fout_Vi_Dpi);
  
  /////////// FINE OUTPUT
  
#endif

  ///////// FINE SPLINE V0, Vi, S0, f0_S per il caso D_to_pi

  

  ///////// SPLINE V0, Vi, S0, f0_S per il caso D_to_K
  
#ifdef SPLINE_V_S_DK
  
  ///////// SPLINE V0, Vi, S0, f0_S per il caso D_to_K_charm
  
#ifdef SPLINE_V_S_DK_charm
  
  ///////// CREO I BOOTSTRAP
  
  double V0_DK_charm[Nmass_charm][Nev+1], Vi_DK_charm[Nmass_charm][Nev+1], S0_DK_charm[Nmass_charm][Nev+1], f0_S_DK_charm[Nmass_charm][Nev+1];
  char file_open_V0_DK_charm[LEN_NAME], file_open_Vi_DK_charm[LEN_NAME], file_open_S0_DK_charm[LEN_NAME], file_open_f0_S_DK_charm[LEN_NAME];
  
  for(int imc = i_mc0; imc < i_mc0 + Nmass_charm; imc++){

    sprintf(file_open_V0_DK_charm, "OUTPUT_SMEAR/%s/%s/V0/Hl_sym_av/%s/%s/V0_fit.%s.m1m2_%s%s.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im].c_str(), m_[imc].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
    sprintf(file_open_Vi_DK_charm, "OUTPUT_SMEAR/%s/%s/Vi/Hl_sym_av/%s/%s/Vi_fit.%s.m1m2_%s%s.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im].c_str(), m_[imc].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
    sprintf(file_open_S0_DK_charm, "OUTPUT_SMEAR/%s/%s/S0/Hl_sym_av/%s/%s/S0_fit.%s.m1m2_%s%s.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im].c_str(), m_[imc].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
    sprintf(file_open_f0_S_DK_charm, "OUTPUT_SMEAR/%s/%s/f0_S/Hl_sym_av/%s/%s/f0_S.%s.m1m2_%s%s.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im].c_str(), m_[imc].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
    
    bootstrap_sampling( V0_DK_charm[imc-i_mc0], file_open_V0_DK_charm, ibeta, imusea, clusterfile);
    bootstrap_sampling( Vi_DK_charm[imc-i_mc0], file_open_Vi_DK_charm, ibeta, imusea, clusterfile);
    bootstrap_sampling( S0_DK_charm[imc-i_mc0], file_open_S0_DK_charm, ibeta, imusea, clusterfile);
    bootstrap_sampling( f0_S_DK_charm[imc-i_mc0], file_open_f0_S_DK_charm, ibeta, imusea, clusterfile);
    
  }// imc  
  
  ///////// FINE CREO I BOOTSTRAP

  ///////// FACCIO LA SPLINE

  temp_x = (double*)malloc(sizeof(double)*(Nmass_charm));
  double *temp_y_V0_DK_charm = (double*)malloc(sizeof(double)*(Nmass_charm)), *temp_y_Vi_DK_charm = (double*)malloc(sizeof(double)*(Nmass_charm));
  double *temp_y_S0_DK_charm = (double*)malloc(sizeof(double)*(Nmass_charm)), *temp_y_f0_S_DK_charm = (double*)malloc(sizeof(double)*(Nmass_charm));
  double *B_V0_DK_charm = (double*)malloc(sizeof(double)*(Nmass_charm)), *B_Vi_DK_charm = (double*)malloc(sizeof(double)*(Nmass_charm));
  double *B_S0_DK_charm = (double*)malloc(sizeof(double)*(Nmass_charm)), *B_f0_S_DK_charm = (double*)malloc(sizeof(double)*(Nmass_charm));
  double *C_V0_DK_charm = (double*)malloc(sizeof(double)*(Nmass_charm)), *C_Vi_DK_charm = (double*)malloc(sizeof(double)*(Nmass_charm));
  double *C_S0_DK_charm = (double*)malloc(sizeof(double)*(Nmass_charm)), *C_f0_S_DK_charm = (double*)malloc(sizeof(double)*(Nmass_charm));
  double *yspline_V0_DK_charm = (double*)malloc(sizeof(double)*(Nev+1)), *yspline_Vi_DK_charm = (double*)malloc(sizeof(double)*(Nev+1));
  double *yspline_S0_DK_charm = (double*)malloc(sizeof(double)*(Nev+1)), *yspline_f0_S_DK_charm = (double*)malloc(sizeof(double)*(Nev+1));
  
  for(int iev = 0; iev <= Nev; iev++){
    for(int imc = i_mc0; imc < i_mc0 + Nmass_charm; imc++){
      
      temp_x[imc-i_mc0] = (mq[ibeta][imc]*ainv[ibeta][0])/Zev[ibeta][0];  // prendo le masse in GeV
      //temp_x[imc-i_mc0] = mq[ibeta][imc];  // prendo le masse in unità di reticolo
      
      temp_y_V0_DK_charm[imc-i_mc0] = V0_DK_charm[imc-i_mc0][iev];
      temp_y_Vi_DK_charm[imc-i_mc0] = Vi_DK_charm[imc-i_mc0][iev];
      temp_y_S0_DK_charm[imc-i_mc0] = S0_DK_charm[imc-i_mc0][iev];
      temp_y_f0_S_DK_charm[imc-i_mc0] = f0_S_DK_charm[imc-i_mc0][iev];
      
    }// imc
    
    phys_mc = mcharm[iev];  // prendo le masse in GeV
    //phys_mc = (mcharm[iev]/ainv[ibeta][iev])*Zev[ibeta][iev];  // prendo le masse in unità di reticolo
    
    spline2(Nmass_charm, temp_x, temp_y_V0_DK_charm, B_V0_DK_charm, C_V0_DK_charm);
    spline2(Nmass_charm, temp_x, temp_y_Vi_DK_charm, B_Vi_DK_charm, C_Vi_DK_charm);
    spline2(Nmass_charm, temp_x, temp_y_S0_DK_charm, B_S0_DK_charm, C_S0_DK_charm);
    spline2(Nmass_charm, temp_x, temp_y_f0_S_DK_charm, B_f0_S_DK_charm, C_f0_S_DK_charm);
    
    yspline_V0_DK_charm[iev] = yspline2(Nmass_charm, temp_x, temp_y_V0_DK_charm, B_V0_DK_charm, C_V0_DK_charm, phys_mc);
    yspline_Vi_DK_charm[iev] = yspline2(Nmass_charm, temp_x, temp_y_Vi_DK_charm, B_Vi_DK_charm, C_Vi_DK_charm, phys_mc);
    yspline_S0_DK_charm[iev] = yspline2(Nmass_charm, temp_x, temp_y_S0_DK_charm, B_S0_DK_charm, C_S0_DK_charm, phys_mc);
    yspline_f0_S_DK_charm[iev] = yspline2(Nmass_charm, temp_x, temp_y_f0_S_DK_charm, B_f0_S_DK_charm, C_f0_S_DK_charm, phys_mc);
    
  }// iev
  
  free(temp_x);
  free(temp_y_V0_DK_charm);
  free(temp_y_Vi_DK_charm);
  free(temp_y_S0_DK_charm);
  free(temp_y_f0_S_DK_charm);
  free(B_V0_DK_charm);
  free(B_Vi_DK_charm);
  free(B_S0_DK_charm);
  free(B_f0_S_DK_charm);
  free(C_V0_DK_charm);
  free(C_Vi_DK_charm);
  free(C_S0_DK_charm);
  free(C_f0_S_DK_charm);
  
  ///////// FINE FACCIO LA SPLINE
  
  /////////// CHECK GRAFICO

  //////// check per l'evento iev_check_V_S_DK_charm = 2
  int iev_check_V_S_DK_charm = 2;
  
  FILE *fout_graph_check_V0_DK_charm, *fout_graph_check_Vi_DK_charm, *fout_graph_check_S0_DK_charm, *fout_graph_check_f0_S_DK_charm;
  char file_out_graph_check_V0_DK_charm[LEN_NAME], file_out_graph_check_Vi_DK_charm[LEN_NAME], file_out_graph_check_S0_DK_charm[LEN_NAME], file_out_graph_check_f0_S_DK_charm[LEN_NAME];
  
  sprintf(file_out_graph_check_V0_DK_charm, "OUTPUT_SMEAR/%s/%s/check_spline/V0/V0_spline_check.%s.%s.m1m2_%sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V_2[ibeta].c_str(), mu_sea_2[imusea].c_str(), m_[im].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  sprintf(file_out_graph_check_Vi_DK_charm, "OUTPUT_SMEAR/%s/%s/check_spline/Vi/Vi_spline_check.%s.%s.m1m2_%sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V_2[ibeta].c_str(), mu_sea_2[imusea].c_str(), m_[im].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  sprintf(file_out_graph_check_S0_DK_charm, "OUTPUT_SMEAR/%s/%s/check_spline/S0/S0_spline_check.%s.%s.m1m2_%sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V_2[ibeta].c_str(), mu_sea_2[imusea].c_str(), m_[im].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  sprintf(file_out_graph_check_f0_S_DK_charm, "OUTPUT_SMEAR/%s/%s/check_spline/f0_S/f0_S_spline_check.%s.%s.m1m2_%sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V_2[ibeta].c_str(), mu_sea_2[imusea].c_str(), m_[im].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  
  if ((fout_graph_check_V0_DK_charm = fopen(file_out_graph_check_V0_DK_charm, "w")) == NULL ){
    printf("Error opening the output file file_out_graph_check_V0_DK_charm\n");
    exit(EXIT_FAILURE);
  }
  if ((fout_graph_check_Vi_DK_charm = fopen(file_out_graph_check_Vi_DK_charm, "w")) == NULL ){
    printf("Error opening the output file file_out_graph_check_Vi_DK_charm\n");
    exit(EXIT_FAILURE);
  }
  if ((fout_graph_check_S0_DK_charm = fopen(file_out_graph_check_S0_DK_charm, "w")) == NULL ){
    printf("Error opening the output file file_out_graph_check_S0_DK_charm\n");
    exit(EXIT_FAILURE);
  }
  if ((fout_graph_check_f0_S_DK_charm = fopen(file_out_graph_check_f0_S_DK_charm, "w")) == NULL ){
    printf("Error opening the output file file_out_graph_check_f0_S_DK_charm\n");
    exit(EXIT_FAILURE);
  }
  
  fprintf(fout_graph_check_V0_DK_charm, "@type xydy\n");
  fprintf(fout_graph_check_Vi_DK_charm, "@type xydy\n");
  fprintf(fout_graph_check_S0_DK_charm, "@type xydy\n");
  fprintf(fout_graph_check_f0_S_DK_charm, "@type xydy\n");
  
  for(int imc = i_mc0; imc < i_mc0 + Nmass_charm; imc++){
    
    //  prendendo le masse in GeV
    fprintf(fout_graph_check_V0_DK_charm, "%f %f %f\n", (mq[ibeta][imc]*ainv[ibeta][iev_check_V_S_DK_charm])/Zev[ibeta][iev_check_V_S_DK_charm], V0_DK_charm[imc-i_mc0][iev_check_V_S_DK_charm], sigma_JK_modified_2(V0_DK_charm[imc-i_mc0], 1, 100, clusterfile));
    fprintf(fout_graph_check_Vi_DK_charm, "%f %f %f\n", (mq[ibeta][imc]*ainv[ibeta][iev_check_V_S_DK_charm])/Zev[ibeta][iev_check_V_S_DK_charm], Vi_DK_charm[imc-i_mc0][iev_check_V_S_DK_charm], sigma_JK_modified_2(Vi_DK_charm[imc-i_mc0], 1, 100, clusterfile));
    fprintf(fout_graph_check_S0_DK_charm, "%f %f %f\n", (mq[ibeta][imc]*ainv[ibeta][iev_check_V_S_DK_charm])/Zev[ibeta][iev_check_V_S_DK_charm], S0_DK_charm[imc-i_mc0][iev_check_V_S_DK_charm], sigma_JK_modified_2(S0_DK_charm[imc-i_mc0], 1, 100, clusterfile));
    fprintf(fout_graph_check_f0_S_DK_charm, "%f %f %f\n", (mq[ibeta][imc]*ainv[ibeta][iev_check_V_S_DK_charm])/Zev[ibeta][iev_check_V_S_DK_charm], f0_S_DK_charm[imc-i_mc0][iev_check_V_S_DK_charm], sigma_JK_modified_2(f0_S_DK_charm[imc-i_mc0], 1, 100, clusterfile));
    
    //  prendendo le masse in unità di reticolo
    //fprintf(fout_graph_check_V0_DK_charm, "%f %f %f\n", mq[ibeta][imc], V0_DK_charm[imc-i_mc0][iev_check_V_S_DK_charm], sigma_JK_modified_2(V0_DK_charm[imc-i_mc0], 1, 100, clusterfile));
    //fprintf(fout_graph_check_Vi_DK_charm, "%f %f %f\n", mq[ibeta][imc], Vi_DK_charm[imc-i_mc0][iev_check_V_S_DK_charm], sigma_JK_modified_2(Vi_DK_charm[imc-i_mc0], 1, 100, clusterfile));
    //fprintf(fout_graph_check_S0_DK_charm, "%f %f %f\n", mq[ibeta][imc], S0_DK_charm[imc-i_mc0][iev_check_V_S_DK_charm], sigma_JK_modified_2(S0_DK_charm[imc-i_mc0], 1, 100, clusterfile));
    //fprintf(fout_graph_check_f0_S_DK_charm, "%f %f %f\n", mq[ibeta][imc], f0_S_DK_charm[imc-i_mc0][iev_check_V_S_DK_charm], sigma_JK_modified_2(f0_S_DK_charm[imc-i_mc0], 1, 100, clusterfile));
    
  }// imc
  
  fprintf(fout_graph_check_V0_DK_charm, "&\n");
  fprintf(fout_graph_check_Vi_DK_charm, "&\n");
  fprintf(fout_graph_check_S0_DK_charm, "&\n");
  fprintf(fout_graph_check_f0_S_DK_charm, "&\n");
  
  
  //  prendendo le masse in GeV
  fprintf(fout_graph_check_V0_DK_charm, "%f %f %f\n", mcharm[iev_check_V_S_DK_charm], yspline_V0_DK_charm[iev_check_V_S_DK_charm], sigma_JK_modified_2(yspline_V0_DK_charm, 1, 100, clusterfile));
  fprintf(fout_graph_check_Vi_DK_charm, "%f %f %f\n", mcharm[iev_check_V_S_DK_charm], yspline_Vi_DK_charm[iev_check_V_S_DK_charm], sigma_JK_modified_2(yspline_Vi_DK_charm, 1, 100, clusterfile));
  fprintf(fout_graph_check_S0_DK_charm, "%f %f %f\n", mcharm[iev_check_V_S_DK_charm], yspline_S0_DK_charm[iev_check_V_S_DK_charm], sigma_JK_modified_2(yspline_S0_DK_charm, 1, 100, clusterfile));
  fprintf(fout_graph_check_f0_S_DK_charm, "%f %f %f\n", mcharm[iev_check_V_S_DK_charm], yspline_f0_S_DK_charm[iev_check_V_S_DK_charm], sigma_JK_modified_2(yspline_f0_S_DK_charm, 1, 100, clusterfile));
  
  //  prendendo le masse in unità di reticolo
  //fprintf(fout_graph_check_V0_DK_charm, "%f %f %f\n", (mcharm[iev_check_V_S_DK_charm]/ainv[ibeta][iev_check_V_S_DK_charm])*Zev[ibeta][iev_check_V_S_DK_charm], yspline_V0_DK_charm[iev_check_V_S_DK_charm], sigma_JK_modified_2(yspline_V0_DK_charm, 1, 100, clusterfile));
  //fprintf(fout_graph_check_Vi_DK_charm, "%f %f %f\n", (mcharm[iev_check_V_S_DK_charm]/ainv[ibeta][iev_check_V_S_DK_charm])*Zev[ibeta][iev_check_V_S_DK_charm], yspline_Vi_DK_charm[iev_check_V_S_DK_charm], sigma_JK_modified_2(yspline_Vi_DK_charm, 1, 100, clusterfile));
  //fprintf(fout_graph_check_S0_DK_charm, "%f %f %f\n", (mcharm[iev_check_V_S_DK_charm]/ainv[ibeta][iev_check_V_S_DK_charm])*Zev[ibeta][iev_check_V_S_DK_charm], yspline_S0_DK_charm[iev_check_V_S_DK_charm], sigma_JK_modified_2(yspline_S0_DK_charm, 1, 100, clusterfile));
  //fprintf(fout_graph_check_f0_S_DK_charm, "%f %f %f\n", (mcharm[iev_check_V_S_DK_charm]/ainv[ibeta][iev_check_V_S_DK_charm])*Zev[ibeta][iev_check_V_S_DK_charm], yspline_f0_S_DK_charm[iev_check_V_S_DK_charm], sigma_JK_modified_2(yspline_f0_S_DK_charm, 1, 100, clusterfile));
  
  fclose(fout_graph_check_V0_DK_charm);
  fclose(fout_graph_check_Vi_DK_charm);
  fclose(fout_graph_check_S0_DK_charm);
  fclose(fout_graph_check_f0_S_DK_charm);
  
  /////////// FINE CHECK GRAFICO  
  
  /////////// OUTPUT

  char file_out_Vi_DK_charm[LEN_NAME], file_out_V0_DK_charm[LEN_NAME], file_out_S0_DK_charm[LEN_NAME], file_out_f0_S_DK_charm[LEN_NAME];
  FILE *fout_Vi_DK_charm, *fout_V0_DK_charm, *fout_S0_DK_charm, *fout_f0_S_DK_charm;

  sprintf(file_out_Vi_DK_charm, "OUTPUT_SMEAR/%s/%s/Vi/Hl_sym_av/%s/%s/Vi.%s.m1m2_%sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  sprintf(file_out_V0_DK_charm, "OUTPUT_SMEAR/%s/%s/V0/Hl_sym_av/%s/%s/V0.%s.m1m2_%sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  sprintf(file_out_S0_DK_charm, "OUTPUT_SMEAR/%s/%s/S0/Hl_sym_av/%s/%s/S0.%s.m1m2_%sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  sprintf(file_out_f0_S_DK_charm, "OUTPUT_SMEAR/%s/%s/f0_S/Hl_sym_av/%s/%s/f0_S.%s.m1m2_%sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  
  if ((fout_Vi_DK_charm = fopen(file_out_Vi_DK_charm, "w")) == NULL ){
    printf("Error opening the output file Vi_DK_charm!!\n");
    exit(EXIT_FAILURE);
  }
  if ((fout_V0_DK_charm = fopen(file_out_V0_DK_charm, "w")) == NULL ){
    printf("Error opening the output file V0_DK_charm!!\n");
    exit(EXIT_FAILURE);
  }
  if ((fout_S0_DK_charm = fopen(file_out_S0_DK_charm, "w")) == NULL ){
    printf("Error opening the output file S0_DK_charm!!\n");
    exit(EXIT_FAILURE);
  }
  if ((fout_f0_S_DK_charm = fopen(file_out_f0_S_DK_charm, "w")) == NULL ){
    printf("Error opening the output file f0_S_DK_charm!!\n");
    exit(EXIT_FAILURE);
  }
  
  
  for(int iev = 0; iev <= Nev; iev++ ){
    
    fprintf(fout_V0_DK_charm, "%+f\n", yspline_V0_DK_charm[iev]);
    fprintf(fout_S0_DK_charm, "%+f\n", yspline_S0_DK_charm[iev]);
    fprintf(fout_f0_S_DK_charm, "%+f\n", yspline_f0_S_DK_charm[iev]);
    fprintf(fout_Vi_DK_charm, "%+f\n", yspline_Vi_DK_charm[iev]);
    
  }// iev

  free(yspline_V0_DK_charm);
  free(yspline_Vi_DK_charm);
  free(yspline_S0_DK_charm);
  free(yspline_f0_S_DK_charm);
  
  fclose(fout_V0_DK_charm);
  fclose(fout_S0_DK_charm);
  fclose(fout_f0_S_DK_charm);
  fclose(fout_Vi_DK_charm);
  
  /////////// FINE OUTPUT
  
#endif

  ///////// FINE SPLINE V0, Vi, S0, f0_S per il caso D_to_K_charm


  ///////// SPLINE V0, Vi, S0, f0_S per il caso D_to_K_strange
  
#ifdef SPLINE_V_S_DK_strange
  
  ///////// CREO I BOOTSTRAP
  
  double V0_DK_strange[Nmass_strange][Nev+1], Vi_DK_strange[Nmass_strange][Nev+1], S0_DK_strange[Nmass_strange][Nev+1], f0_S_DK_strange[Nmass_strange][Nev+1];
  char file_open_V0_DK_strange[LEN_NAME], file_open_Vi_DK_strange[LEN_NAME], file_open_S0_DK_strange[LEN_NAME], file_open_f0_S_DK_strange[LEN_NAME];
  
  for(int ims = i_ms0; ims < i_ms0 + Nmass_strange; ims++){
    
    sprintf(file_open_V0_DK_strange, "OUTPUT_SMEAR/%s/%s/V0/Hl_sym_av/%s/%s/V0.%s.m1m2_%sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[ims].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
    sprintf(file_open_Vi_DK_strange, "OUTPUT_SMEAR/%s/%s/Vi/Hl_sym_av/%s/%s/Vi.%s.m1m2_%sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[ims].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
    sprintf(file_open_S0_DK_strange, "OUTPUT_SMEAR/%s/%s/S0/Hl_sym_av/%s/%s/S0.%s.m1m2_%sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[ims].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
    sprintf(file_open_f0_S_DK_strange, "OUTPUT_SMEAR/%s/%s/f0_S/Hl_sym_av/%s/%s/f0_S.%s.m1m2_%sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[ims].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
    
    bootstrap_sampling( V0_DK_strange[ims-i_ms0], file_open_V0_DK_strange, ibeta, imusea, clusterfile);
    bootstrap_sampling( Vi_DK_strange[ims-i_ms0], file_open_Vi_DK_strange, ibeta, imusea, clusterfile);
    bootstrap_sampling( S0_DK_strange[ims-i_ms0], file_open_S0_DK_strange, ibeta, imusea, clusterfile);
    bootstrap_sampling( f0_S_DK_strange[ims-i_ms0], file_open_f0_S_DK_strange, ibeta, imusea, clusterfile);
    
  }// ims  
  
  ///////// FINE CREO I BOOTSTRAP
  
  
  
  ///////// FACCIO LA SPLINE
  
  temp_x = (double*)malloc(sizeof(double)*(Nmass_strange));
  double *temp_y_V0_DK_strange = (double*)malloc(sizeof(double)*(Nmass_strange)), *temp_y_Vi_DK_strange = (double*)malloc(sizeof(double)*(Nmass_strange));
  double *temp_y_S0_DK_strange = (double*)malloc(sizeof(double)*(Nmass_strange)), *temp_y_f0_S_DK_strange = (double*)malloc(sizeof(double)*(Nmass_strange));
  double *B_V0_DK_strange = (double*)malloc(sizeof(double)*(Nmass_strange)), *B_Vi_DK_strange = (double*)malloc(sizeof(double)*(Nmass_strange));
  double *B_S0_DK_strange = (double*)malloc(sizeof(double)*(Nmass_strange)), *B_f0_S_DK_strange = (double*)malloc(sizeof(double)*(Nmass_strange));
  double *C_V0_DK_strange = (double*)malloc(sizeof(double)*(Nmass_strange)), *C_Vi_DK_strange = (double*)malloc(sizeof(double)*(Nmass_strange));
  double *C_S0_DK_strange = (double*)malloc(sizeof(double)*(Nmass_strange)), *C_f0_S_DK_strange = (double*)malloc(sizeof(double)*(Nmass_strange));
  double *yspline_V0_DK_strange = (double*)malloc(sizeof(double)*(Nev+1)), *yspline_Vi_DK_strange = (double*)malloc(sizeof(double)*(Nev+1));
  double *yspline_S0_DK_strange = (double*)malloc(sizeof(double)*(Nev+1)), *yspline_f0_S_DK_strange = (double*)malloc(sizeof(double)*(Nev+1));
  
  for(int iev = 0; iev <= Nev; iev++){
    for(int ims = i_ms0; ims < i_ms0 + Nmass_strange; ims++){
      
      temp_x[ims-i_ms0] = (mq[ibeta][ims]*ainv[ibeta][0])/Zev[ibeta][0];  // prendo le masse in GeV
      //temp_x[ims-i_ms0] = mq[ibeta][ims];  // prendo le masse in unità di reticolo
      
      temp_y_V0_DK_strange[ims-i_ms0] = V0_DK_strange[ims-i_ms0][iev];
      temp_y_Vi_DK_strange[ims-i_ms0] = Vi_DK_strange[ims-i_ms0][iev];
      temp_y_S0_DK_strange[ims-i_ms0] = S0_DK_strange[ims-i_ms0][iev];
      temp_y_f0_S_DK_strange[ims-i_ms0] = f0_S_DK_strange[ims-i_ms0][iev];
      
    }// ims
    
    phys_ms = mstrange[iev];  // prendo le masse in GeV
    //phys_ms = (mstrange[iev]/ainv[ibeta][iev])*Zev[ibeta][iev];  // prendo le masse in unità di reticolo
    
    spline2(Nmass_strange, temp_x, temp_y_V0_DK_strange, B_V0_DK_strange, C_V0_DK_strange);
    spline2(Nmass_strange, temp_x, temp_y_Vi_DK_strange, B_Vi_DK_strange, C_Vi_DK_strange);
    spline2(Nmass_strange, temp_x, temp_y_S0_DK_strange, B_S0_DK_strange, C_S0_DK_strange);
    spline2(Nmass_strange, temp_x, temp_y_f0_S_DK_strange, B_f0_S_DK_strange, C_f0_S_DK_strange);
    
    yspline_V0_DK_strange[iev] = yspline2(Nmass_strange, temp_x, temp_y_V0_DK_strange, B_V0_DK_strange, C_V0_DK_strange, phys_ms);
    yspline_Vi_DK_strange[iev] = yspline2(Nmass_strange, temp_x, temp_y_Vi_DK_strange, B_Vi_DK_strange, C_Vi_DK_strange, phys_ms);
    yspline_S0_DK_strange[iev] = yspline2(Nmass_strange, temp_x, temp_y_S0_DK_strange, B_S0_DK_strange, C_S0_DK_strange, phys_ms);
    yspline_f0_S_DK_strange[iev] = yspline2(Nmass_strange, temp_x, temp_y_f0_S_DK_strange, B_f0_S_DK_strange, C_f0_S_DK_strange, phys_ms);
    
  }// iev
  
  free(temp_x);
  free(temp_y_V0_DK_strange);
  free(temp_y_Vi_DK_strange);
  free(temp_y_S0_DK_strange);
  free(temp_y_f0_S_DK_strange);
  free(B_V0_DK_strange);
  free(B_Vi_DK_strange);
  free(B_S0_DK_strange);
  free(B_f0_S_DK_strange);
  free(C_V0_DK_strange);
  free(C_Vi_DK_strange);
  free(C_S0_DK_strange);
  free(C_f0_S_DK_strange);
  
  ///////// FINE FACCIO LA SPLINE
  
  /////////// CHECK GRAFICO

  //////// check per l'evento iev_check_V_S_DK_strange = 2
  int iev_check_V_S_DK_strange = 2;
  
  FILE *fout_graph_check_V0_DK_strange, *fout_graph_check_Vi_DK_strange, *fout_graph_check_S0_DK_strange, *fout_graph_check_f0_S_DK_strange;
  char file_out_graph_check_V0_DK_strange[LEN_NAME], file_out_graph_check_Vi_DK_strange[LEN_NAME], file_out_graph_check_S0_DK_strange[LEN_NAME], file_out_graph_check_f0_S_DK_strange[LEN_NAME];
  
  sprintf(file_out_graph_check_V0_DK_strange, "OUTPUT_SMEAR/%s/%s/check_spline/V0/V0_spline_check.%s.%s.m1m2_sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V_2[ibeta].c_str(), mu_sea_2[imusea].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  sprintf(file_out_graph_check_Vi_DK_strange, "OUTPUT_SMEAR/%s/%s/check_spline/Vi/Vi_spline_check.%s.%s.m1m2_sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V_2[ibeta].c_str(), mu_sea_2[imusea].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  sprintf(file_out_graph_check_S0_DK_strange, "OUTPUT_SMEAR/%s/%s/check_spline/S0/S0_spline_check.%s.%s.m1m2_sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V_2[ibeta].c_str(), mu_sea_2[imusea].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  sprintf(file_out_graph_check_f0_S_DK_strange, "OUTPUT_SMEAR/%s/%s/check_spline/f0_S/f0_S_spline_check.%s.%s.m1m2_sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V_2[ibeta].c_str(), mu_sea_2[imusea].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  
  if ((fout_graph_check_V0_DK_strange = fopen(file_out_graph_check_V0_DK_strange, "w")) == NULL ){
    printf("Error opening the output file file_out_graph_check_V0_DK_strange\n");
    exit(EXIT_FAILURE);
  }
  if ((fout_graph_check_Vi_DK_strange = fopen(file_out_graph_check_Vi_DK_strange, "w")) == NULL ){
    printf("Error opening the output file file_out_graph_check_Vi_DK_strange\n");
    exit(EXIT_FAILURE);
  }
  if ((fout_graph_check_S0_DK_strange = fopen(file_out_graph_check_S0_DK_strange, "w")) == NULL ){
    printf("Error opening the output file file_out_graph_check_S0_DK_strange\n");
    exit(EXIT_FAILURE);
  }
  if ((fout_graph_check_f0_S_DK_strange = fopen(file_out_graph_check_f0_S_DK_strange, "w")) == NULL ){
    printf("Error opening the output file file_out_graph_check_f0_S_DK_strange\n");
    exit(EXIT_FAILURE);
  }
  
  fprintf(fout_graph_check_V0_DK_strange, "@type xydy\n");
  fprintf(fout_graph_check_Vi_DK_strange, "@type xydy\n");
  fprintf(fout_graph_check_S0_DK_strange, "@type xydy\n");
  fprintf(fout_graph_check_f0_S_DK_strange, "@type xydy\n");
  
  for(int ims = i_ms0; ims < i_ms0 + Nmass_strange; ims++){
    
    //  prendendo le masse in GeV
    fprintf(fout_graph_check_V0_DK_strange, "%f %f %f\n", (mq[ibeta][ims]*ainv[ibeta][iev_check_V_S_DK_strange])/Zev[ibeta][iev_check_V_S_DK_strange], V0_DK_strange[ims-i_ms0][iev_check_V_S_DK_strange], sigma_JK_modified_2(V0_DK_strange[ims-i_ms0], 1, 100, clusterfile));
    fprintf(fout_graph_check_Vi_DK_strange, "%f %f %f\n", (mq[ibeta][ims]*ainv[ibeta][iev_check_V_S_DK_strange])/Zev[ibeta][iev_check_V_S_DK_strange], Vi_DK_strange[ims-i_ms0][iev_check_V_S_DK_strange], sigma_JK_modified_2(Vi_DK_strange[ims-i_ms0], 1, 100, clusterfile));
    fprintf(fout_graph_check_S0_DK_strange, "%f %f %f\n", (mq[ibeta][ims]*ainv[ibeta][iev_check_V_S_DK_strange])/Zev[ibeta][iev_check_V_S_DK_strange], S0_DK_strange[ims-i_ms0][iev_check_V_S_DK_strange], sigma_JK_modified_2(S0_DK_strange[ims-i_ms0], 1, 100, clusterfile));
    fprintf(fout_graph_check_f0_S_DK_strange, "%f %f %f\n", (mq[ibeta][ims]*ainv[ibeta][iev_check_V_S_DK_strange])/Zev[ibeta][iev_check_V_S_DK_strange], f0_S_DK_strange[ims-i_ms0][iev_check_V_S_DK_strange], sigma_JK_modified_2(f0_S_DK_strange[ims-i_ms0], 1, 100, clusterfile));
    
    //  prendendo le masse in unità di reticolo
    //fprintf(fout_graph_check_V0_DK_strange, "%f %f %f\n", mq[ibeta][ims], V0_DK_strange[ims-i_ms0][iev_check_V_S_DK_strange], sigma_JK_modified_2(V0_DK_strange[ims-i_ms0], 1, 100, clusterfile));
    //fprintf(fout_graph_check_Vi_DK_strange, "%f %f %f\n", mq[ibeta][ims], Vi_DK_strange[ims-i_ms0][iev_check_V_S_DK_strange], sigma_JK_modified_2(Vi_DK_strange[ims-i_ms0], 1, 100, clusterfile));
    //fprintf(fout_graph_check_S0_DK_strange, "%f %f %f\n", mq[ibeta][ims], S0_DK_strange[ims-i_ms0][iev_check_V_S_DK_strange], sigma_JK_modified_2(S0_DK_strange[ims-i_ms0], 1, 100, clusterfile));
    //fprintf(fout_graph_check_f0_S_DK_strange, "%f %f %f\n", mq[ibeta][ims], f0_S_DK_strange[ims-i_ms0][iev_check_V_S_DK_strange], sigma_JK_modified_2(f0_S_DK_strange[ims-i_ms0], 1, 100, clusterfile));
    
  }// ims
  
  fprintf(fout_graph_check_V0_DK_strange, "&\n");
  fprintf(fout_graph_check_Vi_DK_strange, "&\n");
  fprintf(fout_graph_check_S0_DK_strange, "&\n");
  fprintf(fout_graph_check_f0_S_DK_strange, "&\n");
  
  
  //  prendendo le masse in GeV
  fprintf(fout_graph_check_V0_DK_strange, "%f %f %f\n", mstrange[iev_check_V_S_DK_strange], yspline_V0_DK_strange[iev_check_V_S_DK_strange], sigma_JK_modified_2(yspline_V0_DK_strange, 1, 100, clusterfile));
  fprintf(fout_graph_check_Vi_DK_strange, "%f %f %f\n", mstrange[iev_check_V_S_DK_strange], yspline_Vi_DK_strange[iev_check_V_S_DK_strange], sigma_JK_modified_2(yspline_Vi_DK_strange, 1, 100, clusterfile));
  fprintf(fout_graph_check_S0_DK_strange, "%f %f %f\n", mstrange[iev_check_V_S_DK_strange], yspline_S0_DK_strange[iev_check_V_S_DK_strange], sigma_JK_modified_2(yspline_S0_DK_strange, 1, 100, clusterfile));
  fprintf(fout_graph_check_f0_S_DK_strange, "%f %f %f\n", mstrange[iev_check_V_S_DK_strange], yspline_f0_S_DK_strange[iev_check_V_S_DK_strange], sigma_JK_modified_2(yspline_f0_S_DK_strange, 1, 100, clusterfile));
  
  //  prendendo le masse in unità di reticolo
  //fprintf(fout_graph_check_V0_DK_strange, "%f %f %f\n", (mstrange[iev_check_V_S_DK_strange]/ainv[ibeta][iev_check_V_S_DK_strange])*Zev[ibeta][iev_check_V_S_DK_strange], yspline_V0_DK_strange[iev_check_V_S_DK_strange], sigma_JK_modified_2(yspline_V0_DK_strange, 1, 100, clusterfile));
  //fprintf(fout_graph_check_Vi_DK_strange, "%f %f %f\n", (mstrange[iev_check_V_S_DK_strange]/ainv[ibeta][iev_check_V_S_DK_strange])*Zev[ibeta][iev_check_V_S_DK_strange], yspline_Vi_DK_strange[iev_check_V_S_DK_strange], sigma_JK_modified_2(yspline_Vi_DK_strange, 1, 100, clusterfile));
  //fprintf(fout_graph_check_S0_DK_strange, "%f %f %f\n", (mstrange[iev_check_V_S_DK_strange]/ainv[ibeta][iev_check_V_S_DK_strange])*Zev[ibeta][iev_check_V_S_DK_strange], yspline_S0_DK_strange[iev_check_V_S_DK_strange], sigma_JK_modified_2(yspline_S0_DK_strange, 1, 100, clusterfile));
  //fprintf(fout_graph_check_f0_S_DK_strange, "%f %f %f\n", (mstrange[iev_check_V_S_DK_strange]/ainv[ibeta][iev_check_V_S_DK_strange])*Zev[ibeta][iev_check_V_S_DK_strange], yspline_f0_S_DK_strange[iev_check_V_S_DK_strange], sigma_JK_modified_2(yspline_f0_S_DK_strange, 1, 100, clusterfile));
  
  fclose(fout_graph_check_V0_DK_strange);
  fclose(fout_graph_check_Vi_DK_strange);
  fclose(fout_graph_check_S0_DK_strange);
  fclose(fout_graph_check_f0_S_DK_strange);
  
  /////////// FINE CHECK GRAFICO  

  /////////// OUTPUT

  char file_out_Vi_DK_strange[LEN_NAME], file_out_V0_DK_strange[LEN_NAME], file_out_S0_DK_strange[LEN_NAME], file_out_f0_S_DK_strange[LEN_NAME];
  FILE *fout_Vi_DK_strange, *fout_V0_DK_strange, *fout_S0_DK_strange, *fout_f0_S_DK_strange;

  sprintf(file_out_Vi_DK_strange, "OUTPUT_SMEAR/%s/%s/Vi/Hl_sym_av/%s/%s/Vi.%s.m1m2_sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  sprintf(file_out_V0_DK_strange, "OUTPUT_SMEAR/%s/%s/V0/Hl_sym_av/%s/%s/V0.%s.m1m2_sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  sprintf(file_out_S0_DK_strange, "OUTPUT_SMEAR/%s/%s/S0/Hl_sym_av/%s/%s/S0.%s.m1m2_sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  sprintf(file_out_f0_S_DK_strange, "OUTPUT_SMEAR/%s/%s/f0_S/Hl_sym_av/%s/%s/f0_S.%s.m1m2_sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
  
  if ((fout_Vi_DK_strange = fopen(file_out_Vi_DK_strange, "w")) == NULL ){
    printf("Error opening the output file Vi_DK_strange!!\n");
    exit(EXIT_FAILURE);
  }
  if ((fout_V0_DK_strange = fopen(file_out_V0_DK_strange, "w")) == NULL ){
    printf("Error opening the output file V0_DK_strange!!\n");
    exit(EXIT_FAILURE);
  }
  if ((fout_S0_DK_strange = fopen(file_out_S0_DK_strange, "w")) == NULL ){
    printf("Error opening the output file S0_DK_strange!!\n");
    exit(EXIT_FAILURE);
  }
  if ((fout_f0_S_DK_strange = fopen(file_out_f0_S_DK_strange, "w")) == NULL ){
    printf("Error opening the output file f0_S_DK_strange!!\n");
    exit(EXIT_FAILURE);
  }
  
  for(int iev = 0; iev <= Nev; iev++ ){
    
    fprintf(fout_V0_DK_strange, "%+f\n", yspline_V0_DK_strange[iev]);
    fprintf(fout_S0_DK_strange, "%+f\n", yspline_S0_DK_strange[iev]);
    fprintf(fout_f0_S_DK_strange, "%+f\n", yspline_f0_S_DK_strange[iev]);
    fprintf(fout_Vi_DK_strange, "%+f\n", yspline_Vi_DK_strange[iev]);
    
  }// iev

  free(yspline_V0_DK_strange);
  free(yspline_Vi_DK_strange);
  free(yspline_S0_DK_strange);
  free(yspline_f0_S_DK_strange);
  
  fclose(fout_V0_DK_strange);
  fclose(fout_S0_DK_strange);
  fclose(fout_f0_S_DK_strange);
  fclose(fout_Vi_DK_strange);
  
  /////////// FINE OUTPUT
  
#endif  
  
  ///////// FINE SPLINE V0, Vi, S0, f0_S per il caso D_to_K_strange
  
#endif

  ///////// FINE SPLINE V0, Vi, S0, f0_S per il caso D_to_K

  return 0;
  
}

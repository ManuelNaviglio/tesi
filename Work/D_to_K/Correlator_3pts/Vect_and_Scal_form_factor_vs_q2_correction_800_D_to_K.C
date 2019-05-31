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

//////////////////////////////

int f0_with_S = 1;
int f0_without_S = 0;

//////////////////////////////

//////////////////////////////
   
int energia_da_sinhDR = 0;
int energia_da_stdDR = 1;
int energia_dal_fit = 0;

//////////////////////////////

int ibeta, imusea, no_correction;

int main(){

   ////////  START LETTURA Ultimate Input

  double mlight[Nev+1], mstrange[Nev+1], mcharm[Nev+1], a[Nbeta][Nev+1], ainv[Nbeta][Nev+1], r0[Nev+1], Zev[Nbeta][Nev+1], ZTev[Nbeta][Nev+1], f0[Nev+1], B0[Nev+1], fkfpi[Nev+1];
  int  iboot[Nbeta][Nmusea][Nev+1];

  lettura_ultimate_input( mlight, mstrange, mcharm, a, ainv, r0, Zev, iboot, f0, B0, fkfpi, ZTev);
    

   ////////  FINE LETTURA Ultimate Input


  FILE *fin;

   if ((fin = fopen("Input_corr_3pts/file_input_corr_3pts.out", "r")) == NULL ){
     printf("Error opening the input file!!\n");
     exit(EXIT_FAILURE);
   }
  fscanf(fin, "%d %d %d", &ibeta, &imusea, &no_correction);
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


   /////////// CHIAMO LE MASSE

FILE *fr_M_from_2pts_1, *fr_M_from_2pts_2;
char open_M_from_2pts_fit_1[LEN_NAME], open_M_from_2pts_fit_2[LEN_NAME];

double mass_2pts_m1[Nev+1], mass_2pts_m2[Nev+1];

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

       for(int iev = 0; iev <= Nev; iev++){

	 fscanf(fr_M_from_2pts_1,"%lf\n", &mass_2pts_m1[iev]);  // IN lattice unit
	 fscanf(fr_M_from_2pts_2,"%lf\n", &mass_2pts_m2[iev]);  // IN lattice unit

	 mass_2pts_m1[iev] = mass_2pts_m1[iev]/a[ibeta][iev];  // IN GeV
	 mass_2pts_m2[iev] = mass_2pts_m2[iev]/a[ibeta][iev];  // IN GeV

       }// iev

       fclose(fr_M_from_2pts_1);
       fclose(fr_M_from_2pts_2);


   /////////// FINE CHIAMO LE MASSE



   //////// LEGGO I FATTORI DI FORMA E COSTRUISCO q2

   char open_fzero[LEN_NAME], open_fplus[LEN_NAME];
   FILE *fr_fzero, *fr_fplus;

   double q2[Nth1][Nth2][Nev+1];

   double th1_val, th2_val, momentum_1, momentum_2, q0, qi;
   double energy_2pts_m1, energy_2pts_m2;

   double fzero[Nth1][Nth2][Nev+1], fplus[Nth1][Nth2][Nev+1];
   double sigma_fzero[Nth1][Nth2], sigma_fplus[Nth1][Nth2];



    for(int ith1 = 0; ith1 < Nth1; ith1++){
        for(int ith2 = 0; ith2 < Nth2; ith2++){

       if(no_correction == 0){

	 sprintf(open_fzero, "OUTPUT_SMEAR/%s/%s/fzero_corrected/%s/%s/%s/fzero_corrected.%s.m1m2_sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
	 sprintf(open_fplus, "OUTPUT_SMEAR/%s/%s/fplus_corrected/%s/%s/%s/fplus_corrected.%s.m1m2_sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), th_[ith1].c_str(), th_[ith2].c_str());

       }else if(no_correction == 1){

	 sprintf(open_fzero, "OUTPUT_SMEAR/%s/%s/fzero_corrected/%s/%s/%s/fzero_no_corrected.%s.m1m2_sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), th_[ith1].c_str(), th_[ith2].c_str());
	 sprintf(open_fplus, "OUTPUT_SMEAR/%s/%s/fplus_corrected/%s/%s/%s/fplus_no_corrected.%s.m1m2_sc.th_%s%s.out", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), th_[ith1].c_str(), th_[ith2].c_str());

       }

       if ((fr_fzero = fopen(open_fzero, "r")) == NULL ){
	 printf("Error opening the input file: %s\n",open_fzero);
	 exit(EXIT_FAILURE);
       }
       if(ith1 != 3 || ith2 != 3){
	 if ((fr_fplus = fopen(open_fplus, "r")) == NULL ){
	   printf("Error opening the input file: %s\n",open_fplus);
	   exit(EXIT_FAILURE);
	 }
       }


       th1_val = theta_value[ibeta][ith1];
       th2_val = theta_value[ibeta][ith2];


       for(int iev = 0; iev <= Nev; iev++){

	 fscanf(fr_fzero,"%lf\n", &fzero[ith1][ith2][iev]);

	 if(ith1 != 3 || ith2 != 3){
	   fscanf(fr_fplus,"%lf\n", &fplus[ith1][ith2][iev]);
	 }


	 momentum_1 = (th1_val*PI)/(L[ibeta]*a[ibeta][iev]);                                    // In GeV
	 momentum_2 = (th2_val*PI)/(L[ibeta]*a[ibeta][iev]);                                    // In GeV 

	 energy_2pts_m1 = sqrt( pow(mass_2pts_m1[iev] ,2) + 3*pow(momentum_1 ,2) );             // In GeV
	 energy_2pts_m2 = sqrt( pow(mass_2pts_m2[iev] ,2) + 3*pow(momentum_2 ,2) );             // In GeV

	 q0 = (energy_2pts_m2 - energy_2pts_m1);                                                // In GeV
	 qi = (momentum_2 - momentum_1);                                                        // In GeV

	 q2[ith1][ith2][iev] = pow(q0 ,2) - 3*pow(qi ,2);

       }// iev

     fclose(fr_fzero);
     if(ith1 != 3 || ith2 != 3){
	     fclose(fr_fplus);
     }


       sigma_fzero[ith1][ith2] = sigma_bootstrap( fzero[ith1][ith2], analysis_in, analysis_fin, Nev_an, clusterfile);

       if(ith1 != 3 || ith2 != 3){
	      	sigma_fplus[ith1][ith2] = sigma_bootstrap( fplus[ith1][ith2], analysis_in, analysis_fin, Nev_an, clusterfile);
       }


     } // ith2
   } // ith1

   //////// FINE LEGGO I FATTORI DI FORMA E COSTRUISCO q2



   /////////// COSTRUISCO L'OUTPUT

   char file_out_fzero[LEN_NAME], file_out_fplus[LEN_NAME];
   FILE *fout_fzero, *fout_fplus;

   if(no_correction == 0){

     sprintf(file_out_fzero, "OUTPUT_SMEAR/%s/%s/fzero_corrected/%s/fzero_vs_q2/fzero_corrected_vs_q2.%s.%s.m1m2_sc.xmg", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V_2[ibeta].c_str(), mu_sea_2[imusea].c_str());
     sprintf(file_out_fplus, "OUTPUT_SMEAR/%s/%s/fplus_corrected/%s/fplus_vs_q2/fplus_corrected_vs_q2.%s.%s.m1m2_sc.xmg", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V_2[ibeta].c_str(), mu_sea_2[imusea].c_str());

   } else if(no_correction == 1){

     sprintf(file_out_fzero, "OUTPUT_SMEAR/%s/%s/fzero_corrected/%s/fzero_vs_q2/fzero_no_corrected_vs_q2.%s.%s.m1m2_sc.xmg", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V_2[ibeta].c_str(), mu_sea_2[imusea].c_str());
     sprintf(file_out_fplus, "OUTPUT_SMEAR/%s/%s/fplus_corrected/%s/fplus_vs_q2/fplus_no_corrected_vs_q2.%s.%s.m1m2_sc.xmg", dir_S[iS].c_str(), dir_E[iE].c_str(), strategy[istrategy].c_str(), beta_V_2[ibeta].c_str(), mu_sea_2[imusea].c_str());

   }


   if ((fout_fzero = fopen(file_out_fzero, "w")) == NULL ){
     printf("Error opening the input file: %s\n",file_out_fzero);
     exit(EXIT_FAILURE);
   }
   if ((fout_fplus = fopen(file_out_fplus, "w")) == NULL ){
     printf("Error opening the input file: %s\n",file_out_fplus);
     exit(EXIT_FAILURE);
   }

   fprintf(fout_fzero, "@type xydy\n");
   fprintf(fout_fplus, "@type xydy\n");


    for(int ith1 = 0; ith1 < Nth1; ith1++){
        for(int ith2 = 0; ith2 < Nth2; ith2++){

       fprintf(fout_fzero, "%f %f %f\n", q2[ith1][ith2][0], fzero[ith1][ith2][0], sigma_fzero[ith1][ith2]);

       if(ith1 != 3 || ith2 != 3){

	 fprintf(fout_fplus, "%f %f %f\n", q2[ith1][ith2][0], fplus[ith1][ith2][0], sigma_fplus[ith1][ith2]);

       }// if(ith1 != 3 || ith2 != 3)

     }// ith2

     if(ith1 != 3){

       fprintf(fout_fzero, "&\n");
       fprintf(fout_fplus, "&\n");

     }// if(ith1 != 3)

   }// ith1

   fclose(fout_fzero);
   fclose(fout_fplus);

   /////////// FINE COSTRUISCO L'OUTPUT
    
  return 0;

}



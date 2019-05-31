#include <TMinuit.h>
#include <TF1.h>            // questo serve per la funzione cout
#include <fstream>
#include "lettura_ult_inp.h"
#include "algebra_func.h"
#include "randgen.h"
#include "stat_analysis_func.h"
#include "definitions.h"

#define LEN_NAME 1024
#define PI 3.141592653589793
#define GF 0.000011664

using namespace std;

const int Nrand = 1000;      // UGUALI
const int Nev_sint = 1000;   // UGUALI

int nbin = 1000;

const int Nq = 8;
const int Nev_tot = 3200;
int analisys_in = 1, analisys_fin = 32;
const int Nanalisys = 32;

const int npar = 3;

int block_par = 1;

const int nfix = 2;

double fix_par[nfix] = {10, 20};

void chi2_boot( int &, double *, double &, double *, int);
double chi2_boot_num(double *);

double fplus_func( double, double, double *);
double fplus_func_HFAG( double, double, double *);

double P_phi_fplus_func( double, double, double *);
double P_phi_fplus_func_HFAG( double, double, double *);

//////////////////////////////

int energia_da_sinhDR = 0;
int energia_da_stdDR = 1;
int energia_dal_fit = 0;

//////////////////////////////

//////////////////////////////

int f0_with_S = 1;
int f0_without_S = 0;

//////////////////////////////

int num_points = 0, num_par = 0, dof;

////// QUANTITÀ FISICHE

double MD_phys = 1.867, Mpi_phys = 0.135, MK_phys = 0.4942;

double tp = pow( MD_phys + MK_phys , 2), t0 = ( MD_phys +  MK_phys)*pow( sqrt( MD_phys ) - sqrt( MK_phys ) , 2), z0 = ( sqrt( tp ) - sqrt( tp - t0 ) )/( sqrt( tp ) + sqrt( tp - t0 ) );

double tm = pow( MD_phys - MK_phys , 2);

double MD0_pdg = 1.865, MKm_pdg = 0.494, MDp_pdg = 1.870, MK0_pdg = 0.498;

//////////////////////

////// QUANTITÀ GLOBALI FIT

double f_fit[Nq], q2_fit[Nq], sigma_f_fit[Nq], mc_fit;

///////////////////////

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
  
  ////////  LETTURA Ultimate Input
  
  double mlight[Nev+1], mstrange[Nev+1], mcharm[Nev+1], a[Nbeta][Nev+1], ainv[Nbeta][Nev+1], r0[Nev+1], Zev[Nbeta][Nev+1], ZTev[Nbeta][Nev+1], f0[Nev+1], B0[Nev+1], fkfpi[Nev+1];
  int  iboot[Nbeta][Nmusea][Nev+1];
  
  lettura_ultimate_input( mlight, mstrange, mcharm, a, ainv, r0, Zev, iboot, f0, B0, fkfpi, ZTev);
  
  ////////  FINE LETTURA Ultimate Input
  
  //////////////// LEGGO I DATI SINTETICI      
  
  FILE *fin_sint_f0_q2, *fin_sint_fp_q2;   
  char file_sint_f0_q2[LEN_NAME], file_sint_fp_q2[LEN_NAME];
  
  double fzero_sint[Nq][Nev_tot+1], fplus_sint[Nq][Nev_tot+1], sigma_fzero_sint[Nq], sigma_fplus_sint[Nq];
  
  double q2_max = pow((MD_phys - MK_phys), 2), q2_step, av_fzero=0, av_fplus=0;
  
  printf("FFFF\n");
  
  for(int iq = 0; iq < Nq; iq++){
    
    sprintf(file_sint_f0_q2, "OUTPUT_SMEAR/%s/%s/dati_sintetici/fzero/fzero_sint_q2_%d.out", dir_S[iS].c_str(), dir_E[iE].c_str(), iq);
    sprintf(file_sint_fp_q2, "OUTPUT_SMEAR/%s/%s/dati_sintetici/fplus/fplus_sint_q2_%d.out", dir_S[iS].c_str(), dir_E[iE].c_str(), iq);
    
    if ((fin_sint_f0_q2 = fopen( file_sint_f0_q2, "r")) == NULL ){
      printf("Error opening the output file: file_sint_f0_q2\n");
      exit(EXIT_FAILURE);
    }
    if ((fin_sint_fp_q2 = fopen( file_sint_fp_q2, "r")) == NULL ){
      printf("Error opening the output file: file_sint_fp_q2\n");
      exit(EXIT_FAILURE);
    }
    
    for(int iev = 1; iev <= Nev_tot; iev++){
      
      fscanf(fin_sint_f0_q2,"%lf\n", &fzero_sint[iq][iev]);
      fscanf(fin_sint_fp_q2,"%lf\n", &fplus_sint[iq][iev]);

      av_fzero = av_fzero + fzero_sint[iq][iev]/Nev_tot;
      av_fplus = av_fplus + fplus_sint[iq][iev]/Nev_tot;
      
    }// iev
    
    fzero_sint[iq][0] = av_fzero;
    fplus_sint[iq][0] = av_fplus;
    
    sigma_fzero_sint[iq] = sigma_bootstrap( fzero_sint[iq], analisys_in, analisys_fin, Nev_an, clusterfile);
    sigma_fplus_sint[iq] = sigma_bootstrap( fplus_sint[iq], analisys_in, analisys_fin, Nev_an, clusterfile);
    
    printf("iq2=%d  fzero=%f  sigma_fzero=%f  fplus=%f  sigma_fplus=%f\n", iq, fzero_sint[iq][0], sigma_fzero_sint[iq], fplus_sint[iq][0], sigma_fplus_sint[iq]);
    
    fclose(fin_sint_f0_q2);
    fclose(fin_sint_fp_q2);
    
    av_fzero = 0;
    av_fplus = 0;
    
  } // iq
  
  printf("#fzero\n");  
  printf("@type xydy\n");
  
  for(int iq = 0; iq < Nq; iq++){
    
    q2_step = (q2_max/( (double)  (Nq - 1)))*iq;
    
    printf("%f  %f  %f\n", q2_step, fzero_sint[iq][0], sigma_fzero_sint[iq]);  
  }
  
  printf("#fplus\n");  
  printf("@type xydy\n");

  for(int iq = 0; iq < Nq; iq++){
    
    q2_step = (q2_max/( (double)  (Nq - 1)))*iq;
    
    printf("%f  %f  %f\n", q2_step, fplus_sint[iq][0], sigma_fplus_sint[iq]); 
  }
  //////////////// FINE LEGGO I DATI SINTETICI
  
  /////////// COSTRUISCO LE QUANTITÀ CINEMATICHE
  
  double q2[Nq];
  
  for(int iq = 0; iq < Nq; iq++){
    
    q2[iq] = (q2_max/( (double) (Nq-1) ))*iq;
    
  }

  /////////// FINE COSTRUISCO LE QUANTITÀ CINEMATICHE
  
  ////////////////////////
  //                    //
  //   FIT BOOTSTRAP    //
  //                    //
  ////////////////////////
  
  double epsilon = 0.000001, chi2_f;
  
  double outpar[npar], err[npar];
  
  double step[npar] = {0.01, 0.01, 0.01};
  
  double par[npar] = {0.00, 0.00, 0.00};
  
  double min[npar] = {0.00, 0.00, 0.00};
  
  double max[npar] = {0.00, 0.00, 0.00};
  
  string cpar[npar] = {"a0", "r1", "r2"};    
  
  double par_boot[Nev_tot+1][npar], parameter[npar][Nev_tot+1], chi2_par[npar];
  
  for(int iev = 0; iev <= Nev_tot; iev++){
    
    for(int i = 0; i < Nq; i++){          
      
      q2_fit[i] = q2[i];
      f_fit[i] = fplus_sint[i][iev];      
      
      sigma_f_fit[i] = sigma_fplus_sint[i];
      
    }
    
    mc_fit = mcharm[0];
    
    TMinuit minuit(npar);
    
    minuit.SetFCN(chi2_boot);
    
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
      
      par_boot[iev][ipar] = parameter[ipar][iev]; 
      
      printf("GG%dG iev=%d %s=%f\n", ipar, iev, cpar[ipar].c_str(), parameter[ipar][iev]);
      
    }
    
    for(int ipar = 0; ipar < npar; ipar++){
      
      if( fabs(outpar[ipar]) > epsilon){
	num_par = num_par + 1;
      }
      
      chi2_par[ipar] = parameter[ipar][iev];
    }
    
    chi2_f = chi2_boot_num(chi2_par);
    
    printf("BBBB chi2_ev%d = %f\t dof=%d\n", iev, chi2_f, dof);	   
    
    num_points = 0;
    num_par = 0;
    dof = 0;
    
  }// iev
  
  ////////// AGGIUNTO ALLA FINE
  
  for(int ipar = 0; ipar < npar; ipar++){
    parameter[ipar][0] = 0.0;
  }// ipar
  
  for(int ipar = 0; ipar < npar; ipar++){
    for(int iev = 1; iev <= Nev_tot; iev++){
      parameter[ipar][0] = parameter[ipar][0] + parameter[ipar][iev]/Nev_tot;
    }// iev
  }// ipar
  
  ////////// FINE AGGIUNTO ALLA FINE
  
  ///////////////  SIGMA BOOTSTRAP DEI PARAMETRI
  
  double sigma_parameter[npar];
  
  for(int ipar = 0; ipar < npar; ipar++){
    sigma_parameter[ipar] = sigma_bootstrap( parameter[ipar], analisys_in, analisys_fin, Nev_an, clusterfile);
  }// ipar

  for(int ipar = 0; ipar < npar; ipar++){
    printf("PPPP %s = %+.3e (%+.3e)\n", cpar[ipar].c_str(), parameter[ipar][0], sigma_parameter[ipar]);
  }// ipar
  
  ///////////////  FINE SIGMA BOOTSTRAP DEI PARAMETRI
  
  //////////////// FINE FIT BOOTSTRAP

  
  /////////// CALCOLO LA MATRICE DI COVARIANZA DEI PARAMETRI DI f+
  
  double cov_matrix[npar][npar];
  
  for(int i = 0; i < npar; i++){
    for(int j = i; j < npar; j++){
      
      cov_matrix[i][j] = cov_bootstrap( parameter[i], parameter[j], analisys_in, analisys_fin, Nev_an, clusterfile);   
      
      cov_matrix[j][i] = cov_matrix[i][j];
      
    }// j
  }// i
  
  printf("CVCV\n\n");
  
  for(int i = 0; i < npar; i++){
    for(int j = 0; j < npar; j++){
      printf("%+.5e  ", cov_matrix[i][j]);
    }
    printf("\n");
  }
  
  printf("\n\n");
  
  printf("CRCR\n\n");
  
  for(int i = 0; i < npar; i++){
    for(int j = 0; j < npar; j++){
      printf("%+.5e  ", cov_matrix[i][j]/(sigma_parameter[i]*sigma_parameter[j]));
    }
    printf("\n");
  }
  
  /////////// FINE CALCOLO LA MATRICE DI COVARIANZA DEI PARAMETRI DI f+

  
  //////////  DEFINISCO I PARAMETRI E LA COVARIANZA DI HFAG

  //double Vcs_HFAG = 0.97343;
  //double Vfp0_HFAG = 0.7226;
  double sigma_Vfp0_HFAG = 0.0034;

  double fp0_HFAG = 0.7423, sigma_fp0_HFAG = 0.0035;

  double r1_HFAG = -2.38, r2_HFAG = 4.7, sigma_r1_HFAG = 0.13, sigma_r2_HFAG = 3.0;

  double cov_matrix_HFAG[npar][npar], par_conversion = sigma_fp0_HFAG/sigma_Vfp0_HFAG;

  double par_HFAG[npar] = {fp0_HFAG, r1_HFAG, r2_HFAG};
  double sigma_par_HFAG[npar];
  
  double *cov_m_par_HFAG=(double*)malloc(sizeof(double)*(npar)*(npar));
  
  cov_matrix_HFAG[0][0] = pow(sigma_fp0_HFAG, 2); 
  cov_matrix_HFAG[1][1] = pow(sigma_r1_HFAG, 2);
  cov_matrix_HFAG[2][2] = pow(sigma_r2_HFAG, 2);

  cov_matrix_HFAG[0][1] = -0.0000810835*par_conversion;
  cov_matrix_HFAG[1][0] = cov_matrix_HFAG[0][1];

  cov_matrix_HFAG[0][2] = 0.00512929*par_conversion;
  cov_matrix_HFAG[2][0] = cov_matrix_HFAG[0][2];

  cov_matrix_HFAG[1][2] = -0.310805;
  cov_matrix_HFAG[2][1] = cov_matrix_HFAG[1][2];
  

  for(int i = 0; i < npar; i++){
    sigma_par_HFAG[i] = sqrt(cov_matrix_HFAG[i][i]);
  }

  for(int i = 0; i < npar; i++){
    for(int j = 0; j < npar; j++){
      
      cov_m_par_HFAG[j+npar*i] = cov_matrix_HFAG[i][j];
      
    }// j
  }// i
    
  printf("DET = %.20f\n", determinant( cov_m_par_HFAG, npar) );
  
  //////////  FINE DEFINISCO I PARAMETRI E LA COVARIANZA DI HFAG


  /////////////// GENERO LA DISTRIBUZIONE SISTETICA DEI FATTORI DI FORMA

  double par_sint[npar][Nev_sint], par_sint_2[Nev_sint][npar], av_P_phi_fplus_sint_HFAG = 0, av_P_phi_fplus_sint_lattice = 0, av_fplus_sint_lattice = 0;
  
  const int Nstep = 50;
  
  double P_phi_fplus_plot_HFAG[Nstep][Nev_sint+1], sigma_P_phi_fplus_plot_HFAG[Nstep], P_phi_fplus_plot_lattice[Nstep][Nev_tot+1], sigma_P_phi_fplus_plot_lattice[Nstep];
  
  double P_phi_fplus_plot_HFAG_sup[Nstep], P_phi_fplus_plot_HFAG_inf[Nstep], P_phi_fplus_plot_lattice_sup[Nstep], P_phi_fplus_plot_lattice_inf[Nstep];
  
  double fplus_plot_lattice[Nstep][Nev_tot+1], sigma_fplus_plot_lattice[Nstep];
  
  double fplus_plot_lattice_sup[Nstep], fplus_plot_lattice_inf[Nstep];
  
  double z_plot[Nstep];
  
  num_random_multivariata_3par( par_HFAG, sigma_par_HFAG, cov_m_par_HFAG, npar, par_sint[0], par_sint[1], par_sint[2], Nrand);
  
  free(cov_m_par_HFAG);
  
  /////////////////// CHECK GENERAZIONE GAUSSIANA PARAMETRI
  
  printf("SSSS\n");
  
  for(int ipar = 0; ipar < npar; ipar++){
    
    printf("distr%d={", ipar);
    
    for(int i = 0; i < Nrand; i++){
      
      if( i != Nrand-1 ){
	
	printf("%f,", par_sint[ipar][i]);
	
      }else if( i == Nrand-1 ){
	
	printf("%f", par_sint[ipar][i]);
	
	printf("};\n");
	
      }
      
    }// i
  }// ipar
  
  /////////////////// FINE CHECK GENERAZIONE GAUSSIANA PARAMETRI
  
  
  for(int ipar = 0; ipar < npar; ipar++){
    for(int iev = 0; iev < Nev_sint; iev++){
      
      par_sint_2[iev][ipar] = par_sint[ipar][iev];
      
    }// iev
  }// ipar
  
  
  for(int iq = 0; iq < Nstep; iq++){
    
    q2_step = (q2_max/( (double)  (Nstep - 1)))*iq;    
    
    for(int iev = 1; iev <= Nev_sint; iev++){
      
      P_phi_fplus_plot_HFAG[iq][iev] = P_phi_fplus_func_HFAG( q2_step, mc_fit, par_sint_2[iev-1]);
      
      av_P_phi_fplus_sint_HFAG = av_P_phi_fplus_sint_HFAG + P_phi_fplus_plot_HFAG[iq][iev]/Nev_sint;
      
    }// iev
    
    for(int iev = 1; iev <= Nev_tot; iev++){
      
      P_phi_fplus_plot_lattice[iq][iev] = P_phi_fplus_func( q2_step, mc_fit, par_boot[iev]);
      
      av_P_phi_fplus_sint_lattice = av_P_phi_fplus_sint_lattice + P_phi_fplus_plot_lattice[iq][iev]/Nev_tot;
      
    }// iev
    
    P_phi_fplus_plot_HFAG[iq][0] = av_P_phi_fplus_sint_HFAG;
    
    P_phi_fplus_plot_lattice[iq][0] = av_P_phi_fplus_sint_lattice;
    
    av_P_phi_fplus_sint_HFAG = 0;
    
    av_P_phi_fplus_sint_lattice = 0;    
    
    z_plot[iq] = ( sqrt( tp - q2_step ) - sqrt( tp - t0 ) )/( sqrt( tp - q2_step ) + sqrt( tp - t0 ) );
    
    
    /////////////////// SIGMA DEI FATTORI DI FORMA PER IL PLOT
    
    sigma_P_phi_fplus_plot_HFAG[iq] = sigma_std( P_phi_fplus_plot_HFAG[iq], 1, Nev_sint);
    
    sigma_P_phi_fplus_plot_lattice[iq] = sigma_bootstrap( P_phi_fplus_plot_lattice[iq], analisys_in, analisys_fin, Nev_an, clusterfile);
    
    /////////////////////////// FINE SIGMA DEI FATTORI DI FORMA PER IL PLOT
    
    
    P_phi_fplus_plot_HFAG_sup[iq] = P_phi_fplus_plot_HFAG[iq][0] + sigma_P_phi_fplus_plot_HFAG[iq];
    
    P_phi_fplus_plot_HFAG_inf[iq] = P_phi_fplus_plot_HFAG[iq][0] - sigma_P_phi_fplus_plot_HFAG[iq];
    
    
    P_phi_fplus_plot_lattice_sup[iq] = P_phi_fplus_plot_lattice[iq][0] + sigma_P_phi_fplus_plot_lattice[iq];
    
    P_phi_fplus_plot_lattice_inf[iq] = P_phi_fplus_plot_lattice[iq][0] - sigma_P_phi_fplus_plot_lattice[iq];
    
  }// iq
  
  
  //////// bubblesort
  
  double temp_P_phi_fplus_plot_HFAG_sup, temp_P_phi_fplus_plot_HFAG_inf, temp_P_phi_fplus_plot_lattice_sup, temp_P_phi_fplus_plot_lattice_inf, temp_z_plot;
  
  for(int i = 0; i < Nstep-1; i++){
    for(int j = Nstep-1; j > i; j--){
      
      if( z_plot[j-1] > z_plot[j] ){
	
	temp_z_plot = z_plot[j-1];
	temp_P_phi_fplus_plot_HFAG_sup = P_phi_fplus_plot_HFAG_sup[j-1];
	temp_P_phi_fplus_plot_HFAG_inf = P_phi_fplus_plot_HFAG_inf[j-1];
	temp_P_phi_fplus_plot_lattice_sup = P_phi_fplus_plot_lattice_sup[j-1];
	temp_P_phi_fplus_plot_lattice_inf = P_phi_fplus_plot_lattice_inf[j-1];
	
	z_plot[j-1] = z_plot[j];
	P_phi_fplus_plot_HFAG_sup[j-1] = P_phi_fplus_plot_HFAG_sup[j];
	P_phi_fplus_plot_HFAG_inf[j-1] = P_phi_fplus_plot_HFAG_inf[j];
	P_phi_fplus_plot_lattice_sup[j-1] = P_phi_fplus_plot_lattice_sup[j];
	P_phi_fplus_plot_lattice_inf[j-1] = P_phi_fplus_plot_lattice_inf[j];
	
	z_plot[j] = temp_z_plot;
	P_phi_fplus_plot_HFAG_sup[j] = temp_P_phi_fplus_plot_HFAG_sup;
	P_phi_fplus_plot_HFAG_inf[j] = temp_P_phi_fplus_plot_HFAG_inf;
	P_phi_fplus_plot_lattice_sup[j] = temp_P_phi_fplus_plot_lattice_sup;
	P_phi_fplus_plot_lattice_inf[j] = temp_P_phi_fplus_plot_lattice_inf;
	
      }// if( z_plot[j-1] > z_plot[j] )
      
    }// j
  }// i
  
  //////// fine bubblesort
  
  
  printf("TTTT\n");
  
  printf("#P_phi_fplus_lattice_vs_z\n");
  
  printf("@type xy\n");
  
  for(int iq = 0; iq < Nstep; iq++){
    
    printf("%f %f\n", z_plot[iq], P_phi_fplus_plot_lattice_sup[iq]);
    
  }
  
  for(int iq = Nstep-1; iq >= 0; iq--){
    
    printf("%f %f\n", z_plot[iq], P_phi_fplus_plot_lattice_inf[iq]);
    
  }
  
  printf("#P_phi_fplus_HFAG_vs_z\n");
  
  printf("@type xy\n");
  
  for(int iq = 0; iq < Nstep; iq++){
    
    printf("%f %f\n", z_plot[iq], P_phi_fplus_plot_HFAG_sup[iq]);
    
  }
  
  for(int iq = Nstep-1; iq >= 0; iq--){
    
    printf("%f %f\n", z_plot[iq], P_phi_fplus_plot_HFAG_inf[iq]);
    
  }
  
  /////////////// FINE GENERO LA DISTRIBUZIONE SISTETICA DEI FATTORI DI FORMA
  
  
  ////////// PLOT CONFRONTO DATI SINTETICI
  
  for(int iq = 0; iq < Nstep; iq++){
    
    q2_step = (q2_max/( (double)  (Nstep - 1)))*iq;    
    
    for(int iev = 1; iev <= Nev_tot; iev++){
      
      fplus_plot_lattice[iq][iev] = fplus_func( q2_step, mc_fit, par_boot[iev]);
      
      av_fplus_sint_lattice = av_fplus_sint_lattice + fplus_plot_lattice[iq][iev]/Nev_tot;
      
    }// iev
    
    fplus_plot_lattice[iq][0] = av_fplus_sint_lattice;
    
    av_fplus_sint_lattice = 0;    
    
    /////////////////// SIGMA DEI FATTORI DI FORMA PER IL PLOT
    
    sigma_fplus_plot_lattice[iq] = sigma_bootstrap( fplus_plot_lattice[iq], analisys_in, analisys_fin, Nev_an, clusterfile);
    
    /////////////////////////// FINE SIGMA DEI FATTORI DI FORMA PER IL PLOT
    
    fplus_plot_lattice_sup[iq] = fplus_plot_lattice[iq][0] + sigma_fplus_plot_lattice[iq];
    
    fplus_plot_lattice_inf[iq] = fplus_plot_lattice[iq][0] - sigma_fplus_plot_lattice[iq];
    
  }// iq
  
  
  printf("LLLL\n");
  
  printf("#fplus_lattice_vs_q2\n");
  
  printf("@type xy\n");
  
  for(int iq = 0; iq < Nstep; iq++){
    
    q2_step = (q2_max/( (double)  (Nstep - 1)))*iq;    
    
    printf("%f %f\n", q2_step, fplus_plot_lattice_sup[iq]);
    
  }
  
  for(int iq = Nstep-1; iq >= 0; iq--){
    
    q2_step = (q2_max/( (double)  (Nstep - 1)))*iq;  
    
    printf("%f %f\n", q2_step, fplus_plot_lattice_inf[iq]);
    
  }
  
  printf("#fplus_sint\n");
  printf("@type xydy\n");
  
  for(int iq = 0; iq < Nq; iq++){
    
    q2_step = (q2_max/( (double)  (Nq - 1)))*iq;  
    
    printf("%f %f %f\n", q2_step, fplus_sint[iq][0], sigma_fplus_sint[iq]);

  }
  
  ////////// FINE PLOT CONFRONTO DATI SINTETICI
  
  return 0;
}


void chi2_boot( int &npar, double *deriv, double &f, double *par, int iflag){
  
  f = 0;
  double f2 = 0;  

  for(int i = 0; i < Nq; i++){
    
    f = f + pow( (f_fit[i] - fplus_func( q2_fit[i], mc_fit, par)) , 2)/pow( sigma_f_fit[i], 2);
    
  }// i
  
  //f2 = pow( par[0] - 0.765, 2)/pow( 0.031, 2);

  f = f + f2;

}// chi2_boot


double chi2_boot_num(double *parameters){
  
  double f = 0, f2 = 0;

  for(int i = 0; i < Nq; i++){
    
    f = f + pow( (f_fit[i] - fplus_func( q2_fit[i], mc_fit, parameters)) , 2)/pow( sigma_f_fit[i], 2);

    num_points = num_points + 1;
      
  }// i
      
  //f2 = pow( parameters[0] - 0.765, 2)/pow( 0.031, 2);

  f = f + f2;

  dof = num_points - (num_par - block_par);

  return f/dof;

}// chi2_boot_num


double fplus_func( double q2, double mc, double *par){
  
  double r1, r2, z, fp0;
  double P, P0, phi, phi0, fp;
  double alpha, ph1, ph2, ph3, ph1_0, ph2_0, ph3_0;  
  
  double MDsStar2 = 4.46054;
  
  fp0 = par[0];
  r1 = par[1];
  r2 = par[2];
  
  alpha = sqrt((PI*pow( mc, 2))/3.0);
  
  ///////// phi
  
  ph1 = ( sqrt( tp - q2 ) + sqrt( tp - t0 ) );
  
  ph2 = ( ( tp - q2 )/pow( (tp - t0), 0.25) );
  
  ph3 = pow( ( sqrt( tp - q2 ) + sqrt( tp - tm ) ), 1.5)/pow( sqrt( tp - q2 ) + sqrt( tp ), 5);
  
  phi = alpha*ph1*ph2*ph3;
  
  /////////////
  
  ///////// phi0
  
  ph1_0 = ( sqrt( tp ) + sqrt( tp - t0 ) );
  
  ph2_0 = ( ( tp )/pow( (tp - t0), 0.25) );
  
  ph3_0 = pow( ( sqrt( tp ) + sqrt( tp - tm ) ), 1.5)/pow( sqrt( tp ) + sqrt( tp ), 5);
  
  phi0 = alpha*ph1_0*ph2_0*ph3_0;
  
  ////////////
  
  P = ( sqrt( tp - q2 ) - sqrt( tp - MDsStar2 ) )/( sqrt( tp - q2 ) + sqrt( tp - MDsStar2 ) );
  P0 = ( sqrt( tp ) - sqrt( tp - MDsStar2 ) )/( sqrt( tp ) + sqrt( tp - MDsStar2 ) );
  
  z = ( sqrt( tp - q2 ) - sqrt( tp - t0 ) )/( sqrt( tp - q2 ) + sqrt( tp - t0 ) );
  
  fp = ( 1/(P*phi) )*( ( fp0*P0*phi0  )/( 1.0 + r1*z0 + r2*pow( z0, 2)  ) )*( 1.0 + r1*z + r2*pow( z, 2));
  
  return fp;
  
}// fplus_func


double P_phi_fplus_func( double q2, double mc, double *par){
  
  double r1, r2, z, fp0;
  double P0, phi0, P_phi_fp;
  double alpha, ph1_0, ph2_0, ph3_0;  
  
  double MDsStar2 = 4.46054;
  
  fp0 = par[0];
  r1 = par[1];
  r2 = par[2];
  
  alpha = sqrt((PI*pow( mc, 2))/3.0);
  
  ///////// phi0
  
  ph1_0 = ( sqrt( tp ) + sqrt( tp - t0 ) );
  
  ph2_0 = ( ( tp )/pow( (tp - t0), 0.25) );
  
  ph3_0 = pow( ( sqrt( tp ) + sqrt( tp - tm ) ), 1.5)/pow( sqrt( tp ) + sqrt( tp ), 5);
  
  phi0 = alpha*ph1_0*ph2_0*ph3_0;
  
  ////////////
  
  P0 = ( sqrt( tp ) - sqrt( tp - MDsStar2 ) )/( sqrt( tp ) + sqrt( tp - MDsStar2 ) );
  
  z = ( sqrt( tp - q2 ) - sqrt( tp - t0 ) )/( sqrt( tp - q2 ) + sqrt( tp - t0 ) );
  
  P_phi_fp = ( ( fp0*P0*phi0  )/( 1.0 + r1*z0 + r2*pow( z0, 2)  ) )*( 1.0 + r1*z + r2*pow( z, 2));
  
  return P_phi_fp;
  
}// P_phi_fplus_func


double fplus_func_HFAG( double q2, double mc, double *par){
  
  double r1, r2, z;
  double P, phi, fp, P0, phi0, fp0;
  double alpha, ph1, ph2, ph3, ph1_0, ph2_0, ph3_0;  
  
  double MDsStar2 = 4.46054;
  
  fp0 = par[0];
  r1 = par[1];
  r2 = par[2];
  
  alpha = sqrt((PI*pow( mc, 2))/3.0);
  
  ///////// phi
  
  ph1 = ( sqrt( tp - q2 ) + sqrt( tp - t0 ) );
  
  ph2 = ( ( tp - q2 )/pow( (tp - t0), 0.25) );
  
  ph3 = pow( ( sqrt( tp - q2 ) + sqrt( tp - tm ) ), 1.5)/pow( sqrt( tp - q2 ) + sqrt( tp ), 5);
  
  phi = alpha*ph1*ph2*ph3;
  
  /////////////
  
  ///////// phi0
  
  ph1_0 = ( sqrt( tp ) + sqrt( tp - t0 ) );
  
  ph2_0 = ( ( tp )/pow( (tp - t0), 0.25) );
  
  ph3_0 = pow( ( sqrt( tp ) + sqrt( tp - tm ) ), 1.5)/pow( sqrt( tp ) + sqrt( tp ), 5);
  
  phi0 = alpha*ph1_0*ph2_0*ph3_0;
  
  ////////////
  
  P = ( sqrt( tp - q2 ) - sqrt( tp - MDsStar2 ) )/( sqrt( tp - q2 ) + sqrt( tp - MDsStar2 ) );
  P0 = ( sqrt( tp ) - sqrt( tp - MDsStar2 ) )/( sqrt( tp ) + sqrt( tp - MDsStar2 ) );
  
  z = ( sqrt( tp - q2 ) - sqrt( tp - t0 ) )/( sqrt( tp - q2 ) + sqrt( tp - t0 ) );
  
  fp = ( 1/(P*phi) )*( ( fp0*P0*phi0  )/( 1.0 + r1*z0 + r2*pow( z0, 2)  ) )*( 1.0 + r1*z + r2*pow( z, 2));
  
  return fp;
  
}// fplus_func_HFAG


double P_phi_fplus_func_HFAG( double q2, double mc, double *par){

  double r1, r2, z;
  double P_phi_fp, P0, phi0, fp0;
  double alpha, ph1_0, ph2_0, ph3_0;  

  double MDsStar2 = 4.46054;

  fp0 = par[0];
  r1 = par[1];
  r2 = par[2];

  alpha = sqrt((PI*pow( mc, 2))/3.0);

  ///////// phi0

  ph1_0 = ( sqrt( tp ) + sqrt( tp - t0 ) );

  ph2_0 = ( ( tp )/pow( (tp - t0), 0.25) );

  ph3_0 = pow( ( sqrt( tp ) + sqrt( tp - tm ) ), 1.5)/pow( sqrt( tp ) + sqrt( tp ), 5);

  phi0 = alpha*ph1_0*ph2_0*ph3_0;

  ////////////

  P0 = ( sqrt( tp ) - sqrt( tp - MDsStar2 ) )/( sqrt( tp ) + sqrt( tp - MDsStar2 ) );

  z = ( sqrt( tp - q2 ) - sqrt( tp - t0 ) )/( sqrt( tp - q2 ) + sqrt( tp - t0 ) );

  P_phi_fp = ( ( fp0*P0*phi0  )/( 1.0 + r1*z0 + r2*pow( z0, 2)  ) )*( 1.0 + r1*z + r2*pow( z, 2));

  return P_phi_fp;

}// P_phi_fplus_func_HFAG

#include <TMinuit.h>
#include <TF1.h>            // questo serve per la funzione cout
#include <fstream>
#include "algebra_func.h"
#include "randgen.h"
#include "stat_analysis_func.h"
#include "array_manip.h"
#include "synt_data_func.h"
#include "grace_t.h"
#include "fit_Minuit.h"
#include "definitions.h"

#define LEN_NAME 1024
#define PI 3.141592653589793
#define GF 0.000011664

using namespace std;

const int Nq = 8;
const int Nmatrix = 16;
const int Nrand = 500;      // UGUALI
const int Nev_sint = 500;   // UGUALI
const int Nev_tot = 3200;
int analysis_in = 1, analysis_fin = 32;

string form_factors[2] = {"fzero", "fplus"};

const int npar = 3;
int nbin = 1000;

int output = 0; // IL COMANDO PER ACCENDERE L'OUTPUT

int block_par = 0;
const int nfix = 2;
double fix_par[nfix] = {10, 20};

//////////////////////////////

int energia_da_sinhDR = 0;
int energia_da_stdDR = 1;
int energia_dal_fit = 0;

//////////////////////////////

//////////////////////////////

int f0_with_S = 1;
int f0_without_S = 0;

//////////////////////////////


////// QUANTITÀ FISICHE

double MD_phys = 1.867, Mpi_phys = 0.135, MK_phys = 0.4942;

double tp = pow( MD_phys + MK_phys , 2), t0 = ( MD_phys +  MK_phys)*pow( sqrt( MD_phys ) - sqrt( MK_phys ) , 2), z0 = ( sqrt( tp ) - sqrt( tp - t0 ) )/( sqrt( tp ) + sqrt( tp - t0 ) );

double MD0_pdg = 1.865, MKm_pdg = 0.494, MDp_pdg = 1.870, MK0_pdg = 0.498;

//////////////////////

////////// QUANTITÀ GLOBALI PER IL FIT

double f_fit[Nmatrix-1], q2_fit[Nmatrix-1], inv_cov_fit[Nmatrix-1][Nmatrix-1], sigma_f_fit[Nmatrix-1];

//////////////////////////////////////

void chi2_cov( int &, double *, double &, double *, int);
void chi2_boot( int &, double *, double &, double *, int);
double chi2_cov_num(double *, int, int, int*);
double chi2_boot_num(double *, int, int, int*);

double fzero_func( double, double *);
double fplus_func( double, double *);

double p_cube( double, int);
double bin_integral( double, double, int, double *, int);
double Vcs_func( double, double, double *, double, int);

void read_fpVcs( double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, int, int, int, int);
void read_sqr_dGamma( double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, int, int, int, int, int);
void check_read_data( double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, int, int, int, int, int, int);

void Vckm_from_sqrdGamma_boot( double*, int, int, double*, double*, double [Nev_tot+1][npar], double*, int, double (*)( double, double, double*, double, int));
void Vckm_from_fpVckm_boot( double*, int, int, double*, double [Nev_tot+1][npar], double*, double (*)( double, double*));
void sigma_Vckm_from_sqrdGamma_boot( double*, double*, int, int, int, double*, double*, int, int, int);
double cov_Vckm_exp_func( double*, double*, double*, double*, int, int, int, int, int);
double Vckm_const_fit_with_cut( double*, double*, double*, double*, int, int, int, int, int, int);
double Vckm_const_fit_with_cut_corr( double*, double*, double*, double*, int, int, int, int);
double Vckm_single_fit_with_cut( double*, double*, double*, double*, int, int, int, int, int);
double Vckm_single_fit_with_cut_corr( double*, double*, double*, double*, int, int, int, int);
double sigma_Vckm_single_fit_with_cut( double*, double*, double*, int, int, int, int);
double Xi2_with_cut( double, double*, double*, double*, double*, int, int, int, int, int);
double Xi2_with_cut_corr( double, double*, double*, double*, double*, int, int, int, int);

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


  /////////// COSTRUISCO LE QUANTITÀ CINEMATICHE
  double q2_max = pow( MD_phys - MK_phys, 2), q2_min = 0.0;

  for(int iq = 0; iq < Nmatrix-1; iq++){

		if( iq < Nq ){
    	q2_fit[iq] = (q2_max/( (double)  (Nq-1) ))*iq;
		}else{
      q2_fit[iq] = q2_fit[iq-Nq+1];
		}

  }// iq
  /////////// FINE COSTRUISCO LE QUANTITÀ CINEMATICHE


  //////////////// LEGGO I DATI SINTETICI PER fplus

  double *fplus_sint=(double*)malloc(sizeof(double)*(Nq)*(Nev_tot+1)), sigma_fplus_sint[Nq], *temp_fp=(double*)malloc(sizeof(double)*(Nq-1)*(Nev_tot+1));
  
  read_synt_form_factor( fplus_sint, sigma_fplus_sint, Nq, Nev_tot, Nev_an, analysis_in, analysis_fin, q2_max, dir_S, dir_E, form_factors, clusterfile, iS, iE, 1);
  
		// COSTRUISCO L'ARRAY DI fplus SENZA IL VALORE IN q^2 = 0, CHE CORRISPONDE A f_0(q^2=0) 
	for(int iq = 1; iq < Nq; iq++){
		for(int iev = 1; iev <= Nev_tot; iev++){
			temp_fp[(iq-1)*(Nev_tot+1)+iev] = fplus_sint[iq*(Nev_tot+1)+iev];
		}// iev
	}// iq

  //////////////// FINE LEGGO I DATI SINTETICI PER fplus
  
  
  //////////////// LEGGO I DATI SINTETICI PER fzero
  
  double *fzero_sint=(double*)malloc(sizeof(double)*(Nq)*(Nev_tot+1)), sigma_fzero_sint[Nq];
  
  read_synt_form_factor( fzero_sint, sigma_fzero_sint, Nq, Nev_tot, Nev_an, analysis_in, analysis_fin, q2_max, dir_S, dir_E, form_factors, clusterfile, iS, iE, 0);
  //////////////// FINE LEGGO I DATI SINTETICI PER fzero


  //////////////// CALCOLO LA COVARIANZA E L'INVERSA TRA fzero ED fplus
  
  printf("\n\nCZCP Covariance fzero vs fplus\n");

  int Nq_f0p = Nmatrix-1;
  double *f0p=(double*)malloc(sizeof(double)*(Nq_f0p)*(Nev_tot+1));
  
  double *cov_f0p=(double*)malloc(sizeof(double)*(Nq_f0p)*(Nq_f0p)), *var_f0p=(double*)calloc( 1, sizeof(double)), *inv_cov_f0p=(double*)malloc(sizeof(double)*(Nq_f0p)*(Nq_f0p));
  *var_f0p = 5.0*pow(10, -9);
  
	append_two_arrays( f0p, fzero_sint, temp_fp, Nq, Nq-1, Nev_tot);

  covariance_syntetic_data( cov_f0p, Nq_f0p, Nev_tot, Nev_an, *var_f0p, f0p, analysis_in, analysis_fin, clusterfile, 0);

  LU_invert( inv_cov_f0p, cov_f0p, Nmatrix-1);
  
  free(var_f0p);
	free(temp_fp);
  free(fplus_sint);
  free(fzero_sint);

  //////////////// FINE CALCOLO LA COVARIANZA E L'INVERSA TRA fzero ED fplus


  ////////////////////////
  //                    //
  //   FIT COVARIANZA   //
  //                    //
  ////////////////////////
 

  double *mu=(double*)malloc(sizeof(double)*(npar)), *sigma=(double*)malloc(sizeof(double)*(npar));
  double *m_par=(double*)malloc(sizeof(double)*(npar)*(npar));

  double epsilon = 0.000001;
  double step[npar] = {0.01, 0.01, 0.01};
  double par[npar]  = {0.60, 0.00, 0.00};      ///// CONTROLLA QUESTO
  double min[npar]  = {0.00, 0.00, 0.00};
  double max[npar]  = {0.00, 0.00, 0.00};
  std::string cpar[npar] = {"f(0)", "C0", "C+"};
   
	for(int i = 0; i < Nmatrix-1; i++){
		f_fit[i] = f0p[i*(Nev_tot+1)+0];

    for(int j = 0; j < Nmatrix-1; j++){  
      inv_cov_fit[i][j] = inv_cov_f0p[i*(Nmatrix-1)+j];
    }// j

  }// i

	fit_with_Minuit( mu, sigma, m_par, step, par, min, max, cpar, fix_par, epsilon, npar, nfix, block_par, chi2_cov, chi2_cov_num, "CCCC", 0);

  for(int ipar = 0; ipar < npar; ipar++){
    printf("LLLL  %s = %+.3e (%+.3e)\n", cpar[ipar].c_str(), mu[ipar], sigma[ipar]);
  }// ipar

  printf("\nGGGG Covariance Matrix of parameters from Minuit\n\n");
	print_matrix( npar, npar, m_par, 1);

  free(inv_cov_f0p);
  //////////////// FINE FIT COVARIANZA

  
  ////////////////////////
  //                    //
  //   FIT BOOTSTRAP    //
  //                    //
  ////////////////////////

  double par_boot[Nev_tot+1][npar];
 
  double **parameter=(double**)malloc(sizeof(double)*(npar));
  for(int ipar = 0; ipar < npar; ipar++) parameter[ipar]=(double*)malloc(sizeof(double)*(Nev_tot+1));
 
  double *temp_par=(double*)malloc(sizeof(double)*(npar)), *temp_cov_par=(double*)malloc(sizeof(double)*(npar)*(npar));
  double *temp_sigma_par=(double*)malloc(sizeof(double)*(npar));

	//check_string[0].replace( check_string[0].begin(), check_string[0].end(), "BBBB");

	for(int iev = 0; iev <= Nev_tot; iev++){

    for(int i = 0; i < Nmatrix-1; i++){
			f_fit[i] = f0p[i*(Nev_tot+1)+iev];
			sigma_f_fit[i] = sqrt(cov_f0p[i*(Nmatrix-1)+i]);
    }// i 

		fit_with_Minuit( temp_par, temp_sigma_par, temp_cov_par, step, par, min, max, cpar, fix_par, epsilon, npar, nfix, block_par, chi2_boot, chi2_boot_num, "BBBB", iev);

    for(int ipar = 0; ipar < npar; ipar++){
      parameter[ipar][iev] = temp_par[ipar];
      par_boot[iev][ipar] = temp_par[ipar];
    }// ipar

  }// iev  

	for(int ipar = 0; ipar < npar; ipar++){
    parameter[ipar][0] = 0.0;
    par_boot[0][ipar] = 0.0;

		for(int iev = 1; iev <= Nev_tot; iev++){
      parameter[ipar][0] += parameter[ipar][iev]/Nev_tot;
      par_boot[0][ipar] += par_boot[iev][ipar]/Nev_tot;
    }// iev
  }// ipar  

  free(cov_f0p);
	free(temp_par);
	free(temp_sigma_par);
  free(temp_cov_par);

  /// SIGMA BOOTSTRAP DEI PARAMETRI
  double *sigma_parameter=(double*)malloc(sizeof(double)*(npar));
  
  for(int ipar = 0; ipar < npar; ipar++){
    sigma_parameter[ipar] = sigma_bootstrap( parameter[ipar], analysis_in, analysis_fin, Nev_an, clusterfile);
  }

  printf("\n");
  for(int ipar = 0; ipar < npar; ipar++){
    printf("PPPP  %s = %+.3e (%+.3e)\n", cpar[ipar].c_str(), par_boot[0][ipar], sigma_parameter[ipar]);
  }
  free(sigma_parameter);
  /// FINE SIGMA BOOTSTRAP DEI PARAMETRI

  //////////////// CALCOLO LA COVARIANZA BOOTSTRAP E L'INVERSA DEI PARAMETRI
  
  printf("\nBCBC Covariance Matrix of parameters from Bootstrap\n");

  double *boot_cov_par=(double*)malloc(sizeof(double)*(npar)*(npar)), *var_boot_par=(double*)calloc( 1, sizeof(double)), *temp_boot_par=(double*)malloc(sizeof(double)*(npar)*(Nev_tot+1));
	*var_boot_par = 0.0;

  for(int iev = 0; iev <= Nev_tot; iev++){
    for(int ipar = 0; ipar < npar; ipar++){
      temp_boot_par[ipar*(Nev_tot+1)+iev] = parameter[ipar][iev];
    }// ipar
	}// iev
  
  covariance_syntetic_data( boot_cov_par, npar, Nev_tot, Nev_an, *var_boot_par, temp_boot_par, analysis_in, analysis_fin, clusterfile, 1);

  free(var_boot_par);
  free(temp_boot_par);
  free(parameter);
  //////////////// CALCOLO LA COVARIANZA BOOTSTRAP E L'INVERSA DEI PARAMETRI

  ///// FINE FIT BOOTSTRAP

  /////////////////////////////////////////////////////////////////
  //                                                             //
  //   GENERO LA DISTRIBUZIONE SISTETICA DEI FATTORI DI FORMA    //
  //                                                             //
  /////////////////////////////////////////////////////////////////

  const int Nstep = 50;
  
  double **par_sint=(double**)malloc(sizeof(double)*(npar));
  for(int ipar = 0; ipar < npar; ipar++) par_sint[ipar]=(double*)malloc(sizeof(double)*(Nev_sint));
  
  double **par_sint_2=(double**)malloc(sizeof(double)*(Nev_sint));
  for(int iev = 0; iev < Nev_sint; iev++) par_sint_2[iev]=(double*)malloc(sizeof(double)*(npar));  
  
  double **fzero_plot=(double**)malloc(sizeof(double)*(Nstep));
  for(int iq = 0; iq < Nstep; iq++) fzero_plot[iq]=(double*)malloc(sizeof(double)*(Nev_sint+1));
  
  double **fplus_plot=(double**)malloc(sizeof(double)*(Nstep));
  for(int iq = 0; iq < Nstep; iq++) fplus_plot[iq]=(double*)malloc(sizeof(double)*(Nev_sint+1));

  double *temp_fzero_plot=(double*)malloc(sizeof(double)*(Nstep)), *temp_fplus_plot=(double*)malloc(sizeof(double)*(Nstep));
  double *sigma_fzero_plot=(double*)malloc(sizeof(double)*(Nstep)), *sigma_fplus_plot=(double*)malloc(sizeof(double)*(Nstep));
  
  double q2_step;
  

  //num_random_multivariata_3par( mu, sigma, m_par, npar, par_sint[0], par_sint[1], par_sint[2], Nrand);
  //num_random_multivariata_3par( par_boot[0], sigma, boot_cov_par, npar, par_sint[0], par_sint[1], par_sint[2], Nrand);
  num_random_multivariata_3par( mu, sigma, boot_cov_par, npar, par_sint[0], par_sint[1], par_sint[2], Nrand);

  for(int ipar = 0; ipar < npar; ipar++){
    for(int iev = 0; iev < Nev_sint; iev++){
      par_sint_2[iev][ipar] = par_sint[ipar][iev];
    }// iev
  }//ipar
  
  for(int iq = 0; iq < Nstep; iq++){

    q2_step = (q2_max/( (double)  (Nstep - 1)))*iq;    

    for(int iev = 1; iev <= Nev_sint; iev++){
      fzero_plot[iq][iev] = fzero_func( q2_step, par_sint_2[iev-1]);
      fplus_plot[iq][iev] = fplus_func( q2_step, par_sint_2[iev-1]);

      fzero_plot[iq][0] += fzero_plot[iq][iev]/Nev_sint;
      fplus_plot[iq][0] += fplus_plot[iq][iev]/Nev_sint;
    } // iev

    temp_fzero_plot[iq] = fzero_plot[iq][0];
    temp_fplus_plot[iq] = fplus_plot[iq][0];

    sigma_fzero_plot[iq] = sigma_std( fzero_plot[iq], 1, Nev_sint);
    sigma_fplus_plot[iq] = sigma_std( fplus_plot[iq], 1, Nev_sint);
  }// iq

  printf("TTTT\n");

	print_band_grace( 0.0, q2_max, Nstep, "fplus", temp_fplus_plot, sigma_fplus_plot, 0);
	print_band_grace( 0.0, q2_max, Nstep, "fzero", temp_fzero_plot, sigma_fzero_plot, 0);

  free(mu);
  free(sigma);
  free(m_par);
  free(boot_cov_par);  
  free(par_sint);
  free(par_sint_2);
  free(fzero_plot);
  free(fplus_plot);  
  free(temp_fzero_plot); 
  free(temp_fplus_plot); 
  free(sigma_fzero_plot); 
  free(sigma_fplus_plot); 
  //// FINE GENERO LA DISTRIBUZIONE SISTETICA DEI FATTORI DI FORMA


  //////////////////////
  //                  //
  //   CALCOLO Vcs    //
  //                  //
  //////////////////////

  const int Nset_data = 6;

  const int NBaBar   = 10, Nsubgrp_1 = NBaBar;
  const int NCleo_D0 =  9, Nsubgrp_2 = NBaBar+NCleo_D0;
  const int NCleo_Dp =  9, Nsubgrp_3 = NBaBar+NCleo_D0+NCleo_Dp;
  const int NBes3_D0 = 18, Nsubgrp_4 = NBaBar+NCleo_D0+NCleo_Dp+NBes3_D0;
  const int NBes3_Dp =  9, Nsubgrp_5 = NBaBar+NCleo_D0+NCleo_Dp+NBes3_D0+NBes3_Dp;
  const int NBelle =   27, Nsubgrp_6 = NBaBar+NCleo_D0+NCleo_Dp+NBes3_D0+NBes3_Dp+NBelle;

  const int Ntot = NBaBar+NCleo_D0+NCleo_Dp+NBes3_D0+NBes3_Dp+NBelle, Ncov = NBaBar+NCleo_D0+NCleo_Dp+NBes3_D0+NBes3_Dp;
 
#if defined(BELLE)
	const int Nbin_max = Ntot;
#endif
#if defined(NO_BELLE)
	const int Nbin_max = Ncov;
#endif  
 
  ////////// LEGGO I DATI SPERIMENTALI DI q2, fpVcs, sqr_dGamma E LA MATRICE DI COVARIANZA DELLE sqr_dGamma
  
  double *cov_exp_sqr_dGamma=(double*)malloc(sizeof(double)*(Ncov)*(Ncov));
  std::string file_name_cov_matrix[1] = {"Input_corr_3pts/Matrice_covarianza_sqrt_dGamma_exp_D_to_K.out"};

  ////// sqr_dGamma
  double *sqr_dGamma_BaBar=(double*)malloc(sizeof(double)*(NBaBar)),     *sigma_sqr_dGamma_BaBar=(double*)malloc(sizeof(double)*(NBaBar));
  double *sqr_dGamma_Cleo_D0=(double*)malloc(sizeof(double)*(NCleo_D0)), *sigma_sqr_dGamma_Cleo_D0=(double*)malloc(sizeof(double)*(NCleo_D0));
  double *sqr_dGamma_Cleo_Dp=(double*)malloc(sizeof(double)*(NCleo_Dp)), *sigma_sqr_dGamma_Cleo_Dp=(double*)malloc(sizeof(double)*(NCleo_Dp));
  double *sqr_dGamma_Bes3_D0=(double*)malloc(sizeof(double)*(NBes3_D0)), *sigma_sqr_dGamma_Bes3_D0=(double*)malloc(sizeof(double)*(NBes3_D0));
  double *sqr_dGamma_Bes3_Dp=(double*)malloc(sizeof(double)*(NBes3_Dp)), *sigma_sqr_dGamma_Bes3_Dp=(double*)malloc(sizeof(double)*(NBes3_Dp));

  ////// fpVcs
  double *fpVcs_BaBar=(double*)malloc(sizeof(double)*(NBaBar)),     *sigma_fpVcs_BaBar=(double*)malloc(sizeof(double)*(NBaBar));
  double *fpVcs_Cleo_D0=(double*)malloc(sizeof(double)*(NCleo_D0)), *sigma_fpVcs_Cleo_D0=(double*)malloc(sizeof(double)*(NCleo_D0));
  double *fpVcs_Cleo_Dp=(double*)malloc(sizeof(double)*(NCleo_Dp)), *sigma_fpVcs_Cleo_Dp=(double*)malloc(sizeof(double)*(NCleo_Dp));
  double *fpVcs_Belle=(double*)malloc(sizeof(double)*(NBelle)),     *sigma_fpVcs_Belle=(double*)malloc(sizeof(double)*(NBelle));

  ////// q2
  double *q2_inf_BaBar=(double*)malloc(sizeof(double)*(NBaBar)),     *q2_sup_BaBar=(double*)malloc(sizeof(double)*(NBaBar));
  double *q2_inf_Cleo_D0=(double*)malloc(sizeof(double)*(NCleo_D0)), *q2_sup_Cleo_D0=(double*)malloc(sizeof(double)*(NCleo_D0));
  double *q2_inf_Cleo_Dp=(double*)malloc(sizeof(double)*(NCleo_Dp)), *q2_sup_Cleo_Dp=(double*)malloc(sizeof(double)*(NCleo_Dp));
  double *q2_inf_Bes3_D0=(double*)malloc(sizeof(double)*(NBes3_D0)), *q2_sup_Bes3_D0=(double*)malloc(sizeof(double)*(NBes3_D0));
  double *q2_inf_Bes3_Dp=(double*)malloc(sizeof(double)*(NBes3_Dp)), *q2_sup_Bes3_Dp=(double*)malloc(sizeof(double)*(NBes3_Dp));
  
  double *q2_BaBar=(double*)malloc(sizeof(double)*(NBaBar));
  double *q2_Cleo_D0=(double*)malloc(sizeof(double)*(NCleo_D0));
  double *q2_Cleo_Dp=(double*)malloc(sizeof(double)*(NCleo_Dp));
  double *q2_Bes3_D0=(double*)malloc(sizeof(double)*(NBes3_D0));
  double *q2_Bes3_Dp=(double*)malloc(sizeof(double)*(NBes3_Dp));
  double *q2_Belle=(double*)malloc(sizeof(double)*(NBelle));

  ////// Vcs
  double *Vcs_BaBar=(double*)malloc(sizeof(double)*(NBaBar)*(Nev_tot+1));
  double *Vcs_Cleo_D0=(double*)malloc(sizeof(double)*(NCleo_D0)*(Nev_tot+1));
  double *Vcs_Cleo_Dp=(double*)malloc(sizeof(double)*(NCleo_Dp)*(Nev_tot+1));
  double *Vcs_Bes3_D0=(double*)malloc(sizeof(double)*(NBes3_D0)*(Nev_tot+1));
  double *Vcs_Bes3_Dp=(double*)malloc(sizeof(double)*(NBes3_Dp)*(Nev_tot+1));
  double *Vcs_Belle=(double*)malloc(sizeof(double)*(NBelle)*(Nev_tot+1));

  
	read_fpVcs( q2_BaBar, q2_Cleo_D0, q2_Cleo_Dp, q2_Belle, fpVcs_BaBar, fpVcs_Cleo_D0, fpVcs_Cleo_Dp, fpVcs_Belle, sigma_fpVcs_BaBar, sigma_fpVcs_Cleo_D0, sigma_fpVcs_Cleo_Dp, sigma_fpVcs_Belle, NBaBar, NCleo_D0, NCleo_Dp, NBelle);

	read_sqr_dGamma( q2_inf_BaBar, q2_inf_Cleo_D0, q2_inf_Cleo_Dp, q2_inf_Bes3_D0, q2_inf_Bes3_Dp, q2_sup_BaBar, q2_sup_Cleo_D0, q2_sup_Cleo_Dp, q2_sup_Bes3_D0, q2_sup_Bes3_Dp, q2_BaBar, q2_Cleo_D0, q2_Cleo_Dp, q2_Bes3_D0, q2_Bes3_Dp, sqr_dGamma_BaBar, sqr_dGamma_Cleo_D0, sqr_dGamma_Cleo_Dp, sqr_dGamma_Bes3_D0, sqr_dGamma_Bes3_Dp, sigma_sqr_dGamma_BaBar, sigma_sqr_dGamma_Cleo_D0, sigma_sqr_dGamma_Cleo_Dp, sigma_sqr_dGamma_Bes3_D0, sigma_sqr_dGamma_Bes3_Dp, NBaBar, NCleo_D0, NCleo_Dp, NBes3_D0, NBes3_Dp);

	read_matrix( Ncov, Ncov, cov_exp_sqr_dGamma, file_name_cov_matrix);

  //printf("\nH555\n\n");
  //check_read_data( q2_inf_BaBar, q2_inf_Cleo_D0, q2_inf_Cleo_Dp, q2_inf_Bes3_D0, q2_inf_Bes3_Dp, q2_sup_BaBar, q2_sup_Cleo_D0, q2_sup_Cleo_Dp, q2_sup_Bes3_D0, q2_sup_Bes3_Dp, q2_BaBar, q2_Cleo_D0, q2_Cleo_Dp, q2_Bes3_D0, q2_Bes3_Dp, q2_Belle, sqr_dGamma_BaBar, sqr_dGamma_Cleo_D0, sqr_dGamma_Cleo_Dp, sqr_dGamma_Bes3_D0, sqr_dGamma_Bes3_Dp, sigma_sqr_dGamma_BaBar, sigma_sqr_dGamma_Cleo_D0, sigma_sqr_dGamma_Cleo_Dp, sigma_sqr_dGamma_Bes3_D0, sigma_sqr_dGamma_Bes3_Dp, fpVcs_BaBar, fpVcs_Cleo_D0, fpVcs_Cleo_Dp, fpVcs_Belle, sigma_fpVcs_BaBar, sigma_fpVcs_Cleo_D0, sigma_fpVcs_Cleo_Dp, sigma_fpVcs_Belle, NBaBar, NCleo_D0, NCleo_Dp, NBes3_D0, NBes3_Dp, NBelle);

  //printf("\nCVCV\n\n");
  //print_matrix( Ncov, Ncov, cov_exp_sqr_dGamma, 1);

  free(fpVcs_BaBar);
  free(fpVcs_Cleo_D0);
  free(fpVcs_Cleo_Dp);
  free(sigma_fpVcs_BaBar);
  free(sigma_fpVcs_Cleo_D0);
  free(sigma_fpVcs_Cleo_Dp);
  ////////// FINE LEGGO I DATI SPERIMENTALI DI q2, fpVcs, sqr_dGamma E LA MATRICE DI COVARIANZA DELLE sqr_dGamma

  //////// CALCOLO Vcs AI VALORI DI q2 RELATIVI AGLI ESPERIMENTI

  // BABAR
	Vckm_from_sqrdGamma_boot( Vcs_BaBar, Nev_tot, NBaBar, q2_inf_BaBar, q2_sup_BaBar, par_boot, sqr_dGamma_BaBar, 0, Vcs_func);

  // CLEO-D0
	Vckm_from_sqrdGamma_boot( Vcs_Cleo_D0, Nev_tot, NCleo_D0, q2_inf_Cleo_D0, q2_sup_Cleo_D0, par_boot, sqr_dGamma_Cleo_D0, 0, Vcs_func);

  // CLEO-D+
	Vckm_from_sqrdGamma_boot( Vcs_Cleo_Dp, Nev_tot, NCleo_Dp, q2_inf_Cleo_Dp, q2_sup_Cleo_Dp, par_boot, sqr_dGamma_Cleo_Dp, 1, Vcs_func);

  // BES3-D0
	Vckm_from_sqrdGamma_boot( Vcs_Bes3_D0, Nev_tot, NBes3_D0, q2_inf_Bes3_D0, q2_sup_Bes3_D0, par_boot, sqr_dGamma_Bes3_D0, 0, Vcs_func);

  // BES3-D+
	Vckm_from_sqrdGamma_boot( Vcs_Bes3_Dp, Nev_tot, NBes3_Dp, q2_inf_Bes3_Dp, q2_sup_Bes3_Dp, par_boot, sqr_dGamma_Bes3_Dp, 1, Vcs_func);

  // BELLE
	Vckm_from_fpVckm_boot( Vcs_Belle, Nev_tot, NBelle, q2_Belle, par_boot, fpVcs_Belle, fplus_func);

  free(q2_inf_BaBar);
  free(q2_sup_BaBar);
  free(q2_inf_Cleo_D0);
  free(q2_sup_Cleo_D0);
  free(q2_inf_Cleo_Dp);
  free(q2_sup_Cleo_Dp);
  free(q2_inf_Bes3_D0);
  free(q2_sup_Bes3_D0);
  free(q2_inf_Bes3_Dp);
  free(q2_sup_Bes3_Dp);
  //////// FINE CALCOLO Vcs AI VALORI DI q2 RELATIVI AGLI ESPERIMENTI

  //////////// CREO IL VETTORE DI TUTTI I TERMINI

  double *q2_all=(double*)malloc(sizeof(double)*(Ntot)), *q2_all_sort=(double*)malloc(sizeof(double)*(Ntot));
  double *sqr_dGamma_all=(double*)malloc(sizeof(double)*(Ntot)), *sigma_sqr_dGamma_all=(double*)malloc(sizeof(double)*(Ntot));
  double *Vcs_all=(double*)malloc(sizeof(double)*(Ntot)*(Nev_tot+1)), *Vcs_all_sort=(double*)malloc(sizeof(double)*(Ntot));
  
  //q2_all
	append_six_arrays( q2_all, q2_BaBar, q2_Cleo_D0, q2_Cleo_Dp, q2_Bes3_D0, q2_Bes3_Dp, q2_Belle, NBaBar, NCleo_D0, NCleo_Dp, NBes3_D0, NBes3_Dp, NBelle, 0);

  //sqr_dGamma_all
	append_six_arrays( sqr_dGamma_all, sqr_dGamma_BaBar, sqr_dGamma_Cleo_D0, sqr_dGamma_Cleo_Dp, sqr_dGamma_Bes3_D0, sqr_dGamma_Bes3_Dp, fpVcs_Belle, NBaBar, NCleo_D0, NCleo_Dp, NBes3_D0, NBes3_Dp, NBelle, 0);

  //sigma_sqr_dGamma_all
	append_six_arrays( sigma_sqr_dGamma_all, sigma_sqr_dGamma_BaBar, sigma_sqr_dGamma_Cleo_D0, sigma_sqr_dGamma_Cleo_Dp, sigma_sqr_dGamma_Bes3_D0, sigma_sqr_dGamma_Bes3_Dp, sigma_fpVcs_Belle, NBaBar, NCleo_D0, NCleo_Dp, NBes3_D0, NBes3_Dp, NBelle, 0);

  //Vcs_all
	append_six_arrays( Vcs_all, Vcs_BaBar, Vcs_Cleo_D0, Vcs_Cleo_Dp, Vcs_Bes3_D0, Vcs_Bes3_Dp, Vcs_Belle, NBaBar, NCleo_D0, NCleo_Dp, NBes3_D0, NBes3_Dp, NBelle, Nev_tot);

  // ORDINO I VALORI IN q2
  for(int i = 0; i < Ntot; i++){
    q2_all_sort[i] = q2_all[i];
  }// i
#if defined(VCS_COV_LAT_TOT)
	sort_vector( q2_all_sort, Ntot, q2_all);
#endif

  free(Vcs_BaBar);
  free(Vcs_Cleo_D0);
  free(Vcs_Cleo_Dp);
  free(Vcs_Bes3_D0);
  free(Vcs_Bes3_Dp);
  free(Vcs_Belle);
  free(sqr_dGamma_BaBar);
  free(sqr_dGamma_Cleo_D0);
  free(sqr_dGamma_Cleo_Dp);
  free(sqr_dGamma_Bes3_D0);
  free(sqr_dGamma_Bes3_Dp);
  free(sigma_sqr_dGamma_BaBar);
  free(sigma_sqr_dGamma_Cleo_D0);
  free(sigma_sqr_dGamma_Cleo_Dp);
  free(sigma_sqr_dGamma_Bes3_D0);
  free(sigma_sqr_dGamma_Bes3_Dp);
  free(fpVcs_Belle);
  free(sigma_fpVcs_Belle);
  free(q2_BaBar);
  free(q2_Cleo_D0);
  free(q2_Cleo_Dp);
  free(q2_Bes3_D0);
  free(q2_Bes3_Dp);
  free(q2_Belle);
  //////////// FINE CREO IL VETTORE DI TUTTI I TERMINI


  //// MATRICE DI COVARIANZA DOVUTA AI CONTRIBUTI DI RETICOLO PER Vcs
  
  int Size_blks[Nset_data] = { NBaBar, NCleo_D0, NCleo_Dp, NBes3_D0, NBes3_Dp, NBelle};
  
  double *Vcs_cov_lat=(double*)malloc(sizeof(double)*(Ntot)*(Ntot)), *Vcs_cov_lat_block=(double*)malloc(sizeof(double)*(Ntot)*(Ntot)), *var_Vckm=(double*)calloc( 1, sizeof(double));
  *var_Vckm = 0.0;
  
  covariance_syntetic_data( Vcs_cov_lat, Ntot, Nev_tot, Nev_an, *var_Vckm, Vcs_all, analysis_in, analysis_fin, clusterfile, 0);
  free(var_Vckm);

  //// COSTRUISCO LA MATRICE Vcs_cov_lat_block
  for(int i = 0; i < Ntot; i++){
    for(int j = 0; j < Ntot; j++){
      Vcs_cov_lat_block[i*Ntot+j] = Vcs_cov_lat[i*Ntot+j];
    }// j
  }// i

	block_matrix( Vcs_cov_lat_block, Ntot, Nset_data, Size_blks);
  /////////////////////////////////////////////////////////////


#if defined(NO_VCS_COV_LAT)
  ///// PRENDO SOLO LA DIAGONALE
  for(int i = 0; i < Ntot; i++){
    for(int j = 0; j < Ntot; j++){
      Vcs_cov_lat[i*Ntot+j] = 0.0;
    }// j
  }// i
#endif
  
#if defined(VCS_COV_LAT_DIAG)
  ///// PRENDO SOLO LA DIAGONALE
  for(int i = 0; i < Ntot; i++){
    for(int j = 0; j < Ntot; j++){
      if( i != j ){
				Vcs_cov_lat[i*Ntot+j] = 0.0;
      }            
    }// j
  }// i
#endif

#if defined(VCS_COV_LAT_BLOCK)
	block_matrix( Vcs_cov_lat, Ntot, Nset_data, Size_blks);
#endif  

  //printf("\n\nLATC Vcs_cov_lat \n\n");
  //print_matrix_for_mathematica( Ntot, Ntot, Vcs_cov_lat, 0);

  //// FINE MATRICE DI COVARIANZA DOVUTA AI CONTRIBUTI DI RETICOLO PER Vcs


  //////////////////////  MEDIA PESATA PER RICAVARE Vcs
  double *Vcs_cov=(double*)malloc(sizeof(double)*(Ntot)*(Ntot)), *Vcs_cov_sort=(double*)malloc(sizeof(double)*(Ntot)*(Ntot));
  double *Vcs_cov_exp=(double*)malloc(sizeof(double)*(Ntot)*(Ntot)); // La matrice viene calcolata evento per evento quando faccio il fit

  double *Vcs_fit_BaBar=(double*)malloc(sizeof(double)*(Nev_tot+1)),      *sigma_Vcs_fit_BaBar=(double*)malloc(sizeof(double));
  double *Vcs_fit_Cleo_D0=(double*)malloc(sizeof(double)*(Nev_tot+1)),    *sigma_Vcs_fit_Cleo_D0=(double*)malloc(sizeof(double));
  double *Vcs_fit_Cleo_Dp=(double*)malloc(sizeof(double)*(Nev_tot+1)),    *sigma_Vcs_fit_Cleo_Dp=(double*)malloc(sizeof(double));
  double *Vcs_fit_Bes3_D0=(double*)malloc(sizeof(double)*(Nev_tot+1)),    *sigma_Vcs_fit_Bes3_D0=(double*)malloc(sizeof(double));
  double *Vcs_fit_Bes3_Dp=(double*)malloc(sizeof(double)*(Nev_tot+1)),    *sigma_Vcs_fit_Bes3_Dp=(double*)malloc(sizeof(double));
  double *Vcs_fit_Belle=(double*)malloc(sizeof(double)*(Nev_tot+1)),      *sigma_Vcs_fit_Belle=(double*)malloc(sizeof(double));
  double *Vcs_fit_correlated=(double*)malloc(sizeof(double)*(Nev_tot+1)), *sigma_Vcs_fit_correlated=(double*)malloc(sizeof(double));

  const int Ncut = 20;
  double q2_cut_min = 0.1, q2_cut_max = q2_max, *q2_cut=(double*)malloc(sizeof(double)*(Ncut));
  double *Vcs_cut=(double*)malloc(sizeof(double)*(Ncut)), *sigma_Vcs_cut=(double*)malloc(sizeof(double)*(Ncut));

  for(int ic = 0; ic < Ncut; ic++){
    q2_cut[ic] = q2_cut_min + ((q2_cut_max - q2_cut_min)/( (double) (Ncut - 1)))*ic;
  }// ic
  
  for(int ic = 0; ic < Ncut; ic++){
  
  	Vcs_fit_BaBar[0]      = 0.0; 
  	Vcs_fit_Cleo_D0[0]    = 0.0;
  	Vcs_fit_Cleo_Dp[0]    = 0.0;
  	Vcs_fit_Bes3_D0[0]    = 0.0;
  	Vcs_fit_Bes3_Dp[0]    = 0.0;
  	Vcs_fit_Belle[0]      = 0.0;
  	Vcs_fit_correlated[0] = 0.0;
  
  	for(int iev = 1; iev <= Nev_tot; iev++){

  		//// MATRICE DI COVARIANZA DI Vcs DOVUTA AI CONTRIBUTI DELLE DELTA GAMMA SPERIMENTALI
  		for(int i = 0; i < Ntot; i++){
    		for(int j = 0; j < Ntot; j++){
					Vcs_cov_exp[i*Ntot+j] = cov_Vckm_exp_func( Vcs_all, sqr_dGamma_all, sigma_sqr_dGamma_all, cov_exp_sqr_dGamma, i, j, Nev_tot, iev, Ncov);
    		}// j
  		}// i  

    	/////// COSTRUISCO LA MATRICE DI COVARIANZA SORT
    	for(int i = 0; i < Ntot; i++){

      	Vcs_all_sort[i] = Vcs_all[i*(Nev_tot+1)+iev];
      
      	for(int j = 0; j < Ntot; j++){

					Vcs_cov[i*Ntot+j] = Vcs_cov_lat_block[i*Ntot+j] + Vcs_cov_exp[i*Ntot+j];
	
					Vcs_cov_sort[i*Ntot+j] = Vcs_cov_lat[i*Ntot+j] + Vcs_cov_exp[i*Ntot+j];

      	}// j
    	}// i


#if defined(VCS_COV_LAT_TOT)
			////// SORT DI Vcs
			sort_vector( Vcs_all_sort, Ntot, q2_all);

			//// SORT DELLA MATRICE DI COVARIANZA
			sort_matrix( Vcs_cov_sort, Ntot, Ntot, q2_all, q2_all);

			//printf("\n\nWWWW Vcs_cov_sort \n\n");
			//print_matrix( Ntot, Ntot, Vcs_cov_sort, 1);
#endif  

			for(int i = 0; i < Ntot; i++){
				for(int j = 0; j < Ntot; j++){
	
					if( q2_all[i] > q2_cut[ic] || q2_all[j] > q2_cut[ic]){
	  
						if(i == j){
							Vcs_cov[i*Ntot+j] = 1.0;
						}else if(i != j){
							Vcs_cov[i*Ntot+j] = 0.0;
						}
	  
					}// if( q2_all[i] > q2_cut[ic] && q2_all[j] > q2_cut[ic])
	
				}// j
			}// i
    
			for(int i = 0; i < Ntot; i++){
				for(int j = 0; j < Ntot; j++){
	
					if( q2_all_sort[i] > q2_cut[ic] || q2_all_sort[j] > q2_cut[ic]){
	  
						if(i == j){
							Vcs_cov_sort[i*Ntot+j] = 1.0;
						}else if(i != j){
							Vcs_cov_sort[i*Ntot+j] = 0.0;
						}
	  
					}// if
	
				}// j
			}// i  

			Vcs_fit_BaBar[iev]      = Vckm_const_fit_with_cut( Vcs_all, Vcs_cov, q2_all, q2_cut, Ntot, Nev_tot,       0.0, Nsubgrp_1, iev, ic);
			Vcs_fit_Cleo_D0[iev]    = Vckm_const_fit_with_cut( Vcs_all, Vcs_cov, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_1, Nsubgrp_2, iev, ic);
			Vcs_fit_Cleo_Dp[iev]    = Vckm_const_fit_with_cut( Vcs_all, Vcs_cov, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_2, Nsubgrp_3, iev, ic);
			Vcs_fit_Bes3_D0[iev]    = Vckm_const_fit_with_cut( Vcs_all, Vcs_cov, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_3, Nsubgrp_4, iev, ic);
			Vcs_fit_Bes3_Dp[iev]    = Vckm_const_fit_with_cut( Vcs_all, Vcs_cov, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_4, Nsubgrp_5, iev, ic);
			Vcs_fit_Belle[iev]      = Vckm_const_fit_with_cut( Vcs_all, Vcs_cov, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_5, Nsubgrp_6, iev, ic);

			Vcs_fit_correlated[iev] = Vckm_const_fit_with_cut_corr( Vcs_all_sort, Vcs_cov_sort, q2_all_sort, q2_cut, Ntot, 0, Nbin_max, ic);

    	Vcs_fit_BaBar[0]      +=  Vcs_fit_BaBar[iev]/Nev_tot; 
    	Vcs_fit_Cleo_D0[0]    +=  Vcs_fit_Cleo_D0[iev]/Nev_tot; 
    	Vcs_fit_Cleo_Dp[0]    +=  Vcs_fit_Cleo_Dp[iev]/Nev_tot; 
    	Vcs_fit_Bes3_D0[0]    +=  Vcs_fit_Bes3_D0[iev]/Nev_tot; 
    	Vcs_fit_Bes3_Dp[0]    +=  Vcs_fit_Bes3_Dp[iev]/Nev_tot; 
    	Vcs_fit_Belle[0]      +=  Vcs_fit_Belle[iev]/Nev_tot; 
    	Vcs_fit_correlated[0] +=  Vcs_fit_correlated[iev]/Nev_tot; 
 
  	}// iev

    if(ic == Ncut-1){
			*sigma_Vcs_fit_BaBar       =  sigma_bootstrap( Vcs_fit_BaBar,      analysis_in, analysis_fin, Nev_an, clusterfile);
			*sigma_Vcs_fit_Cleo_D0     =  sigma_bootstrap( Vcs_fit_Cleo_D0,    analysis_in, analysis_fin, Nev_an, clusterfile);
			*sigma_Vcs_fit_Cleo_Dp     =  sigma_bootstrap( Vcs_fit_Cleo_Dp,    analysis_in, analysis_fin, Nev_an, clusterfile);
			*sigma_Vcs_fit_Bes3_D0     =  sigma_bootstrap( Vcs_fit_Bes3_D0,    analysis_in, analysis_fin, Nev_an, clusterfile);
			*sigma_Vcs_fit_Bes3_Dp     =  sigma_bootstrap( Vcs_fit_Bes3_Dp,    analysis_in, analysis_fin, Nev_an, clusterfile);
			*sigma_Vcs_fit_Belle       =  sigma_bootstrap( Vcs_fit_Belle,      analysis_in, analysis_fin, Nev_an, clusterfile);
    }// if
	  *sigma_Vcs_fit_correlated    =  sigma_bootstrap( Vcs_fit_correlated, analysis_in, analysis_fin, Nev_an, clusterfile);

		Vcs_cut[ic]     =  Vcs_fit_correlated[0];
		sigma_Vcs_cut[ic] =  *sigma_Vcs_fit_correlated;

  }// ic
    
  printf("RERE\n");
  printf("Vcs_fit_BaBar=%f sigma_Vcs_fit_BaBar=%f\n",           Vcs_fit_BaBar[0],      *sigma_Vcs_fit_BaBar);
  printf("Vcs_fit_Cleo_D0=%f sigma_Vcs_fit_Cleo_D0=%f\n",       Vcs_fit_Cleo_D0[0],    *sigma_Vcs_fit_Cleo_D0);
  printf("Vcs_fit_Cleo_Dp=%f sigma_Vcs_fit_Cleo_Dp=%f\n",       Vcs_fit_Cleo_Dp[0],    *sigma_Vcs_fit_Cleo_Dp);
  printf("Vcs_fit_Bes3_D0=%f sigma_Vcs_fit_Bes3_D0=%f\n",       Vcs_fit_Bes3_D0[0],    *sigma_Vcs_fit_Bes3_D0);
  printf("Vcs_fit_Bes3_Dp=%f sigma_Vcs_fit_Bes3_Dp=%f\n",       Vcs_fit_Bes3_Dp[0],    *sigma_Vcs_fit_Bes3_Dp);
  printf("Vcs_fit_Belle=%f sigma_Vcs_fit_Belle=%f\n",           Vcs_fit_Belle[0],      *sigma_Vcs_fit_Belle);
  printf("Vcs_fit_correlated=%f sigma_Vcs_fit_correlated=%f\n", Vcs_fit_correlated[0], *sigma_Vcs_fit_correlated);
 
  free(Vcs_fit_BaBar);
  free(Vcs_fit_Cleo_D0);
  free(Vcs_fit_Cleo_Dp);
  free(Vcs_fit_Bes3_D0);
  free(Vcs_fit_Bes3_Dp);
  free(Vcs_fit_Belle);
  free(sigma_Vcs_fit_BaBar);
  free(sigma_Vcs_fit_Cleo_D0);
  free(sigma_Vcs_fit_Cleo_Dp);
  free(sigma_Vcs_fit_Bes3_D0);
  free(sigma_Vcs_fit_Bes3_Dp);
  free(sigma_Vcs_fit_Belle);
  //////////////////////  FINE MEDIA PESATA PER RICAVARE Vcs
  //////////////////////  FINE CALCOLO Vcs


  //////////// SCRIVO IL FILE BOOTSTRAP PER Vcs
  char *file_out_Vcs=(char*)malloc(sizeof(char)*(LEN_NAME));
  FILE *fout_Vcs;

  if(output == 1){

    sprintf(file_out_Vcs, "OUTPUT_SMEAR/%s/%s/CKM/Vcs.out", dir_S[iS].c_str(), dir_E[iE].c_str());

    if ((fout_Vcs = fopen(file_out_Vcs, "w")) == NULL ){
      printf("Error opening the input file: file_out_Vcs\n");
      exit(EXIT_FAILURE);
    }
    
    for(int iev = 1; iev <= Nev_tot; iev++){
      fprintf(fout_Vcs, "%f\n", Vcs_fit_correlated[iev]);
    }
    fclose(fout_Vcs);
    free(file_out_Vcs);
  }// if(output == 1)
  //////////// FINE SCRIVO IL FILE BOOTSTRAP PER Vcs

 
  /////////// PLOT GRACE MEDIA PESATA Vcs
  double *temp_Vcs_all=(double*)malloc(sizeof(double)*(Ntot)), *temp_sigma_Vcs_all=(double*)malloc(sizeof(double)*(Ntot));

  for(int i = 0; i < Ntot; i++){
		temp_Vcs_all[i] = Vcs_all[i*(Nev_tot+1)+0];
		temp_sigma_Vcs_all[i] = sqrt(Vcs_cov[i+Ntot*i]);
	}// i

  printf("\n#HHHH\n");

	print_band_grace( q2_cut_min, q2_cut_max, Ncut, "#CUT_STEP", Vcs_cut, sigma_Vcs_cut, 0);

	print_points_grace( q2_cut, Vcs_cut, 0, Ncut-1, "Vcs_cut", 0);

	print_points_with_err_grace( q2_all, temp_Vcs_all, temp_sigma_Vcs_all,       0.0, Nsubgrp_1-1,   "Vcs_BaBar", 0);
	print_points_with_err_grace( q2_all, temp_Vcs_all, temp_sigma_Vcs_all, Nsubgrp_1, Nsubgrp_2-1, "Vcs_Cleo_D0", 1);
	print_points_with_err_grace( q2_all, temp_Vcs_all, temp_sigma_Vcs_all, Nsubgrp_2, Nsubgrp_3-1, "Vcs_Cleo_Dp", 1);
	print_points_with_err_grace( q2_all, temp_Vcs_all, temp_sigma_Vcs_all, Nsubgrp_3, Nsubgrp_4-1, "Vcs_Bes3_D0", 1);
	print_points_with_err_grace( q2_all, temp_Vcs_all, temp_sigma_Vcs_all, Nsubgrp_4, Nsubgrp_5-1, "Vcs_Bes3_Dp", 1);
	print_points_with_err_grace( q2_all, temp_Vcs_all, temp_sigma_Vcs_all, Nsubgrp_5, Nsubgrp_6-1,   "Vcs_Belle", 1);

	print_horizontal_mean_with_err_grace( q2_min, q2_max, Vcs_fit_correlated[0], *sigma_Vcs_fit_correlated, "Vcs Final Result", 0);

  free(temp_Vcs_all);
  free(temp_sigma_Vcs_all);
  free(Vcs_cut);
  free(sigma_Vcs_cut);
  free(Vcs_fit_correlated);
  free(sigma_Vcs_fit_correlated);
  /////////// PLOT GRACE MEDIA PESATA Vcs


  /////////// ESTRAZIONE Vckm DAL FIT SINGOLO
  double Vcs_single_BaBar,   sigma_Vcs_single_BaBar,   Xi2_BaBar,   Xi2_part_BaBar   = 0.0;
  double Vcs_single_Cleo_D0, sigma_Vcs_single_Cleo_D0, Xi2_Cleo_D0, Xi2_part_Cleo_D0 = 0.0;
  double Vcs_single_Cleo_Dp, sigma_Vcs_single_Cleo_Dp, Xi2_Cleo_Dp, Xi2_part_Cleo_Dp = 0.0;
  double Vcs_single_Bes3_D0, sigma_Vcs_single_Bes3_D0, Xi2_Bes3_D0, Xi2_part_Bes3_D0 = 0.0;
  double Vcs_single_Bes3_Dp, sigma_Vcs_single_Bes3_Dp, Xi2_Bes3_Dp, Xi2_part_Bes3_Dp = 0.0;
  double Vcs_single_Belle,   sigma_Vcs_single_Belle,   Xi2_Belle,   Xi2_part_Belle   = 0.0; 
  double Vcs_single_corr,    sigma_Vcs_single_corr,    Xi2_corr,    Xi2_tot          = 0.0;

	//// MATRICE DI COVARIANZA DI Vcs DOVUTA AI CONTRIBUTI DELLE DELTA GAMMA SPERIMENTALI
	for(int i = 0; i < Ntot; i++){
		for(int j = 0; j < Ntot; j++){
			Vcs_cov_exp[i*Ntot+j] = cov_Vckm_exp_func( Vcs_all, sqr_dGamma_all, sigma_sqr_dGamma_all, cov_exp_sqr_dGamma, i, j, Nev_tot, 0, Ncov);
		}// j
	}// i
 
  free(sqr_dGamma_all);
  free(sigma_sqr_dGamma_all);
  free(cov_exp_sqr_dGamma);

  for(int ic = 0; ic < Ncut; ic++){

		/////// COSTRUISCO LA MATRICE DI COVARIANZA SORT
		for(int i = 0; i < Ntot; i++){

			Vcs_all_sort[i] = Vcs_all[i*(Nev_tot+1)+0];
      
			for(int j = 0; j < Ntot; j++){

				Vcs_cov[i*Ntot+j] = Vcs_cov_lat_block[i*Ntot+j] + Vcs_cov_exp[i*Ntot+j];
	
				Vcs_cov_sort[i*Ntot+j] = Vcs_cov_lat[i*Ntot+j] + Vcs_cov_exp[i*Ntot+j];

			}// j
		}// i

#if defined(VCS_COV_LAT_TOT)
		////// SORT DI Vcs
		sort_vector( Vcs_all_sort, Ntot, q2_all);

		//// SORT DELLA MATRICE DI COVARIANZA
		sort_matrix( Vcs_cov_sort, Ntot, Ntot, q2_all, q2_all);

		//printf("\n\nMMMM Vcs_cov_sort \n\n");
		//print_matrix( Ntot, Ntot, Vcs_cov_sort, 1);
#endif  

		for(int i = 0; i < Ntot; i++){
			for(int j = 0; j < Ntot; j++){

				if( q2_all[i] > q2_cut[ic] || q2_all[j] > q2_cut[ic]){
	  
					if(i == j){
						Vcs_cov[i*Ntot+j] = 1.0;
					}else if(i != j){
						Vcs_cov[i*Ntot+j] = 0.0;
					}
	  
				}// if( q2_all[i] <= q2_cut[ic] && q2_all[j] <= q2_cut[ic])

			}// j
		}// i
    
		for(int i = 0; i < Ntot; i++){
			for(int j = 0; j < Ntot; j++){
	
				if( q2_all_sort[i] > q2_cut[ic] || q2_all_sort[j] > q2_cut[ic]){
	  
					if(i == j){
						Vcs_cov_sort[i*Ntot+j] = 1.0;
					}else if(i != j){
						Vcs_cov_sort[i*Ntot+j] = 0.0;
					}
	  
				}// if
	
			}// j
		}// i  

		/////////// MEDIA PESATA
		Vcs_single_BaBar   = Vckm_single_fit_with_cut( Vcs_all, Vcs_cov, q2_all, q2_cut, Ntot, Nev_tot,       0.0, Nsubgrp_1-1, ic);
		Vcs_single_Cleo_D0 = Vckm_single_fit_with_cut( Vcs_all, Vcs_cov, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_1, Nsubgrp_2-1, ic);
		Vcs_single_Cleo_Dp = Vckm_single_fit_with_cut( Vcs_all, Vcs_cov, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_2, Nsubgrp_3-1, ic);
		Vcs_single_Bes3_D0 = Vckm_single_fit_with_cut( Vcs_all, Vcs_cov, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_3, Nsubgrp_4-1, ic);
		Vcs_single_Bes3_Dp = Vckm_single_fit_with_cut( Vcs_all, Vcs_cov, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_4, Nsubgrp_5-1, ic);
		Vcs_single_Belle   = Vckm_single_fit_with_cut( Vcs_all, Vcs_cov, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_5, Nsubgrp_6-1, ic);

		Vcs_single_corr    = Vckm_single_fit_with_cut_corr( Vcs_all_sort, Vcs_cov_sort, q2_all_sort, q2_cut, Ntot, 0.0,  Nbin_max-1, ic);


		sigma_Vcs_single_BaBar   = sigma_Vckm_single_fit_with_cut( Vcs_cov, q2_all, q2_cut, Ntot,       0.0, Nsubgrp_1-1, ic);
		sigma_Vcs_single_Cleo_D0 = sigma_Vckm_single_fit_with_cut( Vcs_cov, q2_all, q2_cut, Ntot, Nsubgrp_1, Nsubgrp_2-1, ic);
		sigma_Vcs_single_Cleo_Dp = sigma_Vckm_single_fit_with_cut( Vcs_cov, q2_all, q2_cut, Ntot, Nsubgrp_2, Nsubgrp_3-1, ic);
		sigma_Vcs_single_Bes3_D0 = sigma_Vckm_single_fit_with_cut( Vcs_cov, q2_all, q2_cut, Ntot, Nsubgrp_3, Nsubgrp_4-1, ic);
		sigma_Vcs_single_Bes3_Dp = sigma_Vckm_single_fit_with_cut( Vcs_cov, q2_all, q2_cut, Ntot, Nsubgrp_4, Nsubgrp_5-1, ic);
		sigma_Vcs_single_Belle   = sigma_Vckm_single_fit_with_cut( Vcs_cov, q2_all, q2_cut, Ntot, Nsubgrp_5, Nsubgrp_6-1, ic);

		sigma_Vcs_single_corr    = sigma_Vckm_single_fit_with_cut( Vcs_cov_sort, q2_all_sort, q2_cut, Ntot, 0.0, Nbin_max-1, ic);

		/////////// CALCOLO I CHI2
		Xi2_BaBar   = Xi2_with_cut( Vcs_single_BaBar,   Vcs_all, Vcs_cov, q2_all, q2_cut, Ntot, Nev_tot,       0.0, Nsubgrp_1-1, ic);
		Xi2_Cleo_D0 = Xi2_with_cut( Vcs_single_Cleo_D0, Vcs_all, Vcs_cov, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_1, Nsubgrp_2-1, ic);
		Xi2_Cleo_Dp = Xi2_with_cut( Vcs_single_Cleo_Dp, Vcs_all, Vcs_cov, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_2, Nsubgrp_3-1, ic);
		Xi2_Bes3_D0 = Xi2_with_cut( Vcs_single_Bes3_D0, Vcs_all, Vcs_cov, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_3, Nsubgrp_4-1, ic);
		Xi2_Bes3_Dp = Xi2_with_cut( Vcs_single_Bes3_Dp, Vcs_all, Vcs_cov, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_4, Nsubgrp_5-1, ic);
		Xi2_Belle   = Xi2_with_cut( Vcs_single_Belle,   Vcs_all, Vcs_cov, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_5, Nsubgrp_6-1, ic);

		Xi2_corr    = Xi2_with_cut_corr( Vcs_single_corr, Vcs_all_sort, Vcs_cov_sort, q2_all_sort, q2_cut, Ntot, 0.0, Nbin_max-1, ic);

		///// CONTRIBUTI PARZIALI DEGLI ESPERIMENTI AL Xi2 TOTALE
#if defined(VCS_COV_LAT_BLOCK) || defined(NO_VCS_COV_LAT) || defined(VCS_COV_LAT_DIAG)

		if(ic == Ncut-1){
			Xi2_part_BaBar   = (Xi2_with_cut( Vcs_single_corr, Vcs_all, Vcs_cov, q2_all, q2_cut, Ntot, Nev_tot,       0.0, Nsubgrp_1-1, ic)*  (NBaBar-1.0))/(Nbin_max-1.0);
			Xi2_part_Cleo_D0 = (Xi2_with_cut( Vcs_single_corr, Vcs_all, Vcs_cov, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_1, Nsubgrp_2-1, ic)*(NCleo_D0-1.0))/(Nbin_max-1.0);
			Xi2_part_Cleo_Dp = (Xi2_with_cut( Vcs_single_corr, Vcs_all, Vcs_cov, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_2, Nsubgrp_3-1, ic)*(NCleo_Dp-1.0))/(Nbin_max-1.0);
			Xi2_part_Bes3_D0 = (Xi2_with_cut( Vcs_single_corr, Vcs_all, Vcs_cov, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_3, Nsubgrp_4-1, ic)*(NBes3_D0-1.0))/(Nbin_max-1.0);
			Xi2_part_Bes3_Dp = (Xi2_with_cut( Vcs_single_corr, Vcs_all, Vcs_cov, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_4, Nsubgrp_5-1, ic)*(NBes3_Dp-1.0))/(Nbin_max-1.0);

		#if defined(BELLE)
			Xi2_part_Belle   = (Xi2_with_cut( Vcs_single_corr, Vcs_all, Vcs_cov, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_5, Nsubgrp_6-1, ic)*(NBelle-1.0))/(Nbin_max-1.0);
		#endif
		#if defined(NO_BELLE)
			Xi2_part_Belle   = 0.0;
		#endif

			Xi2_tot = Xi2_part_BaBar + Xi2_part_Cleo_D0 + Xi2_part_Cleo_Dp + Xi2_part_Bes3_D0 + Xi2_part_Bes3_Dp + Xi2_part_Belle;

		}// if(ic == Ncut-1)
#endif
  /////////// FINE CALCOLO I CHI2

		printf("\n\nVCKM ic=%d q2_cut=%.3f\n", ic, q2_cut[ic]);

		printf("Vcs_BaBar   = %f  %f  Xi2= %f\n", Vcs_single_BaBar,   sigma_Vcs_single_BaBar,   Xi2_BaBar);
		printf("Vcs_Cleo_D0 = %f  %f  Xi2= %f\n", Vcs_single_Cleo_D0, sigma_Vcs_single_Cleo_D0, Xi2_Cleo_D0);
		printf("Vcs_Cleo_Dp = %f  %f  Xi2= %f\n", Vcs_single_Cleo_Dp, sigma_Vcs_single_Cleo_Dp, Xi2_Cleo_Dp);    
		printf("Vcs_Bes3_D0 = %f  %f  Xi2= %f\n", Vcs_single_Bes3_D0, sigma_Vcs_single_Bes3_D0, Xi2_Bes3_D0);
		printf("Vcs_Bes3_Dp = %f  %f  Xi2= %f\n", Vcs_single_Bes3_Dp, sigma_Vcs_single_Bes3_Dp, Xi2_Bes3_Dp);
		printf("Vcs_Belle   = %f  %f  Xi2= %f\n", Vcs_single_Belle,   sigma_Vcs_single_Belle,   Xi2_Belle);  
		printf("Vcs_corr    = %f  %f  Xi2= %f\n", Vcs_single_corr,    sigma_Vcs_single_corr,    Xi2_corr);  
  
#if defined(VCS_COV_LAT_BLOCK) || defined(NO_VCS_COV_LAT) || defined(VCS_COV_LAT_DIAG)
		if(ic == Ncut-1){
			printf("Xi2_tot = %f\nXi2_frac_BaBar = %f\nXi2_frac_Cleo_D0 = %f\nXi2_frac_Cleo_Dp = %f\nXi2_frac_Bes3_D0 = %f\nXi2_frac_Bes3_Dp = %f\nXi2_frac_Belle = %f\n", Xi2_tot, Xi2_part_BaBar/Xi2_tot, Xi2_part_Cleo_D0/Xi2_tot, Xi2_part_Cleo_Dp/Xi2_tot, Xi2_part_Bes3_D0/Xi2_tot, Xi2_part_Bes3_Dp/Xi2_tot, Xi2_part_Belle/Xi2_tot);
		}
#endif

  }// ic
  /////////// FINE ESTRAZIONE Vcs DAL FIT SINGOLO

  free(q2_cut);
  free(q2_all);
  free(q2_all_sort);
  free(Vcs_all);
  free(Vcs_all_sort);
  free(Vcs_cov);
  free(Vcs_cov_sort);
  free(Vcs_cov_exp);
  free(Vcs_cov_lat);
  free(Vcs_cov_lat_block);


  return 0;
 
}



void chi2_cov( int &npar, double *deriv, double &f, double *par, int iflag){

  double vfit[Nmatrix-1];
  
  f = 0;
  
  for(int i = 0; i < Nmatrix-1; i++){
    if( i < Nq ){
      vfit[i] = f_fit[i] - fzero_func( q2_fit[i], par);
    }else if( i >= Nq ){
      vfit[i] = f_fit[i] - fplus_func( q2_fit[i], par);
    }
  }// i

  for(int i = 0; i < Nmatrix-1; i++){
    for(int j = 0; j < Nmatrix-1; j++){
      f = f + vfit[i]*inv_cov_fit[i][j]*vfit[j];
    }// j
  }// i

}// chi2_cov


double chi2_cov_num(double *parameters, int num_par, int block_par, int* dof){

  double vfit[Nmatrix-1], f = 0;
  int num_points = 0;
  
  for(int i = 0; i < Nmatrix-1; i++){
    if( i < Nq ){
      vfit[i] = f_fit[i] - fzero_func( q2_fit[i], parameters);
    }else if( i >= Nq ){
      vfit[i] = f_fit[i] - fplus_func( q2_fit[i], parameters);
    }
  }// i

  for(int i = 0; i < Nmatrix-1; i++){

    num_points = num_points + 1;
    
    for(int j = 0; j < Nmatrix-1; j++){

      f = f + vfit[i]*inv_cov_fit[i][j]*vfit[j];
      
    }// j
  }// i

  *dof = num_points - num_par - block_par;

  return f/(*dof);
  
}// chi2_cov_num


void chi2_boot( int &npar, double *deriv, double &f, double *par, int iflag){

  double vfit[Nmatrix-1];
  
  f = 0;
  
  for(int i = 0; i < Nmatrix-1; i++){

    if( i < Nq ){

      vfit[i] = f_fit[i] - fzero_func( q2_fit[i], par);

    }else if( i >= Nq ){

      vfit[i] = f_fit[i] - fplus_func( q2_fit[i], par);
      
    }

    
  }// i


  for(int i = 0; i < Nmatrix-1; i++){
    
      f = f + pow( vfit[i] , 2)/pow( sigma_f_fit[i], 2);

  }// i

  
}// chi2_boot


double chi2_boot_num(double *parameters, int num_par, int block_par, int* dof){

  double vfit[Nmatrix-1];
  int num_points = 0;
  
  double f = 0;
  
  for(int i = 0; i < Nmatrix-1; i++){

    if( i < Nq ){

      vfit[i] = f_fit[i] - fzero_func( q2_fit[i], parameters);

    }else if( i >= Nq ){

      vfit[i] = f_fit[i] - fplus_func( q2_fit[i], parameters);
      
    }

    
  }// i


  for(int i = 0; i < Nmatrix-1; i++){
    
    f = f + pow( vfit[i] , 2)/pow( sigma_f_fit[i], 2);

    num_points = num_points + 1;
      
  }// i
      

  *dof = num_points - num_par - block_par;

  return f/(*dof);

  
}// chi2_boot_num






double fzero_func( double q2, double *par){

  double res, A_fzero, z, z_term;
  double f0;
  
  res = par[0];
  A_fzero = par[1];
  
  z = ( sqrt( tp - q2 ) - sqrt( tp - t0 ) )/( sqrt( tp - q2 ) + sqrt( tp - t0 ) );
  z_term = ( z - z0 )*( 1 + ( z + z0 )/2  );
  
  f0 = ( res + A_fzero*z_term );

  return f0;

}


double fplus_func( double q2, double *par){

  double res, A_fplus, polo_fplus = 0.224188, z, z_term; // fissato alla massa del Ds*
  double fp;
  
  res = par[0];
  A_fplus = par[2];
  
  z = ( sqrt( tp - q2 ) - sqrt( tp - t0 ) )/( sqrt( tp - q2 ) + sqrt( tp - t0 ) );
  z_term = ( z - z0 )*( 1 + ( z + z0 )/2  );
  
  fp = ( res + A_fplus*z_term )/( 1 - polo_fplus*q2  );

  return fp;

}


double p_cube( double q2, int isospin){

  double p3;

  if(isospin == 0){
    p3 = pow( pow( ( pow(MD0_pdg,2) + pow(MKm_pdg,2) - q2 )/(2*MD0_pdg) , 2) - pow(MKm_pdg,2) , 1.5);
  }else if(isospin == 1){
    p3 = pow( pow( ( pow(MDp_pdg,2) + pow(MK0_pdg,2) - q2 )/(2*MDp_pdg) , 2) - pow(MK0_pdg,2) , 1.5);
  }
  
  return p3;

}


double bin_integral( double a, double b, int nbin, double *par, int isospin){

  double inte = 0;
  double h, i;
   
  h=(b-a)/nbin;
  
  for(i = a; i <= b-h; i = i+h){
    inte = inte + h*( p_cube(i, isospin)*pow(fplus_func( i, par) ,2) + p_cube(i+h, isospin)*pow(fplus_func( i+h, par) ,2) )/2.;
  }

  return inte;
  
}



double Vcs_func( double a, double b, double *par, double sqr_dGamma, int isospin){

  double C, V;

  C = sqrt( (24.*pow(PI,3))/(pow(GF,2)*bin_integral( a, b, nbin, par, isospin)) );
  
  if(isospin == 0){
    V = C*sqr_dGamma;
  }else if(isospin == 1){
    V = C*sqr_dGamma;
  }
  
  return V;
  
}






void read_fpVcs( double* q2_BaBar, double* q2_Cleo_D0, double* q2_Cleo_Dp, double* q2_Belle, double* fpVcs_BaBar, double* fpVcs_Cleo_D0, double* fpVcs_Cleo_Dp, double* fpVcs_Belle, double* sigma_fpVcs_BaBar, double* sigma_fpVcs_Cleo_D0, double* sigma_fpVcs_Cleo_Dp, double* sigma_fpVcs_Belle, int NBaBar, int NCleo_D0, int NCleo_Dp, int NBelle){

  FILE *fr_fpVcs;
  int t;

  if ((fr_fpVcs = fopen( "Input_corr_3pts/fpVcs_sperimentale_DK.out", "r")) == NULL ){
    printf("Error opening the file to read: open_fpVcs\n");
    exit(EXIT_FAILURE);
  }

  for(int i = 0; i < NBaBar; i++){
    t = fscanf(fr_fpVcs,"%lf %lf %lf\n", &q2_BaBar[i], &fpVcs_BaBar[i], &sigma_fpVcs_BaBar[i]);
  }
  for(int i = 0; i < NCleo_D0; i++){
    t = fscanf(fr_fpVcs,"%lf %lf %lf\n", &q2_Cleo_D0[i], &fpVcs_Cleo_D0[i], &sigma_fpVcs_Cleo_D0[i]); 
  }
  for(int i = 0; i < NCleo_Dp; i++){
    t = fscanf(fr_fpVcs,"%lf %lf %lf\n", &q2_Cleo_Dp[i], &fpVcs_Cleo_Dp[i], &sigma_fpVcs_Cleo_Dp[i]);
  }
  for(int i = 0; i < NBelle; i++){
    t = fscanf(fr_fpVcs,"%lf %lf %lf\n", &q2_Belle[i], &fpVcs_Belle[i], &sigma_fpVcs_Belle[i]);
  }
  
  fclose(fr_fpVcs);
}// read_fpVcs



void read_sqr_dGamma( double* q2_inf_BaBar, double* q2_inf_Cleo_D0, double* q2_inf_Cleo_Dp, double* q2_inf_Bes3_D0, double* q2_inf_Bes3_Dp, double* q2_sup_BaBar, double* q2_sup_Cleo_D0, double* q2_sup_Cleo_Dp, double* q2_sup_Bes3_D0, double* q2_sup_Bes3_Dp, double* q2_BaBar, double* q2_Cleo_D0, double* q2_Cleo_Dp, double* q2_Bes3_D0, double* q2_Bes3_Dp, double* sqr_dGamma_BaBar, double* sqr_dGamma_Cleo_D0, double* sqr_dGamma_Cleo_Dp, double* sqr_dGamma_Bes3_D0, double* sqr_dGamma_Bes3_Dp, double* sigma_sqr_dGamma_BaBar, double* sigma_sqr_dGamma_Cleo_D0, double* sigma_sqr_dGamma_Cleo_Dp, double* sigma_sqr_dGamma_Bes3_D0, double* sigma_sqr_dGamma_Bes3_Dp, int NBaBar, int NCleo_D0, int NCleo_Dp, int NBes3_D0, int NBes3_Dp){

  FILE *fr_sqr_dGamma;
  int t;
  
  if ((fr_sqr_dGamma = fopen("Input_corr_3pts/sqrt_delta_gamma_sperimentali_D_to_K.out", "r")) == NULL ){
    printf("Error opening the file to read: open_sqr_dGamma\n");
    exit(EXIT_FAILURE);
  }

  for(int i = 0; i < NBaBar; i++){
    t = fscanf(fr_sqr_dGamma,"%lf %lf %lf %lf %lf\n", &q2_BaBar[i], &q2_inf_BaBar[i], &q2_sup_BaBar[i], &sqr_dGamma_BaBar[i], &sigma_sqr_dGamma_BaBar[i]);
  }
  for(int i = 0; i < NCleo_D0; i++){
    t = fscanf(fr_sqr_dGamma,"%lf %lf %lf %lf %lf\n", &q2_Cleo_D0[i], &q2_inf_Cleo_D0[i], &q2_sup_Cleo_D0[i], &sqr_dGamma_Cleo_D0[i], &sigma_sqr_dGamma_Cleo_D0[i]); 
  }
  for(int i = 0; i < NCleo_Dp; i++){
    t = fscanf(fr_sqr_dGamma,"%lf %lf %lf %lf %lf\n", &q2_Cleo_Dp[i], &q2_inf_Cleo_Dp[i], &q2_sup_Cleo_Dp[i], &sqr_dGamma_Cleo_Dp[i], &sigma_sqr_dGamma_Cleo_Dp[i]);
  }
  for(int i = 0; i < NBes3_D0; i++){
    t = fscanf(fr_sqr_dGamma,"%lf %lf %lf %lf %lf\n", &q2_Bes3_D0[i], &q2_inf_Bes3_D0[i], &q2_sup_Bes3_D0[i], &sqr_dGamma_Bes3_D0[i], &sigma_sqr_dGamma_Bes3_D0[i]);
  }
  for(int i = 0; i < NBes3_Dp; i++){
    t = fscanf(fr_sqr_dGamma,"%lf %lf %lf %lf %lf\n", &q2_Bes3_Dp[i], &q2_inf_Bes3_Dp[i], &q2_sup_Bes3_Dp[i], &sqr_dGamma_Bes3_Dp[i], &sigma_sqr_dGamma_Bes3_Dp[i]);
  }
  fclose(fr_sqr_dGamma);

}// read_sqr_dGamma


void check_read_data( double* q2_inf_BaBar, double* q2_inf_Cleo_D0, double* q2_inf_Cleo_Dp, double* q2_inf_Bes3_D0, double* q2_inf_Bes3_Dp, double* q2_sup_BaBar, double* q2_sup_Cleo_D0, double* q2_sup_Cleo_Dp, double* q2_sup_Bes3_D0, double* q2_sup_Bes3_Dp, double* q2_BaBar, double* q2_Cleo_D0, double* q2_Cleo_Dp, double* q2_Bes3_D0, double* q2_Bes3_Dp, double* q2_Belle, double* sqr_dGamma_BaBar, double* sqr_dGamma_Cleo_D0, double* sqr_dGamma_Cleo_Dp, double* sqr_dGamma_Bes3_D0, double* sqr_dGamma_Bes3_Dp, double* sigma_sqr_dGamma_BaBar, double* sigma_sqr_dGamma_Cleo_D0, double* sigma_sqr_dGamma_Cleo_Dp, double* sigma_sqr_dGamma_Bes3_D0, double* sigma_sqr_dGamma_Bes3_Dp, double* fpVcs_BaBar, double* fpVcs_Cleo_D0, double* fpVcs_Cleo_Dp, double* fpVcs_Belle, double* sigma_fpVcs_BaBar, double* sigma_fpVcs_Cleo_D0, double* sigma_fpVcs_Cleo_Dp, double* sigma_fpVcs_Belle, int NBaBar, int NCleo_D0, int NCleo_Dp, int NBes3_D0, int NBes3_Dp, int NBelle){

  //BABAR
  for(int i = 0; i < NBaBar; i++){
    printf("%f %f %f %f %f %.12f %.12f\n", q2_BaBar[i], q2_inf_BaBar[i], q2_sup_BaBar[i], fpVcs_BaBar[i], sigma_fpVcs_BaBar[i], sqr_dGamma_BaBar[i], sigma_sqr_dGamma_BaBar[i]);
  }
  printf("&\n");

  //CLEO-D0  
  for(int i = 0; i < NCleo_D0; i++){
    printf("%f %f %f %f %f %.12f %.12f\n", q2_Cleo_D0[i], q2_inf_Cleo_D0[i], q2_sup_Cleo_D0[i], fpVcs_Cleo_D0[i], sigma_fpVcs_Cleo_D0[i], sqr_dGamma_Cleo_D0[i], sigma_sqr_dGamma_Cleo_D0[i]);
  }
  printf("&\n");
  
  //CLEO-Dp
  for(int i = 0; i < NCleo_Dp; i++){
    printf("%f %f %f %f %f %.12f %.12f\n", q2_Cleo_Dp[i], q2_inf_Cleo_Dp[i], q2_sup_Cleo_Dp[i], fpVcs_Cleo_Dp[i], sigma_fpVcs_Cleo_Dp[i], sqr_dGamma_Cleo_Dp[i], sigma_sqr_dGamma_Cleo_Dp[i]);
  }
  printf("&\n");
  
  //BES3-D0
  for(int i = 0; i < NBes3_D0; i++){
    printf("%f %f %f %.12f %.12f\n", q2_Bes3_D0[i], q2_inf_Bes3_D0[i], q2_sup_Bes3_D0[i], sqr_dGamma_Bes3_D0[i], sigma_sqr_dGamma_Bes3_D0[i]);
  }
  printf("&\n");

  //BES3-D+
  for(int i = 0; i < NBes3_Dp; i++){
    printf("%f %f %f %.12f %.12f\n", q2_Bes3_Dp[i], q2_inf_Bes3_Dp[i], q2_sup_Bes3_Dp[i], sqr_dGamma_Bes3_Dp[i], sigma_sqr_dGamma_Bes3_Dp[i]);
  }
  printf("&\n");
  
  //BELLE
  for(int i = 0; i < NBelle; i++){
    printf("%f %f %f\n", q2_Belle[i], fpVcs_Belle[i], sigma_fpVcs_Belle[i]);
  }

}// check_read_data


void Vckm_from_sqrdGamma_boot( double* Vckm, int Nev, int Nexp, double* q2_inf, double* q2_sup, double par_boot[Nev_tot+1][npar], double* sqr_dGamma, int isospin, double (*Vckm_func)( double, double, double*, double, int)){

  for(int i = 0; i < Nexp; i++){

    Vckm[i*(Nev+1)+0] = 0.0;
    
    for(int iev = 1; iev <= Nev; iev++){
      
      Vckm[i*(Nev+1)+iev] = Vckm_func( q2_inf[i], q2_sup[i], par_boot[iev], sqr_dGamma[i], isospin);

      Vckm[i*(Nev+1)+0] += Vckm[i*(Nev+1)+iev]/Nev;
 		
    }// iev
  }// i

}// Vckm_from_sqrdGamma_boot


void Vckm_from_fpVckm_boot( double* Vckm, int Nev, int Nexp, double* q2, double par_boot[Nev_tot+1][npar], double* fpVckm, double (*fplus_func)( double, double*)){

  double **fp_synt=(double**)malloc(sizeof(double)*(Nexp));
  for(int i = 0; i < Nexp; i++) fp_synt[i]=(double*)malloc(sizeof(double)*(Nev+1));

  for(int i = 0; i < Nexp; i++){
    fp_synt[i][0] = 0.0;
    Vckm[i*(Nev+1)+0] = 0.0;

    for(int iev = 1; iev <= Nev; iev++){
    
      fp_synt[i][iev] = fplus_func( q2[i], par_boot[iev]);      
      fp_synt[i][0] += fp_synt[i][iev]/Nev;

      Vckm[i*(Nev+1)+iev] = fpVckm[i]/fp_synt[i][iev];
      Vckm[i*(Nev+1)+0] += Vckm[i*(Nev+1)+iev]/Nev;
    }// iev
  }// i

  free(fp_synt);
}// Vckm_from_fpVckm_boot


void sigma_Vckm_from_sqrdGamma_boot( double* sigma_Vckm, double* Vckm, int Nev, int Nev_an, int Nbin, double* sqr_dGamma, double* sigma_sqr_dGamma, int analysis_in, int analysis_fin, int clusterfile){

	double *temp_Vckm=(double*)calloc( (Nbin)*(Nev+1), sizeof(double));
  double *sigma_Vckm_sqr_dGamma=(double*)malloc(sizeof(double)*(Nbin)), *sigma_Vckm_lattice=(double*)malloc(sizeof(double)*(Nbin));

  for(int i = 0; i < Nbin; i++){

		for(int iev = 1; iev <= Nev; iev++){
			temp_Vckm[iev] = Vckm[i*(Nev+1)+iev];
      temp_Vckm[0] += temp_Vckm[iev]/Nev; 
    }// iev

		sigma_Vckm_sqr_dGamma[i] = (temp_Vckm[0]/sqr_dGamma[i])*sigma_sqr_dGamma[i];

    sigma_Vckm_lattice[i] = sigma_bootstrap( temp_Vckm, analysis_in, analysis_fin, Nev_an, clusterfile);

    sigma_Vckm[i] = sqrt( pow(sigma_Vckm_sqr_dGamma[i], 2) + pow(sigma_Vckm_lattice[i], 2) );

		temp_Vckm[0] = 0.0; 

  }// i (bin)

	free(temp_Vckm);
	free(sigma_Vckm_sqr_dGamma);
	free(sigma_Vckm_lattice);

}// sigma_Vckm_from_sqrdGamma_boot


double cov_Vckm_exp_func( double* Vckm, double* sqr_dGamma, double* sigma_sqr_dGamma, double* cov_exp, int ibin_1, int ibin_2, int Nev, int iev, int Ncov){

      if( ibin_1 < Ncov && ibin_2 < Ncov ){
	
        double c_i, c_j;

	  		c_i = Vckm[ibin_1*(Nev+1)+iev]/sqr_dGamma[ibin_1];
				c_j = Vckm[ibin_2*(Nev+1)+iev]/sqr_dGamma[ibin_2];

        return (c_i)*(c_j)*cov_exp[ibin_1*Ncov+ibin_2];

      }else if( ibin_1 >= Ncov && ibin_2 <  Ncov ){

				return 0;

      }else if( ibin_1 <  Ncov && ibin_2 >= Ncov ){

        return 0;

      }else if( ibin_1 >= Ncov && ibin_2 >= Ncov && ibin_1 == ibin_2){

        return pow( (Vckm[ibin_1*(Nev+1)+iev]/sqr_dGamma[ibin_1])*sigma_sqr_dGamma[ibin_1], 2);
	  
      }else if( ibin_1 >= Ncov && ibin_2 >= Ncov && ibin_1 != ibin_2){

				return 0;

			}// if

}// cov_Vckm_exp_func


double Vckm_const_fit_with_cut( double* Vckm, double* Vckm_cov, double* q2_all, double* q2_cut, int Ntot, int Nev_tot, int i_start, int i_end, int iev, int ic){

  double num = 0.0, den = 0.0, *inv_Vckm_cov=(double*)malloc(sizeof(double)*(Ntot)*(Ntot));

	LU_invert( inv_Vckm_cov, Vckm_cov, Ntot);    

	for(int i = i_start; i < i_end; i++){
		for(int j = i_start; j < i_end; j++){

			if( q2_all[i] <= q2_cut[ic] && q2_all[j] <= q2_cut[ic]){
	  
				num += inv_Vckm_cov[i*Ntot+j]*Vckm[j*(Nev_tot+1)+iev];
				den += inv_Vckm_cov[i*Ntot+j];
	  
			}// if
		}// j
	}// i

  free(inv_Vckm_cov);

	return num/den;

}// Vckm_const_fit_with_cut



double Vckm_const_fit_with_cut_corr( double* Vckm, double* Vckm_cov, double* q2_all, double* q2_cut, int Ntot, int i_start, int i_end, int ic){

  double num = 0.0, den = 0.0, *inv_Vckm_cov=(double*)malloc(sizeof(double)*(Ntot)*(Ntot));

	LU_invert( inv_Vckm_cov, Vckm_cov, Ntot);    

	for(int i = i_start; i < i_end; i++){
		for(int j = i_start; j < i_end; j++){

			if( q2_all[i] <= q2_cut[ic] && q2_all[j] <= q2_cut[ic]){
	  
				num += inv_Vckm_cov[i*Ntot+j]*Vckm[j];
				den += inv_Vckm_cov[i*Ntot+j];
	  
			}// if
		}// j
	}// i

  free(inv_Vckm_cov);

	return num/den;

}// Vckm_const_fit_with_cut_corr 


double Vckm_single_fit_with_cut( double* Vckm, double* Vckm_cov, double* q2_all, double* q2_cut, int Ntot, int Nev_tot, int i_start, int i_end, int ic){

  double Vckm_fit = 0.0, *temp_coeff=(double*)calloc( 1, sizeof(double)), *inv_Vckm_cov=(double*)malloc(sizeof(double)*(Ntot)*(Ntot));

  double *coeff_average=(double*)calloc( i_end-i_start+1, sizeof(double));

	LU_invert( inv_Vckm_cov, Vckm_cov, Ntot);    

	for(int i = i_start; i <= i_end; i++){
		for(int j = i_start; j <= i_end; j++){
			if( q2_all[i] <= q2_cut[ic] && q2_all[j] <= q2_cut[ic]){
				*temp_coeff += inv_Vckm_cov[i*Ntot+j];
			}
		}// j
	}// i

	for(int i = i_start; i <= i_end; i++){
		for(int j = i_start; j <= i_end; j++){
			if( q2_all[i] <= q2_cut[ic] && q2_all[j] <= q2_cut[ic]){
				coeff_average[i-i_start] += inv_Vckm_cov[i*Ntot+j]/(*temp_coeff);
			}
		}// j
	}// i
	free(temp_coeff);
  free(inv_Vckm_cov);

  for(int i = i_start; i <= i_end; i++){
    if( q2_all[i] <= q2_cut[ic]){
      Vckm_fit += coeff_average[i-i_start]*Vckm[i*(Nev_tot+1)+0];
    }
  }// i
  free(coeff_average);

  return Vckm_fit;

}// Vckm_single_fit_with_cut




double Vckm_single_fit_with_cut_corr( double* Vckm, double* Vckm_cov, double* q2_all, double* q2_cut, int Ntot, int i_start, int i_end, int ic){

  double Vckm_fit = 0.0, *temp_coeff=(double*)calloc( 1, sizeof(double)), *inv_Vckm_cov=(double*)malloc(sizeof(double)*(Ntot)*(Ntot));

  double *coeff_average=(double*)calloc( i_end-i_start+1, sizeof(double));

	LU_invert( inv_Vckm_cov, Vckm_cov, Ntot);    

	for(int i = i_start; i <= i_end; i++){
		for(int j = i_start; j <= i_end; j++){
			if( q2_all[i] <= q2_cut[ic] && q2_all[j] <= q2_cut[ic]){
				*temp_coeff += inv_Vckm_cov[i*Ntot+j];
			}
		}// j
	}// i

	for(int i = i_start; i <= i_end; i++){
		for(int j = i_start; j <= i_end; j++){
			if( q2_all[i] <= q2_cut[ic] && q2_all[j] <= q2_cut[ic]){
				coeff_average[i-i_start] += inv_Vckm_cov[i*Ntot+j]/(*temp_coeff);
			}
		}// j
	}// i
	free(temp_coeff);
  free(inv_Vckm_cov);

  for(int i = i_start; i <= i_end; i++){
    if( q2_all[i] <= q2_cut[ic]){
      Vckm_fit += coeff_average[i-i_start]*Vckm[i];
    }
  }// i
  free(coeff_average);

  return Vckm_fit;

}// Vckm_single_fit_with_cut_corr


double sigma_Vckm_single_fit_with_cut( double* Vckm_cov, double* q2_all, double* q2_cut, int Ntot, int i_start, int i_end, int ic){

  double var_Vckm_fit = 0.0, *temp_coeff=(double*)calloc( 1, sizeof(double)), *inv_Vckm_cov=(double*)malloc(sizeof(double)*(Ntot)*(Ntot));
  double *coeff_average=(double*)calloc( i_end-i_start+1, sizeof(double));

	LU_invert( inv_Vckm_cov, Vckm_cov, Ntot);    

	for(int i = i_start; i <= i_end; i++){
		for(int j = i_start; j <= i_end; j++){
			if( q2_all[i] <= q2_cut[ic] && q2_all[j] <= q2_cut[ic]){
				*temp_coeff += inv_Vckm_cov[i*Ntot+j];
			}
		}// j
	}// i

	for(int i = i_start; i <= i_end; i++){
		for(int j = i_start; j <= i_end; j++){
			if( q2_all[i] <= q2_cut[ic] && q2_all[j] <= q2_cut[ic]){
				coeff_average[i-i_start] += inv_Vckm_cov[i*Ntot+j]/(*temp_coeff);
			}
		}// j
	}// i
	free(temp_coeff);
  free(inv_Vckm_cov);

  for(int i = i_start; i <= i_end; i++){
    for(int j = i_start; j <= i_end; j++){
      if( q2_all[i] <= q2_cut[ic] && q2_all[j] <= q2_cut[ic]){
				var_Vckm_fit += coeff_average[i-i_start]*coeff_average[j-i_start]*Vckm_cov[i*Ntot+j];
      }
    }// j
  }// i
  free(coeff_average);

  return sqrt(var_Vckm_fit);

}// sigma_Vckm_single_fit_with_cut


double Xi2_with_cut( double Vckm_av, double* Vckm, double* Vckm_cov, double* q2_all, double* q2_cut, int Ntot, int Nev_tot, int i_start, int i_end, int ic){

  int num_pts = 0;
  double Xi2 =0.0, *inv_Vckm_cov=(double*)malloc(sizeof(double)*(Ntot)*(Ntot));

	LU_invert( inv_Vckm_cov, Vckm_cov, Ntot);    

  for(int i = i_start; i <= i_end; i++){
    for(int j = i_start; j <= i_end; j++){

      if( q2_all[i] <= q2_cut[ic] && q2_all[j] <= q2_cut[ic]){
				Xi2 += (Vckm_av - Vckm[i*(Nev_tot+1)+0])*inv_Vckm_cov[i*Ntot+j]*(Vckm_av - Vckm[j*(Nev_tot+1)+0]);
      }
    }// j
    if( q2_all[i] <= q2_cut[ic]){
      num_pts += 1;
    }
  }// i
  free(inv_Vckm_cov);

  return Xi2/(num_pts-1);

}// Xi2_with_cut


double Xi2_with_cut_corr( double Vckm_av, double* Vckm, double* Vckm_cov, double* q2_all, double* q2_cut, int Ntot, int i_start, int i_end, int ic){

  int num_pts = 0;
  double Xi2 =0.0, *inv_Vckm_cov=(double*)malloc(sizeof(double)*(Ntot)*(Ntot));

	LU_invert( inv_Vckm_cov, Vckm_cov, Ntot);    

  for(int i = i_start; i <= i_end; i++){
    for(int j = i_start; j <= i_end; j++){

      if( q2_all[i] <= q2_cut[ic] && q2_all[j] <= q2_cut[ic]){
				Xi2 += (Vckm_av - Vckm[i])*inv_Vckm_cov[i*Ntot+j]*(Vckm_av - Vckm[j]);
      }
    }// j
    if( q2_all[i] <= q2_cut[ic]){
      num_pts += 1;
    }
  }// i
  free(inv_Vckm_cov);

  return Xi2/(num_pts-1);

}// Xi2_with_cut_corr


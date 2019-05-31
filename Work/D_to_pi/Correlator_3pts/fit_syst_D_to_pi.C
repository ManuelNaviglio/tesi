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
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include "definitions.h"

#define LEN_NAME 1024
#define PI 3.141592653589793
#define GF 0.000011664

using namespace std;

struct params_ffs
{
  double res;
  double A_fzero;
  double A_fplus;
  double polo_fzero;
  double polo_fplus;
  params_ffs( double res, double A_fzero, double A_fplus, double polo_fzero, double polo_fplus) : res(res), A_fzero(A_fzero), A_fplus(A_fplus), polo_fzero(polo_fzero), polo_fplus(polo_fplus) {}
};

struct params_Vckm
{
  double res;
  double A_fzero;
  double A_fplus;
  double polo_fzero;
  double polo_fplus;
  double MH;
  double Ml;
  params_Vckm( double res, double A_fzero, double A_fplus, double polo_fzero, double polo_fplus, double MH, double Ml) : res(res), A_fzero(A_fzero), A_fplus(A_fplus), polo_fzero(polo_fzero), polo_fplus(polo_fplus), MH(MH), Ml(Ml) {}
};

struct params_Gamma
{
  double res;
  double A_fzero;
  double A_fplus;
  double polo_fzero;
  double polo_fplus;
  double MH;
  double Ml;
  double mell;
  params_Gamma( double res, double A_fzero, double A_fplus, double polo_fzero, double polo_fplus, double MH, double Ml, double mell) : res(res), A_fzero(A_fzero), A_fplus(A_fplus), polo_fzero(polo_fzero), polo_fplus(polo_fplus), MH(MH), Ml(Ml), mell(mell) {}
};

const int Nq = 8;
const int Nmatrix = 16;
const int Nrand = 500;      // UGUALI
const int Nev_sint = 500;   // UGUALI
const int Nev_tot = 6400;
int analysis_in = 1, analysis_fin = 64;

string form_factors[2] = {"fzero", "fplus"};

const int npar = 5;
int nbin = 1000;

int output = 1; // IL COMANDO PER ACCENDERE L'OUTPUT

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

double tp = pow( MD_phys + Mpi_phys , 2), t0 = ( MD_phys +  Mpi_phys)*pow( sqrt( MD_phys ) - sqrt( Mpi_phys ) , 2), z0 = ( sqrt( tp ) - sqrt( tp - t0 ) )/( sqrt( tp ) + sqrt( tp - t0 ) );

double MD0_pdg = 1.865, Mpim_pdg = 0.140, MDp_pdg = 1.870, Mpi0_pdg = 0.135;

//////////////////////

////////// QUANTITÀ GLOBALI PER IL FIT

double f_fit[Nmatrix-1], q2_fit[Nmatrix-1], inv_cov_fit[Nmatrix-1][Nmatrix-1], sigma_f_fit[Nmatrix-1];

//////////////////////////////////////

void chi2_cov( int &, double *, double &, double *, int);
void chi2_boot( int &, double *, double &, double *, int);
double chi2_cov_num(double *, int, int, int*);
double chi2_boot_num(double *, int, int, int*);

double fzero_func( double, void *);
double fplus_func( double, void *);
double p_cube( double, void*);
double dGamma_without_const( double, void*);

void read_fpVcd( double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, int, int, int, int);
void read_sqr_dGamma( double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, int, int, int, int, int);
void check_read_data( double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, int, int, int, int, int, int);

void Vckm_from_sqrdGamma_boot_new( double*, int, int, double*, double*, double [Nev_tot+1][npar], double*, double*, int, int, int, int, double (*)( double, void*));
void Vckm_from_fpVckm_boot_new( double*, int, int, double*, double [Nev_tot+1][npar], double*, double*, int, int, double (*)( double, void*));

// NON USATA
double cov_Vckm_exp_func( double*, double*, double*, double*, int, int, int, int, int);

double Vckm_const_fit_with_cut( double*, double*, double*, double*, int, int, int, int, int, int);
double Vckm_const_fit_with_cut_corr( double*, double*, double*, double*, int, int, int, int);
double Vckm_const_fit_with_cut_reverse( double*, double*, double*, double*, int, int, int, int, int, int);
double Vckm_const_fit_with_cut_corr_reverse( double*, double*, double*, double*, int, int, int, int);

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
  double q2_max = 3.0, q2_min = 0.0;
  //double q2_step;

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
  *var_f0p = 5.0*pow(10, -6);
  
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
  double step[npar] = {0.01, 0.01, 0.01, 0.01, 0.01};
  double par[npar] = {0.60, 0.00, 0.00, 0.03, 0.1};      ///// CONTROLLA QUESTO
  double min[npar] = {0.00, 0.00, 0.00, 0.00, 0.00};
  double max[npar] = {0.00, 0.00, 0.00, 0.00, 0.00};
  std::string cpar[npar] = {"f(0)", "C0", "C+", "P0", "P+"};
   
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
  
  covariance_syntetic_data( boot_cov_par, npar, Nev_tot, Nev_an, *var_boot_par, temp_boot_par, analysis_in, analysis_fin, clusterfile, 0);

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

  double **fzero_plot=(double**)malloc(sizeof(double)*(Nstep));
  for(int iq = 0; iq < Nstep; iq++) fzero_plot[iq]=(double*)malloc(sizeof(double)*(Nev_sint+1));
  
  double **fplus_plot=(double**)malloc(sizeof(double)*(Nstep));
  for(int iq = 0; iq < Nstep; iq++) fplus_plot[iq]=(double*)malloc(sizeof(double)*(Nev_sint+1));

  double *temp_fzero_plot=(double*)malloc(sizeof(double)*(Nstep)), *temp_fplus_plot=(double*)malloc(sizeof(double)*(Nstep));
  double *sigma_fzero_plot=(double*)malloc(sizeof(double)*(Nstep)), *sigma_fplus_plot=(double*)malloc(sizeof(double)*(Nstep));
  
  double q2_step;

  params_ffs params( 1.0, 2.0, 3.0, 4.0, 5.0);

  gsl_vector* mu_new = gsl_vector_alloc(npar);
  gsl_vector* random_result = gsl_vector_alloc(npar);
  gsl_matrix* boot_cov_par_new = gsl_matrix_alloc( npar, npar);

  for(int ipar = 0; ipar < npar; ipar++){
    gsl_vector_set( mu_new, ipar, mu[ipar]);
  }// ipar
  //gsl_vector_fprintf( stdout, mu_new, "%+.3e");

/*
  for(int irow = 0; irow < npar; irow++){
    for(int icol = 0; icol < npar; icol++){
      gsl_matrix_set( boot_cov_par_new, irow, icol, boot_cov_par[irow*npar + icol]);
    }// irow
  }// icol
	//gsl_matrix_fprintf( stdout, boot_cov_par_new, "%+.5e");
*/

  for(int irow = 0; irow < npar; irow++){
    for(int icol = 0; icol < npar; icol++){
      gsl_matrix_set( boot_cov_par_new, irow, icol, m_par[irow*npar + icol]);
    }// irow
  }// icol

	gsl_linalg_cholesky_decomp1( boot_cov_par_new);

	const gsl_rng_type* ran_T;
	gsl_rng* ran_gen;

	gsl_rng_env_setup();

	ran_T = gsl_rng_default;
	ran_gen = gsl_rng_alloc(ran_T);

	gsl_rng_set( ran_gen, time(NULL));

  for(int iq = 0; iq < Nstep; iq++){

    q2_step = (q2_max/( (double)  (Nstep - 1)))*iq;    

    for(int iev = 1; iev <= Nev_sint; iev++){

			gsl_ran_multivariate_gaussian( ran_gen, mu_new, boot_cov_par_new, random_result);

			params.res        = gsl_vector_get( random_result, 0);
			params.A_fzero    = gsl_vector_get( random_result, 1);
			params.A_fplus    = gsl_vector_get( random_result, 2);
			params.polo_fzero = gsl_vector_get( random_result, 3);
			params.polo_fplus = gsl_vector_get( random_result, 4);

      fzero_plot[iq][iev] = fzero_func( q2_step, &params);
      fplus_plot[iq][iev] = fplus_func( q2_step, &params);

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
  free(fzero_plot);
  free(fplus_plot);  
  free(temp_fzero_plot); 
  free(temp_fplus_plot); 
  free(sigma_fzero_plot); 
  free(sigma_fplus_plot);

  gsl_vector_free(mu_new);
  gsl_vector_free(random_result);
  gsl_matrix_free(boot_cov_par_new);

	gsl_rng_free(ran_gen);
  //// FINE GENERO LA DISTRIBUZIONE SISTETICA DEI FATTORI DI FORMA


  //////////////////////
  //                  //
  //   CALCOLO Vcd    //
  //                  //
  //////////////////////

  const int Nset_data = 6;

  const int NBaBar   = 10, Nsubgrp_1 = NBaBar;
  const int NCleo_D0 =  7, Nsubgrp_2 = NBaBar+NCleo_D0;
  const int NCleo_Dp =  7, Nsubgrp_3 = NBaBar+NCleo_D0+NCleo_Dp;
  const int NBes3_D0 = 14, Nsubgrp_4 = NBaBar+NCleo_D0+NCleo_Dp+NBes3_D0;
  const int NBes3_Dp =  7, Nsubgrp_5 = NBaBar+NCleo_D0+NCleo_Dp+NBes3_D0+NBes3_Dp;
  const int NBelle =   10, Nsubgrp_6 = NBaBar+NCleo_D0+NCleo_Dp+NBes3_D0+NBes3_Dp+NBelle;

  const int Ntot = NBaBar+NCleo_D0+NCleo_Dp+NBes3_D0+NBes3_Dp+NBelle, Ncov = NBaBar+NCleo_D0+NCleo_Dp+NBes3_D0+NBes3_Dp;
 
#if defined(BELLE)
	const int Nbin_max = Ntot;
#endif
#if defined(NO_BELLE)
	const int Nbin_max = Ncov;
#endif  
 
  ////////// LEGGO I DATI SPERIMENTALI DI q2, fpVcd, sqr_dGamma E LA MATRICE DI COVARIANZA DELLE sqr_dGamma
  
  double *cov_exp_sqr_dGamma=(double*)malloc(sizeof(double)*(Ncov)*(Ncov));
  std::string file_name_cov_matrix[1] = {"Input_corr_3pts/Matrice_covarianza_sqrt_dGamma_exp_D_to_pi.out"};

  double *cov_exp=(double*)calloc( (Ntot)*(Ntot), sizeof(double));

  ////// sqr_dGamma
  double *sqr_dGamma_BaBar=(double*)malloc(sizeof(double)*(NBaBar)),     *sigma_sqr_dGamma_BaBar=(double*)malloc(sizeof(double)*(NBaBar));
  double *sqr_dGamma_Cleo_D0=(double*)malloc(sizeof(double)*(NCleo_D0)), *sigma_sqr_dGamma_Cleo_D0=(double*)malloc(sizeof(double)*(NCleo_D0));
  double *sqr_dGamma_Cleo_Dp=(double*)malloc(sizeof(double)*(NCleo_Dp)), *sigma_sqr_dGamma_Cleo_Dp=(double*)malloc(sizeof(double)*(NCleo_Dp));
  double *sqr_dGamma_Bes3_D0=(double*)malloc(sizeof(double)*(NBes3_D0)), *sigma_sqr_dGamma_Bes3_D0=(double*)malloc(sizeof(double)*(NBes3_D0));
  double *sqr_dGamma_Bes3_Dp=(double*)malloc(sizeof(double)*(NBes3_Dp)), *sigma_sqr_dGamma_Bes3_Dp=(double*)malloc(sizeof(double)*(NBes3_Dp));

  ////// fpVcd
  double *fpVcd_BaBar=(double*)malloc(sizeof(double)*(NBaBar)),     *sigma_fpVcd_BaBar=(double*)malloc(sizeof(double)*(NBaBar));
  double *fpVcd_Cleo_D0=(double*)malloc(sizeof(double)*(NCleo_D0)), *sigma_fpVcd_Cleo_D0=(double*)malloc(sizeof(double)*(NCleo_D0));
  double *fpVcd_Cleo_Dp=(double*)malloc(sizeof(double)*(NCleo_Dp)), *sigma_fpVcd_Cleo_Dp=(double*)malloc(sizeof(double)*(NCleo_Dp));
  double *fpVcd_Belle=(double*)malloc(sizeof(double)*(NBelle)),     *sigma_fpVcd_Belle=(double*)malloc(sizeof(double)*(NBelle));

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

  ////// Vcd
  double *Vcd_BaBar=(double*)malloc(sizeof(double)*(NBaBar)*(Nev_tot+1));
  double *Vcd_Cleo_D0=(double*)malloc(sizeof(double)*(NCleo_D0)*(Nev_tot+1));
  double *Vcd_Cleo_Dp=(double*)malloc(sizeof(double)*(NCleo_Dp)*(Nev_tot+1));
  double *Vcd_Bes3_D0=(double*)malloc(sizeof(double)*(NBes3_D0)*(Nev_tot+1));
  double *Vcd_Bes3_Dp=(double*)malloc(sizeof(double)*(NBes3_Dp)*(Nev_tot+1));
  double *Vcd_Belle=(double*)malloc(sizeof(double)*(NBelle)*(Nev_tot+1));

  
	read_fpVcd( q2_BaBar, q2_Cleo_D0, q2_Cleo_Dp, q2_Belle, fpVcd_BaBar, fpVcd_Cleo_D0, fpVcd_Cleo_Dp, fpVcd_Belle, sigma_fpVcd_BaBar, sigma_fpVcd_Cleo_D0, sigma_fpVcd_Cleo_Dp, sigma_fpVcd_Belle, NBaBar, NCleo_D0, NCleo_Dp, NBelle);

	read_sqr_dGamma( q2_inf_BaBar, q2_inf_Cleo_D0, q2_inf_Cleo_Dp, q2_inf_Bes3_D0, q2_inf_Bes3_Dp, q2_sup_BaBar, q2_sup_Cleo_D0, q2_sup_Cleo_Dp, q2_sup_Bes3_D0, q2_sup_Bes3_Dp, q2_BaBar, q2_Cleo_D0, q2_Cleo_Dp, q2_Bes3_D0, q2_Bes3_Dp, sqr_dGamma_BaBar, sqr_dGamma_Cleo_D0, sqr_dGamma_Cleo_Dp, sqr_dGamma_Bes3_D0, sqr_dGamma_Bes3_Dp, sigma_sqr_dGamma_BaBar, sigma_sqr_dGamma_Cleo_D0, sigma_sqr_dGamma_Cleo_Dp, sigma_sqr_dGamma_Bes3_D0, sigma_sqr_dGamma_Bes3_Dp, NBaBar, NCleo_D0, NCleo_Dp, NBes3_D0, NBes3_Dp);

  // MATRICI DI COVARIANZA: "cov_exp_sqr_dGamma" SI RIFERISCE SOLO ALLE sqr_dGamma, MENTRE 
  //                        "cov_exp" CONTIENE ANCHE LE fpVckm NEL BLOCCO RELATIVO A BELLE
	read_matrix( Ncov, Ncov, cov_exp_sqr_dGamma, file_name_cov_matrix);

  for(int i = 0; i < Ncov; i++){
    for(int j = 0; j < Ncov; j++){
      cov_exp[i*Ntot+j] = cov_exp_sqr_dGamma[i*Ncov+j];
    }// j
  }// i

  for(int i = Ncov; i < Ntot; i++){
		cov_exp[i*Ntot+i] = pow( sigma_fpVcd_Belle[i-Ncov], 2);
	}// i

  //printf("\nH555\n\n");
  //check_read_data( q2_inf_BaBar, q2_inf_Cleo_D0, q2_inf_Cleo_Dp, q2_inf_Bes3_D0, q2_inf_Bes3_Dp, q2_sup_BaBar, q2_sup_Cleo_D0, q2_sup_Cleo_Dp, q2_sup_Bes3_D0, q2_sup_Bes3_Dp, q2_BaBar, q2_Cleo_D0, q2_Cleo_Dp, q2_Bes3_D0, q2_Bes3_Dp, q2_Belle, sqr_dGamma_BaBar, sqr_dGamma_Cleo_D0, sqr_dGamma_Cleo_Dp, sqr_dGamma_Bes3_D0, sqr_dGamma_Bes3_Dp, sigma_sqr_dGamma_BaBar, sigma_sqr_dGamma_Cleo_D0, sigma_sqr_dGamma_Cleo_Dp, sigma_sqr_dGamma_Bes3_D0, sigma_sqr_dGamma_Bes3_Dp, fpVcd_BaBar, fpVcd_Cleo_D0, fpVcd_Cleo_Dp, fpVcd_Belle, sigma_fpVcd_BaBar, sigma_fpVcd_Cleo_D0, sigma_fpVcd_Cleo_Dp, sigma_fpVcd_Belle, NBaBar, NCleo_D0, NCleo_Dp, NBes3_D0, NBes3_Dp, NBelle);

  //printf("\nCVCV\n\n");
  //print_matrix( Ncov, Ncov, cov_exp_sqr_dGamma, 1);

  free(fpVcd_BaBar);
  free(fpVcd_Cleo_D0);
  free(fpVcd_Cleo_Dp);
  free(sigma_fpVcd_BaBar);
  free(sigma_fpVcd_Cleo_D0);
  free(sigma_fpVcd_Cleo_Dp);
  ////////// FINE LEGGO I DATI SPERIMENTALI DI q2, fpVcd, sqr_dGamma E LA MATRICE DI COVARIANZA DELLE sqr_dGamma

  //////// CALCOLO Vcd AI VALORI DI q2 RELATIVI AGLI ESPERIMENTI

	Vckm_from_sqrdGamma_boot_new(   Vcd_BaBar, Nev_tot, Ntot,   q2_inf_BaBar,   q2_sup_BaBar, par_boot,   sqr_dGamma_BaBar, cov_exp,         0, Nsubgrp_1-1, 0, clusterfile, p_cube);
	Vckm_from_sqrdGamma_boot_new( Vcd_Cleo_D0, Nev_tot, Ntot, q2_inf_Cleo_D0, q2_sup_Cleo_D0, par_boot, sqr_dGamma_Cleo_D0, cov_exp, Nsubgrp_1, Nsubgrp_2-1, 0, clusterfile, p_cube);
	Vckm_from_sqrdGamma_boot_new( Vcd_Cleo_Dp, Nev_tot, Ntot, q2_inf_Cleo_Dp, q2_sup_Cleo_Dp, par_boot, sqr_dGamma_Cleo_Dp, cov_exp, Nsubgrp_2, Nsubgrp_3-1, 1, clusterfile, p_cube);
	Vckm_from_sqrdGamma_boot_new( Vcd_Bes3_D0, Nev_tot, Ntot, q2_inf_Bes3_D0, q2_sup_Bes3_D0, par_boot, sqr_dGamma_Bes3_D0, cov_exp, Nsubgrp_3, Nsubgrp_4-1, 0, clusterfile, p_cube);
	Vckm_from_sqrdGamma_boot_new( Vcd_Bes3_Dp, Nev_tot, Ntot, q2_inf_Bes3_Dp, q2_sup_Bes3_Dp, par_boot, sqr_dGamma_Bes3_Dp, cov_exp, Nsubgrp_4, Nsubgrp_5-1, 1, clusterfile, p_cube);
	Vckm_from_fpVckm_boot_new(      Vcd_Belle, Nev_tot, Ntot,       q2_Belle,                 par_boot,        fpVcd_Belle, cov_exp, Nsubgrp_5, Nsubgrp_6-1, fplus_func);

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

  free(cov_exp_sqr_dGamma);
  free(cov_exp);
  //////// FINE CALCOLO Vcd AI VALORI DI q2 RELATIVI AGLI ESPERIMENTI

  //////////// CREO IL VETTORE DI TUTTI I TERMINI
  double *q2_all=(double*)malloc(sizeof(double)*(Ntot)), *q2_all_sort=(double*)malloc(sizeof(double)*(Ntot));
  double *Vcd_all=(double*)malloc(sizeof(double)*(Ntot)*(Nev_tot+1)), *Vcd_all_sort=(double*)malloc(sizeof(double)*(Ntot));
  
  //q2_all
	append_six_arrays( q2_all, q2_BaBar, q2_Cleo_D0, q2_Cleo_Dp, q2_Bes3_D0, q2_Bes3_Dp, q2_Belle, NBaBar, NCleo_D0, NCleo_Dp, NBes3_D0, NBes3_Dp, NBelle, 0);

  //Vcd_all
	append_six_arrays( Vcd_all, Vcd_BaBar, Vcd_Cleo_D0, Vcd_Cleo_Dp, Vcd_Bes3_D0, Vcd_Bes3_Dp, Vcd_Belle, NBaBar, NCleo_D0, NCleo_Dp, NBes3_D0, NBes3_Dp, NBelle, Nev_tot);

  // ORDINO I VALORI IN q2
  for(int i = 0; i < Ntot; i++){
    q2_all_sort[i] = q2_all[i];
  }// i

// QUI VA CAMBIATA LA FLAG
#if defined(VCD_COV_LAT_TOT)
	////// SORT DI q2
	sort_vector( q2_all_sort, Ntot, q2_all);
#endif

  free(Vcd_BaBar);
  free(Vcd_Cleo_D0);
  free(Vcd_Cleo_Dp);
  free(Vcd_Bes3_D0);
  free(Vcd_Bes3_Dp);
  free(Vcd_Belle);
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
  free(fpVcd_Belle);
  free(sigma_fpVcd_Belle);
  free(q2_BaBar);
  free(q2_Cleo_D0);
  free(q2_Cleo_Dp);
  free(q2_Bes3_D0);
  free(q2_Bes3_Dp);
  free(q2_Belle);
  //////////// FINE CREO IL VETTORE DI TUTTI I TERMINI

  //// MATRICE DI COVARIANZA DEI Vcd
  int Size_blks[Nset_data] = { NBaBar, NCleo_D0, NCleo_Dp, NBes3_D0, NBes3_Dp, NBelle};
  
  double *Vcd_cov=(double*)malloc(sizeof(double)*(Ntot)*(Ntot)), *Vcd_cov_block=(double*)malloc(sizeof(double)*(Ntot)*(Ntot)), *var_Vckm=(double*)calloc( 1, sizeof(double));
  *var_Vckm = 0.0;
  
  covariance_syntetic_data( Vcd_cov, Ntot, Nev_tot, Nev_an, *var_Vckm, Vcd_all, analysis_in, analysis_fin, clusterfile, 0);
  free(var_Vckm);

  //// COSTRUISCO LA MATRICE Vcd_cov_block
  for(int i = 0; i < Ntot; i++){
    for(int j = 0; j < Ntot; j++){
      Vcd_cov_block[i*Ntot+j] = Vcd_cov[i*Ntot+j];
    }// j
  }// i

	block_matrix( Vcd_cov_block, Ntot, Nset_data, Size_blks);
  /////////////////////////////////////////////////////////

// QUI VA CAMBIATA LA FLAG
#if defined(VCD_COV_LAT_BLOCK)
	block_matrix( Vcd_cov, Ntot, Nset_data, Size_blks);
#endif  

// QUI VA CAMBIATA LA FLAG
#if defined(VCD_COV_LAT_TOT)
	//// SORT DELLA MATRICE DI COVARIANZA
	sort_matrix( Vcd_cov, Ntot, Ntot, q2_all, q2_all);

	//printf("\n\nWWWW Vcd_cov_sort \n\n");
	//print_matrix( Ntot, Ntot, Vcd_cov, 1);
#endif
  //// FINE MATRICE DI COVARIANZA DEI Vcd


  //////////////////////  MEDIA PESATA PER RICAVARE Vcd
  double *temp_Vcd_cov=(double*)malloc(sizeof(double)*(Ntot)*(Ntot)), *temp_Vcd_cov_block=(double*)malloc(sizeof(double)*(Ntot)*(Ntot));

  double *Vcd_fit_BaBar=(double*)malloc(sizeof(double)*(Nev_tot+1)),      *sigma_Vcd_fit_BaBar=(double*)malloc(sizeof(double));
  double *Vcd_fit_Cleo_D0=(double*)malloc(sizeof(double)*(Nev_tot+1)),    *sigma_Vcd_fit_Cleo_D0=(double*)malloc(sizeof(double));
  double *Vcd_fit_Cleo_Dp=(double*)malloc(sizeof(double)*(Nev_tot+1)),    *sigma_Vcd_fit_Cleo_Dp=(double*)malloc(sizeof(double));
  double *Vcd_fit_Bes3_D0=(double*)malloc(sizeof(double)*(Nev_tot+1)),    *sigma_Vcd_fit_Bes3_D0=(double*)malloc(sizeof(double));
  double *Vcd_fit_Bes3_Dp=(double*)malloc(sizeof(double)*(Nev_tot+1)),    *sigma_Vcd_fit_Bes3_Dp=(double*)malloc(sizeof(double));
  double *Vcd_fit_Belle=(double*)malloc(sizeof(double)*(Nev_tot+1)),      *sigma_Vcd_fit_Belle=(double*)malloc(sizeof(double));
  double *Vcd_fit_correlated=(double*)malloc(sizeof(double)*(Nev_tot+1)), *sigma_Vcd_fit_correlated=(double*)malloc(sizeof(double));

  const int Ncut = 15;
  double q2_cut_min = 0.2, q2_cut_max = q2_max, *q2_cut=(double*)malloc(sizeof(double)*(Ncut));
  double *Vcd_cut=(double*)malloc(sizeof(double)*(Ncut)), *sigma_Vcd_cut=(double*)malloc(sizeof(double)*(Ncut));

  for(int ic = 0; ic < Ncut; ic++){
    q2_cut[ic] = q2_cut_min + ((q2_cut_max - q2_cut_min)/( (double) (Ncut - 1)))*ic;
  }// ic
  
  for(int ic = 0; ic < Ncut; ic++){
  
  	Vcd_fit_BaBar[0]      = 0.0; 
  	Vcd_fit_Cleo_D0[0]    = 0.0;
  	Vcd_fit_Cleo_Dp[0]    = 0.0;
  	Vcd_fit_Bes3_D0[0]    = 0.0;
  	Vcd_fit_Bes3_Dp[0]    = 0.0;
  	Vcd_fit_Belle[0]      = 0.0;
  	Vcd_fit_correlated[0] = 0.0;
 
		for(int i = 0; i < Ntot; i++){
			for(int j = 0; j < Ntot; j++){
	      temp_Vcd_cov[i*Ntot+j] = Vcd_cov[i*Ntot+j];
				temp_Vcd_cov_block[i*Ntot+j] = Vcd_cov_block[i*Ntot+j];
			}// j
		}// i

		for(int i = 0; i < Ntot; i++){
			for(int j = 0; j < Ntot; j++){
				if( q2_all[i] > q2_cut[ic] || q2_all[j] > q2_cut[ic]){
					if(i == j){
						temp_Vcd_cov_block[i*Ntot+j] = 1.0;
					}else if(i != j){
						temp_Vcd_cov_block[i*Ntot+j] = 0.0;
					}
				}// if( q2_all[i] > q2_cut[ic] && q2_all[j] > q2_cut[ic])
			}// j
		}// i
    
		for(int i = 0; i < Ntot; i++){
			for(int j = 0; j < Ntot; j++){
				if( q2_all_sort[i] > q2_cut[ic] || q2_all_sort[j] > q2_cut[ic]){
					if(i == j){
						temp_Vcd_cov[i*Ntot+j] = 1.0;
					}else if(i != j){
						temp_Vcd_cov[i*Ntot+j] = 0.0;
					}
				}// if
			}// j
		}// i  
 
  	for(int iev = 1; iev <= Nev_tot; iev++){

    	/////// COSTRUISCO Vcd_all_sort
    	for(int i = 0; i < Ntot; i++){
      	Vcd_all_sort[i] = Vcd_all[i*(Nev_tot+1)+iev];
    	}// i

// QUI VA CAMBIATA LA FLAG
#if defined(VCD_COV_LAT_TOT)
			////// SORT DI Vcd
			sort_vector( Vcd_all_sort, Ntot, q2_all);
#endif  

			Vcd_fit_BaBar[iev]      = Vckm_const_fit_with_cut( Vcd_all, temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nev_tot,       0.0, Nsubgrp_1, iev, ic);
			Vcd_fit_Cleo_D0[iev]    = Vckm_const_fit_with_cut( Vcd_all, temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_1, Nsubgrp_2, iev, ic);
			Vcd_fit_Cleo_Dp[iev]    = Vckm_const_fit_with_cut( Vcd_all, temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_2, Nsubgrp_3, iev, ic);
			Vcd_fit_Bes3_D0[iev]    = Vckm_const_fit_with_cut( Vcd_all, temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_3, Nsubgrp_4, iev, ic);
			Vcd_fit_Bes3_Dp[iev]    = Vckm_const_fit_with_cut( Vcd_all, temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_4, Nsubgrp_5, iev, ic);
			Vcd_fit_Belle[iev]      = Vckm_const_fit_with_cut( Vcd_all, temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_5, Nsubgrp_6, iev, ic);

			Vcd_fit_correlated[iev] = Vckm_const_fit_with_cut_corr( Vcd_all_sort, temp_Vcd_cov, q2_all_sort, q2_cut, Ntot, 0, Nbin_max, ic);

    	Vcd_fit_BaBar[0]      +=  Vcd_fit_BaBar[iev]/Nev_tot; 
    	Vcd_fit_Cleo_D0[0]    +=  Vcd_fit_Cleo_D0[iev]/Nev_tot; 
    	Vcd_fit_Cleo_Dp[0]    +=  Vcd_fit_Cleo_Dp[iev]/Nev_tot; 
    	Vcd_fit_Bes3_D0[0]    +=  Vcd_fit_Bes3_D0[iev]/Nev_tot; 
    	Vcd_fit_Bes3_Dp[0]    +=  Vcd_fit_Bes3_Dp[iev]/Nev_tot; 
    	Vcd_fit_Belle[0]      +=  Vcd_fit_Belle[iev]/Nev_tot; 
    	Vcd_fit_correlated[0] +=  Vcd_fit_correlated[iev]/Nev_tot; 
 
  	}// iev

    if(ic == Ncut-1){
      *sigma_Vcd_fit_BaBar       =  sigma_bootstrap( Vcd_fit_BaBar,      analysis_in, analysis_fin, Nev_an, clusterfile);
      *sigma_Vcd_fit_Cleo_D0     =  sigma_bootstrap( Vcd_fit_Cleo_D0,    analysis_in, analysis_fin, Nev_an, clusterfile);
      *sigma_Vcd_fit_Cleo_Dp     =  sigma_bootstrap( Vcd_fit_Cleo_Dp,    analysis_in, analysis_fin, Nev_an, clusterfile);
      *sigma_Vcd_fit_Bes3_D0     =  sigma_bootstrap( Vcd_fit_Bes3_D0,    analysis_in, analysis_fin, Nev_an, clusterfile);
      *sigma_Vcd_fit_Bes3_Dp     =  sigma_bootstrap( Vcd_fit_Bes3_Dp,    analysis_in, analysis_fin, Nev_an, clusterfile);
      *sigma_Vcd_fit_Belle       =  sigma_bootstrap( Vcd_fit_Belle,      analysis_in, analysis_fin, Nev_an, clusterfile);
    }// if
    *sigma_Vcd_fit_correlated    =  sigma_bootstrap( Vcd_fit_correlated, analysis_in, analysis_fin, Nev_an, clusterfile);

    Vcd_cut[ic]     =  Vcd_fit_correlated[0];
    sigma_Vcd_cut[ic] =  *sigma_Vcd_fit_correlated;

  }// ic
    
  printf("RERE\n");
  printf("Vcd_fit_BaBar=%f sigma_Vcd_fit_BaBar=%f\n",           Vcd_fit_BaBar[0],      *sigma_Vcd_fit_BaBar);
  printf("Vcd_fit_Cleo_D0=%f sigma_Vcd_fit_Cleo_D0=%f\n",       Vcd_fit_Cleo_D0[0],    *sigma_Vcd_fit_Cleo_D0);
  printf("Vcd_fit_Cleo_Dp=%f sigma_Vcd_fit_Cleo_Dp=%f\n",       Vcd_fit_Cleo_Dp[0],    *sigma_Vcd_fit_Cleo_Dp);
  printf("Vcd_fit_Bes3_D0=%f sigma_Vcd_fit_Bes3_D0=%f\n",       Vcd_fit_Bes3_D0[0],    *sigma_Vcd_fit_Bes3_D0);
  printf("Vcd_fit_Bes3_Dp=%f sigma_Vcd_fit_Bes3_Dp=%f\n",       Vcd_fit_Bes3_Dp[0],    *sigma_Vcd_fit_Bes3_Dp);
  printf("Vcd_fit_Belle=%f sigma_Vcd_fit_Belle=%f\n",           Vcd_fit_Belle[0],      *sigma_Vcd_fit_Belle);
  printf("Vcd_fit_correlated=%f sigma_Vcd_fit_correlated=%f\n", Vcd_fit_correlated[0], *sigma_Vcd_fit_correlated);
 
  //////////////////////  FINE MEDIA PESATA PER RICAVARE Vcd
  //////////////////////  FINE CALCOLO Vcd

 
  //////////// SCRIVO IL FILE BOOTSTRAP PER Vcd
  char *file_out_Vcd=(char*)malloc(sizeof(char)*(LEN_NAME));
  FILE *fout_Vcd;

  if(output == 1){

    sprintf(file_out_Vcd, "OUTPUT_SMEAR/%s/%s/CKM/Vcd.out", dir_S[iS].c_str(), dir_E[iE].c_str());

    if ((fout_Vcd = fopen(file_out_Vcd, "w")) == NULL ){
      printf("Error opening the input file: file_out_Vcd\n");
      exit(EXIT_FAILURE);
    }
    
    for(int iev = 1; iev <= Nev_tot; iev++){
      fprintf(fout_Vcd, "%f\n", Vcd_fit_correlated[iev]);
    }
    fclose(fout_Vcd);
    free(file_out_Vcd);
  }// if(output == 1)
  //////////// FINE SCRIVO IL FILE BOOTSTRAP PER Vcd

  
  /////////// PLOT GRACE MEDIA PESATA Vcd
  double *temp_Vcd_all=(double*)malloc(sizeof(double)*(Ntot)), *temp_sigma_Vcd_all=(double*)malloc(sizeof(double)*(Ntot));

  for(int i = 0; i < Ntot; i++){
		temp_Vcd_all[i] = Vcd_all[i*(Nev_tot+1)+0];
		temp_sigma_Vcd_all[i] = sqrt(Vcd_cov_block[i+Ntot*i]);
	}// i

  printf("\n#HHHH\n");

	print_band_grace( q2_cut_min, q2_cut_max, Ncut, "#CUT_STEP", Vcd_cut, sigma_Vcd_cut, 0);

	print_points_grace( q2_cut, Vcd_cut, 0, Ncut-1, "Vcd_cut", 0);

	print_points_with_err_grace( q2_all, temp_Vcd_all, temp_sigma_Vcd_all,       0.0, Nsubgrp_1-1,   "Vcd_BaBar", 0);
	print_points_with_err_grace( q2_all, temp_Vcd_all, temp_sigma_Vcd_all, Nsubgrp_1, Nsubgrp_2-1, "Vcd_Cleo_D0", 1);
	print_points_with_err_grace( q2_all, temp_Vcd_all, temp_sigma_Vcd_all, Nsubgrp_2, Nsubgrp_3-1, "Vcd_Cleo_Dp", 1);
	print_points_with_err_grace( q2_all, temp_Vcd_all, temp_sigma_Vcd_all, Nsubgrp_3, Nsubgrp_4-1, "Vcd_Bes3_D0", 1);
	print_points_with_err_grace( q2_all, temp_Vcd_all, temp_sigma_Vcd_all, Nsubgrp_4, Nsubgrp_5-1, "Vcd_Bes3_Dp", 1);
	print_points_with_err_grace( q2_all, temp_Vcd_all, temp_sigma_Vcd_all, Nsubgrp_5, Nsubgrp_6-1,   "Vcd_Belle", 1);

	print_horizontal_mean_with_err_grace( q2_min, q2_max, Vcd_fit_correlated[0], *sigma_Vcd_fit_correlated, "Vcd Final Result", 0);

  free(temp_Vcd_all);
  free(temp_sigma_Vcd_all);
  /////////// PLOT GRACE MEDIA PESATA Vcd


  /////////// ESTRAZIONE Vckm DAL FIT SINGOLO
  double Vcd_single_BaBar,   sigma_Vcd_single_BaBar,   Xi2_BaBar,   Xi2_part_BaBar   = 0.0;
  double Vcd_single_Cleo_D0, sigma_Vcd_single_Cleo_D0, Xi2_Cleo_D0, Xi2_part_Cleo_D0 = 0.0;
  double Vcd_single_Cleo_Dp, sigma_Vcd_single_Cleo_Dp, Xi2_Cleo_Dp, Xi2_part_Cleo_Dp = 0.0;
  double Vcd_single_Bes3_D0, sigma_Vcd_single_Bes3_D0, Xi2_Bes3_D0, Xi2_part_Bes3_D0 = 0.0;
  double Vcd_single_Bes3_Dp, sigma_Vcd_single_Bes3_Dp, Xi2_Bes3_Dp, Xi2_part_Bes3_Dp = 0.0;
  double Vcd_single_Belle,   sigma_Vcd_single_Belle,   Xi2_Belle,   Xi2_part_Belle   = 0.0; 
  double Vcd_single_corr,    sigma_Vcd_single_corr,    Xi2_corr,    Xi2_tot          = 0.0;

 
  for(int ic = 0; ic < Ncut; ic++){

		/////// COSTRUISCO Vcd_all_sort
		for(int i = 0; i < Ntot; i++){
			Vcd_all_sort[i] = Vcd_all[i*(Nev_tot+1)+0];
		}// i

// QUI VA CAMBIATA LA FLAG
#if defined(VCD_COV_LAT_TOT)
		////// SORT DI Vcd
		sort_vector( Vcd_all_sort, Ntot, q2_all);
#endif  

		for(int i = 0; i < Ntot; i++){
			for(int j = 0; j < Ntot; j++){
	      temp_Vcd_cov[i*Ntot+j] = Vcd_cov[i*Ntot+j];
				temp_Vcd_cov_block[i*Ntot+j] = Vcd_cov_block[i*Ntot+j];
			}// j
		}// i

		for(int i = 0; i < Ntot; i++){
			for(int j = 0; j < Ntot; j++){
				if( q2_all[i] > q2_cut[ic] || q2_all[j] > q2_cut[ic]){
					if(i == j){
						temp_Vcd_cov_block[i*Ntot+j] = 1.0;
					}else if(i != j){
						temp_Vcd_cov_block[i*Ntot+j] = 0.0;
					}
				}// if( q2_all[i] > q2_cut[ic] && q2_all[j] > q2_cut[ic])
			}// j
		}// i
    
		for(int i = 0; i < Ntot; i++){
			for(int j = 0; j < Ntot; j++){
				if( q2_all_sort[i] > q2_cut[ic] || q2_all_sort[j] > q2_cut[ic]){
					if(i == j){
						temp_Vcd_cov[i*Ntot+j] = 1.0;
					}else if(i != j){
						temp_Vcd_cov[i*Ntot+j] = 0.0;
					}
				}// if
			}// j
		}// i  

		/////////// MEDIA PESATA
		Vcd_single_BaBar   = Vckm_single_fit_with_cut( Vcd_all, temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nev_tot,       0.0, Nsubgrp_1-1, ic);
		Vcd_single_Cleo_D0 = Vckm_single_fit_with_cut( Vcd_all, temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_1, Nsubgrp_2-1, ic);
		Vcd_single_Cleo_Dp = Vckm_single_fit_with_cut( Vcd_all, temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_2, Nsubgrp_3-1, ic);
		Vcd_single_Bes3_D0 = Vckm_single_fit_with_cut( Vcd_all, temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_3, Nsubgrp_4-1, ic);
		Vcd_single_Bes3_Dp = Vckm_single_fit_with_cut( Vcd_all, temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_4, Nsubgrp_5-1, ic);
		Vcd_single_Belle   = Vckm_single_fit_with_cut( Vcd_all, temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_5, Nsubgrp_6-1, ic);

		Vcd_single_corr    = Vckm_single_fit_with_cut_corr( Vcd_all_sort, temp_Vcd_cov, q2_all_sort, q2_cut, Ntot, 0.0,  Nbin_max-1, ic);


		sigma_Vcd_single_BaBar   = sigma_Vckm_single_fit_with_cut( temp_Vcd_cov_block, q2_all, q2_cut, Ntot,       0.0, Nsubgrp_1-1, ic);
		sigma_Vcd_single_Cleo_D0 = sigma_Vckm_single_fit_with_cut( temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nsubgrp_1, Nsubgrp_2-1, ic);
		sigma_Vcd_single_Cleo_Dp = sigma_Vckm_single_fit_with_cut( temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nsubgrp_2, Nsubgrp_3-1, ic);
		sigma_Vcd_single_Bes3_D0 = sigma_Vckm_single_fit_with_cut( temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nsubgrp_3, Nsubgrp_4-1, ic);
		sigma_Vcd_single_Bes3_Dp = sigma_Vckm_single_fit_with_cut( temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nsubgrp_4, Nsubgrp_5-1, ic);
		sigma_Vcd_single_Belle   = sigma_Vckm_single_fit_with_cut( temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nsubgrp_5, Nsubgrp_6-1, ic);

		sigma_Vcd_single_corr    = sigma_Vckm_single_fit_with_cut( temp_Vcd_cov, q2_all_sort, q2_cut, Ntot, 0.0, Nbin_max-1, ic);


		/////////// CALCOLO I CHI2
		Xi2_BaBar   = Xi2_with_cut( Vcd_single_BaBar,   Vcd_all, temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nev_tot,       0.0, Nsubgrp_1-1, ic);
		Xi2_Cleo_D0 = Xi2_with_cut( Vcd_single_Cleo_D0, Vcd_all, temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_1, Nsubgrp_2-1, ic);
		Xi2_Cleo_Dp = Xi2_with_cut( Vcd_single_Cleo_Dp, Vcd_all, temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_2, Nsubgrp_3-1, ic);
		Xi2_Bes3_D0 = Xi2_with_cut( Vcd_single_Bes3_D0, Vcd_all, temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_3, Nsubgrp_4-1, ic);
		Xi2_Bes3_Dp = Xi2_with_cut( Vcd_single_Bes3_Dp, Vcd_all, temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_4, Nsubgrp_5-1, ic);
		Xi2_Belle   = Xi2_with_cut( Vcd_single_Belle,   Vcd_all, temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_5, Nsubgrp_6-1, ic);

		Xi2_corr    = Xi2_with_cut_corr( Vcd_single_corr, Vcd_all_sort, temp_Vcd_cov, q2_all_sort, q2_cut, Ntot, 0.0, Nbin_max-1, ic);

		///// CONTRIBUTI PARZIALI DEGLI ESPERIMENTI AL Xi2 TOTALE
#if defined(VCD_COV_LAT_BLOCK) || defined(NO_VCD_COV_LAT) || defined(VCD_COV_LAT_DIAG)

		if(ic == Ncut-1){
			Xi2_part_BaBar   = (Xi2_with_cut( Vcd_single_corr, Vcd_all, temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nev_tot,       0.0, Nsubgrp_1-1, ic)*  (NBaBar-1.0))/(Nbin_max-1.0);
			Xi2_part_Cleo_D0 = (Xi2_with_cut( Vcd_single_corr, Vcd_all, temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_1, Nsubgrp_2-1, ic)*(NCleo_D0-1.0))/(Nbin_max-1.0);
			Xi2_part_Cleo_Dp = (Xi2_with_cut( Vcd_single_corr, Vcd_all, temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_2, Nsubgrp_3-1, ic)*(NCleo_Dp-1.0))/(Nbin_max-1.0);
			Xi2_part_Bes3_D0 = (Xi2_with_cut( Vcd_single_corr, Vcd_all, temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_3, Nsubgrp_4-1, ic)*(NBes3_D0-1.0))/(Nbin_max-1.0);
			Xi2_part_Bes3_Dp = (Xi2_with_cut( Vcd_single_corr, Vcd_all, temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_4, Nsubgrp_5-1, ic)*(NBes3_Dp-1.0))/(Nbin_max-1.0);

		#if defined(BELLE)
			Xi2_part_Belle   = (Xi2_with_cut( Vcd_single_corr, Vcd_all, temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_5, Nsubgrp_6-1, ic)*(NBelle-1.0))/(Nbin_max-1.0);
		#endif
		#if defined(NO_BELLE)
			Xi2_part_Belle   = 0.0;
		#endif

			Xi2_tot = Xi2_part_BaBar + Xi2_part_Cleo_D0 + Xi2_part_Cleo_Dp + Xi2_part_Bes3_D0 + Xi2_part_Bes3_Dp + Xi2_part_Belle;

		}// if(ic == Ncut-1)
#endif
  /////////// FINE CALCOLO I CHI2

		printf("\n\nVCKM ic=%d q2_cut=%.3f\n", ic, q2_cut[ic]);

		printf("Vcd_BaBar   = %f  %f  Xi2= %f\n", Vcd_single_BaBar,   sigma_Vcd_single_BaBar,   Xi2_BaBar);
		printf("Vcd_Cleo_D0 = %f  %f  Xi2= %f\n", Vcd_single_Cleo_D0, sigma_Vcd_single_Cleo_D0, Xi2_Cleo_D0);
		printf("Vcd_Cleo_Dp = %f  %f  Xi2= %f\n", Vcd_single_Cleo_Dp, sigma_Vcd_single_Cleo_Dp, Xi2_Cleo_Dp);    
		printf("Vcd_Bes3_D0 = %f  %f  Xi2= %f\n", Vcd_single_Bes3_D0, sigma_Vcd_single_Bes3_D0, Xi2_Bes3_D0);
		printf("Vcd_Bes3_Dp = %f  %f  Xi2= %f\n", Vcd_single_Bes3_Dp, sigma_Vcd_single_Bes3_Dp, Xi2_Bes3_Dp);
		printf("Vcd_Belle   = %f  %f  Xi2= %f\n", Vcd_single_Belle,   sigma_Vcd_single_Belle,   Xi2_Belle);  
		printf("Vcd_corr    = %f  %f  Xi2= %f\n", Vcd_single_corr,    sigma_Vcd_single_corr,    Xi2_corr);  
  
#if defined(VCD_COV_LAT_BLOCK) || defined(NO_VCD_COV_LAT) || defined(VCD_COV_LAT_DIAG)
		if(ic == Ncut-1){
			printf("Xi2_tot = %f\nXi2_frac_BaBar = %f\nXi2_frac_Cleo_D0 = %f\nXi2_frac_Cleo_Dp = %f\nXi2_frac_Bes3_D0 = %f\nXi2_frac_Bes3_Dp = %f\nXi2_frac_Belle = %f\n", Xi2_tot, Xi2_part_BaBar/Xi2_tot, Xi2_part_Cleo_D0/Xi2_tot, Xi2_part_Cleo_Dp/Xi2_tot, Xi2_part_Bes3_D0/Xi2_tot, Xi2_part_Bes3_Dp/Xi2_tot, Xi2_part_Belle/Xi2_tot);
		}
#endif

  }// ic
  /////////// FINE ESTRAZIONE Vcd DAL FIT SINGOLO


  //////////////////////  CALCOLO DI Vcd DA q2_max IN GIÙ
  q2_cut_min = 0.0;
  q2_cut_max = q2_max - 0.2;

  for(int ic = 0; ic < Ncut; ic++){
    q2_cut[ic] = q2_cut_min + ((q2_cut_max - q2_cut_min)/( (double) (Ncut - 1)))*ic;
  }// ic

  for(int ic = Ncut-1; ic >= 0; ic--){
  
  	Vcd_fit_BaBar[0]      = 0.0; 
  	Vcd_fit_Cleo_D0[0]    = 0.0;
  	Vcd_fit_Cleo_Dp[0]    = 0.0;
  	Vcd_fit_Bes3_D0[0]    = 0.0;
  	Vcd_fit_Bes3_Dp[0]    = 0.0;
  	Vcd_fit_Belle[0]      = 0.0;
  	Vcd_fit_correlated[0] = 0.0;

		Vcd_cut[ic]  = 0.0;

		for(int i = 0; i < Ntot; i++){
			for(int j = 0; j < Ntot; j++){
	      temp_Vcd_cov[i*Ntot+j] = Vcd_cov[i*Ntot+j];
				temp_Vcd_cov_block[i*Ntot+j] = Vcd_cov_block[i*Ntot+j];
			}// j
		}// i

		for(int i = 0; i < Ntot; i++){
			for(int j = 0; j < Ntot; j++){
				if( q2_all[i] <= q2_cut[ic] || q2_all[j] <= q2_cut[ic]){
					if(i == j){
						temp_Vcd_cov_block[i*Ntot+j] = 1.0;
					}else if(i != j){
						temp_Vcd_cov_block[i*Ntot+j] = 0.0;
					}
				}// if
			}// j
		}// i
    
		for(int i = 0; i < Ntot; i++){
			for(int j = 0; j < Ntot; j++){
				if( q2_all_sort[i] < q2_cut[ic] || q2_all_sort[j] < q2_cut[ic]){
					if(i == j){
						temp_Vcd_cov[i*Ntot+j] = 1.0;
					}else if(i != j){
						temp_Vcd_cov[i*Ntot+j] = 0.0;
					}
				}// if
			}// j
		}// i  
 
  	for(int iev = 1; iev <= Nev_tot; iev++){

    	/////// COSTRUISCO Vcd_all_sort
    	for(int i = 0; i < Ntot; i++){
      	Vcd_all_sort[i] = Vcd_all[i*(Nev_tot+1)+iev];
    	}// i

// QUI VA CAMBIATA LA FLAG
#if defined(VCD_COV_LAT_TOT)
			////// SORT DI Vcd
			sort_vector( Vcd_all_sort, Ntot, q2_all);
#endif  

			Vcd_fit_BaBar[iev]      = Vckm_const_fit_with_cut_reverse( Vcd_all, temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nev_tot,       0.0, Nsubgrp_1, iev, ic);
			Vcd_fit_Cleo_D0[iev]    = Vckm_const_fit_with_cut_reverse( Vcd_all, temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_1, Nsubgrp_2, iev, ic);
			Vcd_fit_Cleo_Dp[iev]    = Vckm_const_fit_with_cut_reverse( Vcd_all, temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_2, Nsubgrp_3, iev, ic);
			Vcd_fit_Bes3_D0[iev]    = Vckm_const_fit_with_cut_reverse( Vcd_all, temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_3, Nsubgrp_4, iev, ic);
			Vcd_fit_Bes3_Dp[iev]    = Vckm_const_fit_with_cut_reverse( Vcd_all, temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_4, Nsubgrp_5, iev, ic);
			Vcd_fit_Belle[iev]      = Vckm_const_fit_with_cut_reverse( Vcd_all, temp_Vcd_cov_block, q2_all, q2_cut, Ntot, Nev_tot, Nsubgrp_5, Nsubgrp_6, iev, ic);

			Vcd_fit_correlated[iev] = Vckm_const_fit_with_cut_corr_reverse( Vcd_all_sort, temp_Vcd_cov, q2_all_sort, q2_cut, Ntot, 0, Nbin_max, ic);

    	Vcd_fit_BaBar[0]      +=  Vcd_fit_BaBar[iev]/Nev_tot; 
    	Vcd_fit_Cleo_D0[0]    +=  Vcd_fit_Cleo_D0[iev]/Nev_tot; 
    	Vcd_fit_Cleo_Dp[0]    +=  Vcd_fit_Cleo_Dp[iev]/Nev_tot; 
    	Vcd_fit_Bes3_D0[0]    +=  Vcd_fit_Bes3_D0[iev]/Nev_tot; 
    	Vcd_fit_Bes3_Dp[0]    +=  Vcd_fit_Bes3_Dp[iev]/Nev_tot; 
    	Vcd_fit_Belle[0]      +=  Vcd_fit_Belle[iev]/Nev_tot; 
    	Vcd_fit_correlated[0] +=  Vcd_fit_correlated[iev]/Nev_tot; 
 
  	}// iev

    if(ic == 0){
			*sigma_Vcd_fit_BaBar       =  sigma_bootstrap( Vcd_fit_BaBar,      analysis_in, analysis_fin, Nev_an, clusterfile);
			*sigma_Vcd_fit_Cleo_D0     =  sigma_bootstrap( Vcd_fit_Cleo_D0,    analysis_in, analysis_fin, Nev_an, clusterfile);
			*sigma_Vcd_fit_Cleo_Dp     =  sigma_bootstrap( Vcd_fit_Cleo_Dp,    analysis_in, analysis_fin, Nev_an, clusterfile);
			*sigma_Vcd_fit_Bes3_D0     =  sigma_bootstrap( Vcd_fit_Bes3_D0,    analysis_in, analysis_fin, Nev_an, clusterfile);
			*sigma_Vcd_fit_Bes3_Dp     =  sigma_bootstrap( Vcd_fit_Bes3_Dp,    analysis_in, analysis_fin, Nev_an, clusterfile);
			*sigma_Vcd_fit_Belle       =  sigma_bootstrap( Vcd_fit_Belle,      analysis_in, analysis_fin, Nev_an, clusterfile);
    }// if
	  *sigma_Vcd_fit_correlated    =  sigma_bootstrap( Vcd_fit_correlated, analysis_in, analysis_fin, Nev_an, clusterfile);

		Vcd_cut[ic]     =  Vcd_fit_correlated[0];
		sigma_Vcd_cut[ic] =  *sigma_Vcd_fit_correlated;

  }// ic
    
  printf("RERE\n");
  printf("Vcd_fit_BaBar=%f sigma_Vcd_fit_BaBar=%f\n",           Vcd_fit_BaBar[0],      *sigma_Vcd_fit_BaBar);
  printf("Vcd_fit_Cleo_D0=%f sigma_Vcd_fit_Cleo_D0=%f\n",       Vcd_fit_Cleo_D0[0],    *sigma_Vcd_fit_Cleo_D0);
  printf("Vcd_fit_Cleo_Dp=%f sigma_Vcd_fit_Cleo_Dp=%f\n",       Vcd_fit_Cleo_Dp[0],    *sigma_Vcd_fit_Cleo_Dp);
  printf("Vcd_fit_Bes3_D0=%f sigma_Vcd_fit_Bes3_D0=%f\n",       Vcd_fit_Bes3_D0[0],    *sigma_Vcd_fit_Bes3_D0);
  printf("Vcd_fit_Bes3_Dp=%f sigma_Vcd_fit_Bes3_Dp=%f\n",       Vcd_fit_Bes3_Dp[0],    *sigma_Vcd_fit_Bes3_Dp);
  printf("Vcd_fit_Belle=%f sigma_Vcd_fit_Belle=%f\n",           Vcd_fit_Belle[0],      *sigma_Vcd_fit_Belle);
  printf("Vcd_fit_correlated=%f sigma_Vcd_fit_correlated=%f\n", Vcd_fit_correlated[0], *sigma_Vcd_fit_correlated);
 
  free(Vcd_fit_BaBar);
  free(Vcd_fit_Cleo_D0);
  free(Vcd_fit_Cleo_Dp);
  free(Vcd_fit_Bes3_D0);
  free(Vcd_fit_Bes3_Dp);
  free(Vcd_fit_Belle);
  free(Vcd_fit_correlated);
  free(sigma_Vcd_fit_BaBar);
  free(sigma_Vcd_fit_Cleo_D0);
  free(sigma_Vcd_fit_Cleo_Dp);
  free(sigma_Vcd_fit_Bes3_D0);
  free(sigma_Vcd_fit_Bes3_Dp);
  free(sigma_Vcd_fit_Belle);
  free(sigma_Vcd_fit_correlated);

  //////////////////////  FINE CALCOLO DI Vcd DA q2_max IN GIÙ

  /////////// PLOT GRACE MEDIA PESATA Vcd ORDINE INVERSO
  double *temp_Vcd_all_rev=(double*)malloc(sizeof(double)*(Ntot)), *temp_sigma_Vcd_all_rev=(double*)malloc(sizeof(double)*(Ntot));

  for(int i = 0; i < Ntot; i++){
		temp_Vcd_all_rev[i] = Vcd_all[i*(Nev_tot+1)+0];
		temp_sigma_Vcd_all_rev[i] = sqrt(Vcd_cov_block[i+Ntot*i]);
	}// i

  printf("\n#HRHR\n");

	print_band_grace( q2_cut_min, q2_cut_max, Ncut, "#CUT_STEP", Vcd_cut, sigma_Vcd_cut, 0);

	print_points_grace( q2_cut, Vcd_cut, 0, Ncut-1, "Vcd_cut", 0);

	print_points_with_err_grace( q2_all, temp_Vcd_all_rev, temp_sigma_Vcd_all_rev,       0.0, Nsubgrp_1-1,   "Vcd_BaBar", 0);
	print_points_with_err_grace( q2_all, temp_Vcd_all_rev, temp_sigma_Vcd_all_rev, Nsubgrp_1, Nsubgrp_2-1, "Vcd_Cleo_D0", 1);
	print_points_with_err_grace( q2_all, temp_Vcd_all_rev, temp_sigma_Vcd_all_rev, Nsubgrp_2, Nsubgrp_3-1, "Vcd_Cleo_Dp", 1);
	print_points_with_err_grace( q2_all, temp_Vcd_all_rev, temp_sigma_Vcd_all_rev, Nsubgrp_3, Nsubgrp_4-1, "Vcd_Bes3_D0", 1);
	print_points_with_err_grace( q2_all, temp_Vcd_all_rev, temp_sigma_Vcd_all_rev, Nsubgrp_4, Nsubgrp_5-1, "Vcd_Bes3_Dp", 1);
	print_points_with_err_grace( q2_all, temp_Vcd_all_rev, temp_sigma_Vcd_all_rev, Nsubgrp_5, Nsubgrp_6-1,   "Vcd_Belle", 1);

	print_horizontal_mean_with_err_grace( q2_min, q2_max, Vcd_fit_correlated[0], *sigma_Vcd_fit_correlated, "Vcd Final Result", 0);

  free(temp_Vcd_all_rev);
  free(temp_sigma_Vcd_all_rev);
  /////////// PLOT GRACE MEDIA PESATA Vcd ORDINE INVERSO

  free(q2_cut);
  free(q2_all);
  free(q2_all_sort);
  free(Vcd_all);
  free(Vcd_all_sort);
  free(Vcd_cov);
  free(Vcd_cov_block);
  free(temp_Vcd_cov);
  free(temp_Vcd_cov_block);
  free(Vcd_cut);
  free(sigma_Vcd_cut);


  ////////////////// CALCOLO RLU e-mu

  double  m_el = 0.5109989461e-3,  m_mu = 0.1056583745;
  double m_el2 = pow( m_el, 2),   m_mu2 = pow( m_mu,2);
  double q2_max_RLU = pow( MD_phys - Mpi_phys, 2);

  double Gamma_without_const_e, Gamma_without_const_mu;
  double *RLU_e_mu=(double*)calloc( Nev_tot+1, sizeof(double)), sigma_RLU_e_mu;

  //DICHIARAZIONI PER L'INTEGRALE
  int workspace_size=1000;
  gsl_integration_workspace *workspace=gsl_integration_workspace_alloc(workspace_size);

  gsl_function gsl_dGamma_without_const;
  gsl_dGamma_without_const.function = dGamma_without_const;

  int method=1;
  double abserr, epsabs = 0, epsrel = 1e-6;

  params_Gamma Gamma_params( 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0);

  Gamma_params.MH = MD_phys;
  Gamma_params.Ml = Mpi_phys;

  for(int iev = 1; iev <= Nev_tot; iev++){

    Gamma_params.res        = par_boot[iev][0];
    Gamma_params.A_fzero    = par_boot[iev][1];
    Gamma_params.A_fplus    = par_boot[iev][2];
    Gamma_params.polo_fzero = par_boot[iev][3];
    Gamma_params.polo_fplus = par_boot[iev][4];

    Gamma_params.mell       = m_el;
    gsl_dGamma_without_const.params = &Gamma_params;
    gsl_integration_qag( &gsl_dGamma_without_const, m_el2, q2_max_RLU, epsabs, epsrel, workspace_size, method, workspace, &Gamma_without_const_e, &abserr);

    Gamma_params.mell       = m_mu;
    gsl_dGamma_without_const.params = &Gamma_params;
    gsl_integration_qag( &gsl_dGamma_without_const, m_mu2, q2_max_RLU, epsabs, epsrel, workspace_size, method, workspace, &Gamma_without_const_mu, &abserr);

    RLU_e_mu[iev] = Gamma_without_const_mu/Gamma_without_const_e;
    RLU_e_mu[0]  += RLU_e_mu[iev]/Nev_tot;
  }// iev

  sigma_RLU_e_mu = sigma_bootstrap( RLU_e_mu, analysis_in, analysis_fin, Nev_an, clusterfile);

  printf("Ratio RLU_e_mu\n\n");

  printf("RLU = %.4f (%.2e)\n", RLU_e_mu[0], sigma_RLU_e_mu);

  gsl_integration_workspace_free(workspace);
  free(RLU_e_mu);

  ////////////////// FINE CALCOLO RLU e-mu

  return 0;
 
}// main


void chi2_cov( int &npar, double *deriv, double &f, double *par, int iflag){

  double vfit[Nmatrix-1];
 
	params_ffs parametri( par[0], par[1], par[2], par[3], par[4]);
 
  f = 0;
  
  for(int i = 0; i < Nmatrix-1; i++){
    if( i < Nq ){
      vfit[i] = f_fit[i] - fzero_func( q2_fit[i], &parametri);
    }else if( i >= Nq ){
      vfit[i] = f_fit[i] - fplus_func( q2_fit[i], &parametri);
    }
  }// i

  for(int i = 0; i < Nmatrix-1; i++){
    for(int j = 0; j < Nmatrix-1; j++){
      f += vfit[i]*inv_cov_fit[i][j]*vfit[j];
    }// j
  }// i

}// chi2_cov


double chi2_cov_num(double *parameters, int num_par, int block_par, int* dof){

  double vfit[Nmatrix-1], f = 0;
  int num_points = 0;
 
	params_ffs parametri( parameters[0], parameters[1], parameters[2], parameters[3], parameters[4]);
 
  for(int i = 0; i < Nmatrix-1; i++){
    if( i < Nq ){
      vfit[i] = f_fit[i] - fzero_func( q2_fit[i], &parametri);
    }else if( i >= Nq ){
      vfit[i] = f_fit[i] - fplus_func( q2_fit[i], &parametri);
    }
  }// i

  for(int i = 0; i < Nmatrix-1; i++){

    num_points = num_points + 1;
    
    for(int j = 0; j < Nmatrix-1; j++){

      f += vfit[i]*inv_cov_fit[i][j]*vfit[j];
      
    }// j
  }// i

  *dof = num_points - num_par - block_par;

  return f/(*dof);
  
}// chi2_cov_num


void chi2_boot( int &npar, double *deriv, double &f, double *par, int iflag){

  double vfit[Nmatrix-1];
  
	params_ffs parametri( par[0], par[1], par[2], par[3], par[4]);

  f = 0;
  
  for(int i = 0; i < Nmatrix-1; i++){

    if( i < Nq ){

      vfit[i] = f_fit[i] - fzero_func( q2_fit[i], &parametri);

    }else if( i >= Nq ){

      vfit[i] = f_fit[i] - fplus_func( q2_fit[i], &parametri);
      
    }

    
  }// i


  for(int i = 0; i < Nmatrix-1; i++){
    
      f += pow( vfit[i] , 2)/pow( sigma_f_fit[i], 2);

  }// i

  
}// chi2_boot


double chi2_boot_num(double *parameters, int num_par, int block_par, int* dof){

  double vfit[Nmatrix-1];
  int num_points = 0;
 
	params_ffs parametri( parameters[0], parameters[1], parameters[2], parameters[3], parameters[4]);
 
  double f = 0;
  
  for(int i = 0; i < Nmatrix-1; i++){

    if( i < Nq ){

      vfit[i] = f_fit[i] - fzero_func( q2_fit[i], &parametri);

    }else if( i >= Nq ){

      vfit[i] = f_fit[i] - fplus_func( q2_fit[i], &parametri);
      
    }

    
  }// i


  for(int i = 0; i < Nmatrix-1; i++){
    
    f += pow( vfit[i] , 2)/pow( sigma_f_fit[i], 2);

    num_points = num_points + 1;
      
  }// i
      

  *dof = num_points - num_par - block_par;

  return f/(*dof);

  
}// chi2_boot_num


double fzero_func( double q2, void *arg){

  params_ffs* params = (params_ffs*)arg;
  double res        = params -> res;
  double A_fzero    = params -> A_fzero;
  double polo_fzero = params -> polo_fzero;

  double z, z_term;
  double f0;
  
  z = ( sqrt( tp - q2 ) - sqrt( tp - t0 ) )/( sqrt( tp - q2 ) + sqrt( tp - t0 ) );
  z_term = ( z - z0 )*( 1 + ( z + z0 )/2  );
  
  f0 = ( res + A_fzero*z_term )/( 1 - polo_fzero*q2  );

  return f0;

}// fzero_func


double fplus_func( double q2, void *arg){

  params_ffs* params = (params_ffs*)arg;
  double res        = params -> res;
  double A_fplus    = params -> A_fplus;
  double polo_fplus = params -> polo_fplus;

  double z, z_term;
  double fp;
  
  z = ( sqrt( tp - q2 ) - sqrt( tp - t0 ) )/( sqrt( tp - q2 ) + sqrt( tp - t0 ) );
  z_term = ( z - z0 )*( 1 + ( z + z0 )/2  );
  
  fp = ( res + A_fplus*z_term )/( 1 - polo_fplus*q2  );

  return fp;

}// fplus_func


double p_cube( double q2, void *arg){

  params_Vckm* params   = (params_Vckm*)arg;
  double res        = params -> res;
  double A_fplus    = params -> A_fplus;
  double polo_fplus = params -> polo_fplus;
  double MH         = params -> MH;
  double Ml         = params -> Ml;

  double z, z_term;
  double fp, p3, result;

  z = ( sqrt( tp - q2 ) - sqrt( tp - t0 ) )/( sqrt( tp - q2 ) + sqrt( tp - t0 ) );
  z_term = ( z - z0 )*( 1 + ( z + z0 )/2  );
  
  fp = ( res + A_fplus*z_term )/( 1 - polo_fplus*q2  );
	p3 = pow( pow( ( pow( MH, 2) + pow( Ml, 2) - q2 )/(2*MH), 2) - pow( Ml, 2) , 1.5);

  result = p3*pow( fp, 2);

	return result;

}// p_cube


double dGamma_without_const( double q2, void *arg){

  params_Gamma* params   = (params_Gamma*)arg;
  double res        = params -> res;
  double A_fzero    = params -> A_fzero;
  double A_fplus    = params -> A_fplus;
  double polo_fzero = params -> polo_fzero;
  double polo_fplus = params -> polo_fplus;
  double MH         = params -> MH;
  double Ml         = params -> Ml;
  double mell       = params -> mell;

  double z, z_term;
  double mell2, fp2, f02, p, p3, result;

  mell2 = pow( mell, 2);

  z = ( sqrt( tp - q2 ) - sqrt( tp - t0 ) )/( sqrt( tp - q2 ) + sqrt( tp - t0 ) );
  z_term = ( z - z0 )*( 1 + ( z + z0 )/2  );
  
  fp2 = pow( ( res + A_fplus*z_term )/( 1 - polo_fplus*q2  ), 2);
  f02 = pow( ( res + A_fzero*z_term )/( 1 - polo_fzero*q2  ), 2);

  p  = pow( pow( ( pow( MH, 2) + pow( Ml, 2) - q2 )/(2*MH), 2) - pow( Ml, 2) , 0.5);
  p3 = pow( p, 3);

  result = pow( (q2 - mell2)/(q2), 2)*( (1.0 + mell2/(2*q2))*p3*fp2 + ((3*mell2)/(8*q2))*pow( (pow(MH,2) - pow(Ml,2))/MH, 2)*p*f02 );

  return result;

}// dGamma_without_const


void read_fpVcd( double* q2_BaBar, double* q2_Cleo_D0, double* q2_Cleo_Dp, double* q2_Belle, double* fpVcd_BaBar, double* fpVcd_Cleo_D0, double* fpVcd_Cleo_Dp, double* fpVcd_Belle, double* sigma_fpVcd_BaBar, double* sigma_fpVcd_Cleo_D0, double* sigma_fpVcd_Cleo_Dp, double* sigma_fpVcd_Belle, int NBaBar, int NCleo_D0, int NCleo_Dp, int NBelle){

  FILE *fr_fpVcd;
  int t;

  if ((fr_fpVcd = fopen( "Input_corr_3pts/fpVcd_sperimentale_Dpi.out", "r")) == NULL ){
    printf("Error opening the file to read: open_fpVcd\n");
    exit(EXIT_FAILURE);
  }

  for(int i = 0; i < NBaBar; i++){
    t = fscanf(fr_fpVcd,"%lf %lf %lf\n", &q2_BaBar[i], &fpVcd_BaBar[i], &sigma_fpVcd_BaBar[i]);
  }
  for(int i = 0; i < NCleo_D0; i++){
    t = fscanf(fr_fpVcd,"%lf %lf %lf\n", &q2_Cleo_D0[i], &fpVcd_Cleo_D0[i], &sigma_fpVcd_Cleo_D0[i]); 
  }
  for(int i = 0; i < NCleo_Dp; i++){
    t = fscanf(fr_fpVcd,"%lf %lf %lf\n", &q2_Cleo_Dp[i], &fpVcd_Cleo_Dp[i], &sigma_fpVcd_Cleo_Dp[i]);
  }
  for(int i = 0; i < NBelle; i++){
    t = fscanf(fr_fpVcd,"%lf %lf %lf\n", &q2_Belle[i], &fpVcd_Belle[i], &sigma_fpVcd_Belle[i]);
  }
  
  fclose(fr_fpVcd);
}// read_fpVcd



void read_sqr_dGamma( double* q2_inf_BaBar, double* q2_inf_Cleo_D0, double* q2_inf_Cleo_Dp, double* q2_inf_Bes3_D0, double* q2_inf_Bes3_Dp, double* q2_sup_BaBar, double* q2_sup_Cleo_D0, double* q2_sup_Cleo_Dp, double* q2_sup_Bes3_D0, double* q2_sup_Bes3_Dp, double* q2_BaBar, double* q2_Cleo_D0, double* q2_Cleo_Dp, double* q2_Bes3_D0, double* q2_Bes3_Dp, double* sqr_dGamma_BaBar, double* sqr_dGamma_Cleo_D0, double* sqr_dGamma_Cleo_Dp, double* sqr_dGamma_Bes3_D0, double* sqr_dGamma_Bes3_Dp, double* sigma_sqr_dGamma_BaBar, double* sigma_sqr_dGamma_Cleo_D0, double* sigma_sqr_dGamma_Cleo_Dp, double* sigma_sqr_dGamma_Bes3_D0, double* sigma_sqr_dGamma_Bes3_Dp, int NBaBar, int NCleo_D0, int NCleo_Dp, int NBes3_D0, int NBes3_Dp){

  FILE *fr_sqr_dGamma;
  int t;
  
  if ((fr_sqr_dGamma = fopen("Input_corr_3pts/sqrt_delta_gamma_sperimentali_D_to_pi.out", "r")) == NULL ){
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


void check_read_data( double* q2_inf_BaBar, double* q2_inf_Cleo_D0, double* q2_inf_Cleo_Dp, double* q2_inf_Bes3_D0, double* q2_inf_Bes3_Dp, double* q2_sup_BaBar, double* q2_sup_Cleo_D0, double* q2_sup_Cleo_Dp, double* q2_sup_Bes3_D0, double* q2_sup_Bes3_Dp, double* q2_BaBar, double* q2_Cleo_D0, double* q2_Cleo_Dp, double* q2_Bes3_D0, double* q2_Bes3_Dp, double* q2_Belle, double* sqr_dGamma_BaBar, double* sqr_dGamma_Cleo_D0, double* sqr_dGamma_Cleo_Dp, double* sqr_dGamma_Bes3_D0, double* sqr_dGamma_Bes3_Dp, double* sigma_sqr_dGamma_BaBar, double* sigma_sqr_dGamma_Cleo_D0, double* sigma_sqr_dGamma_Cleo_Dp, double* sigma_sqr_dGamma_Bes3_D0, double* sigma_sqr_dGamma_Bes3_Dp, double* fpVcd_BaBar, double* fpVcd_Cleo_D0, double* fpVcd_Cleo_Dp, double* fpVcd_Belle, double* sigma_fpVcd_BaBar, double* sigma_fpVcd_Cleo_D0, double* sigma_fpVcd_Cleo_Dp, double* sigma_fpVcd_Belle, int NBaBar, int NCleo_D0, int NCleo_Dp, int NBes3_D0, int NBes3_Dp, int NBelle){

  //BABAR
  for(int i = 0; i < NBaBar; i++){
    printf("%f %f %f %f %f %.12f %.12f\n", q2_BaBar[i], q2_inf_BaBar[i], q2_sup_BaBar[i], fpVcd_BaBar[i], sigma_fpVcd_BaBar[i], sqr_dGamma_BaBar[i], sigma_sqr_dGamma_BaBar[i]);
  }
  printf("&\n");

  //CLEO-D0  
  for(int i = 0; i < NCleo_D0; i++){
    printf("%f %f %f %f %f %.12f %.12f\n", q2_Cleo_D0[i], q2_inf_Cleo_D0[i], q2_sup_Cleo_D0[i], fpVcd_Cleo_D0[i], sigma_fpVcd_Cleo_D0[i], sqr_dGamma_Cleo_D0[i], sigma_sqr_dGamma_Cleo_D0[i]);
  }
  printf("&\n");
  
  //CLEO-Dp
  for(int i = 0; i < NCleo_Dp; i++){
    printf("%f %f %f %f %f %.12f %.12f\n", q2_Cleo_Dp[i], q2_inf_Cleo_Dp[i], q2_sup_Cleo_Dp[i], fpVcd_Cleo_Dp[i], sigma_fpVcd_Cleo_Dp[i], sqr_dGamma_Cleo_Dp[i], sigma_sqr_dGamma_Cleo_Dp[i]);
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
    printf("%f %f %f\n", q2_Belle[i], fpVcd_Belle[i], sigma_fpVcd_Belle[i]);
  }

}// check_read_data


void Vckm_from_sqrdGamma_boot_new( double* Vckm, int Nev, int Ncov, double* q2_inf, double* q2_sup, double par_boot[Nev_tot+1][npar], double* sqr_dGamma, double* cov_exp_sqr_dGamma, int i_start, int i_end, int isospin, int clusterfile, double (*p_cube)( double, void*)){

  int Nexp = i_end-i_start+1; 

  //DICHIARO gsl_vectors E gsl_matrix per la generazione random delle sqrdGamma
  gsl_vector* sqr_dGamma_mu = gsl_vector_alloc(Nexp);
  gsl_vector* sqr_dGamma_random = gsl_vector_alloc(Nexp);
  gsl_matrix* gsl_cov_sqr_dGamma = gsl_matrix_alloc( Nexp, Nexp);

  for(int i = 0; i < Nexp; i++){
    gsl_vector_set( sqr_dGamma_mu, i, sqr_dGamma[i]);
  }// i

  for(int i = i_start; i <= i_end; i++){
    for(int j = i_start; j <= i_end; j++){
      gsl_matrix_set( gsl_cov_sqr_dGamma, i-i_start, j-i_start, cov_exp_sqr_dGamma[i*Ncov+j]/(clusterfile - 2.0));
    }// i
  }// j
	gsl_linalg_cholesky_decomp1( gsl_cov_sqr_dGamma);

  //DICHIARO IL GENERATORE RANDOM
	const gsl_rng_type* ran_T;
	gsl_rng* ran_gen;

	gsl_rng_env_setup();

	ran_T = gsl_rng_default;
	ran_gen = gsl_rng_alloc(ran_T);

	gsl_rng_set( ran_gen, time(NULL));

  //DICHIARAZIONI PER L'INTEGRALE
  int workspace_size=1000;
	gsl_integration_workspace *workspace=gsl_integration_workspace_alloc(workspace_size);

  gsl_function integrand;
	integrand.function = p_cube;

  int method=1;
  double int_result, abserr, epsabs = 0, epsrel = 1e-6, MH, Ml, C;

	params_Vckm parameters( 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0);

  if(isospin == 0){
    MH = MD0_pdg;
    Ml = Mpim_pdg;
    C  = sqrt(24.0*pow(PI,3))/GF;
  }else if(isospin == 1){
    MH = MDp_pdg;
    Ml = Mpi0_pdg;
    C  = sqrt(48.0*pow(PI,3))/GF;
  }
	parameters.MH = MH;
	parameters.Ml = Ml;

  for(int i = 0; i < Nexp; i++){
    Vckm[i*(Nev+1)+0] = 0.0;
  }// i

	for(int iev = 1; iev <= Nev; iev++){

		gsl_ran_multivariate_gaussian( ran_gen, sqr_dGamma_mu, gsl_cov_sqr_dGamma, sqr_dGamma_random);

		parameters.res        = par_boot[iev][0];
		parameters.A_fzero    = par_boot[iev][1];
		parameters.A_fplus    = par_boot[iev][2];
		parameters.polo_fzero = par_boot[iev][3];
		parameters.polo_fplus = par_boot[iev][4];

		integrand.params = &parameters;

		for(int i = 0; i < Nexp; i++){

      gsl_integration_qag( &integrand, q2_inf[i], q2_sup[i], epsabs, epsrel, workspace_size, method, workspace, &int_result, &abserr);

      Vckm[i*(Nev+1)+iev] = (C*gsl_vector_get( sqr_dGamma_random, i))/sqrt(int_result);
      Vckm[i*(Nev+1)+0] += Vckm[i*(Nev+1)+iev]/Nev;

    }// i
  }// iev

	gsl_integration_workspace_free(workspace);
  gsl_vector_free(sqr_dGamma_mu);
  gsl_vector_free(sqr_dGamma_random);
  gsl_matrix_free(gsl_cov_sqr_dGamma);
	gsl_rng_free(ran_gen);

}// Vckm_from_sqrdGamma_boot_new


void Vckm_from_fpVckm_boot_new( double* Vckm, int Nev, int Ncov, double* q2, double par_boot[Nev_tot+1][npar], double* fpVckm, double* cov_exp_fpVckm, int i_start, int i_end, double (*fplus_func)( double, void*)){

  int Nexp = i_end-i_start+1; 

  //DICHIARO gsl_vectors E gsl_matrix per la generazione random dei fpVckm
  gsl_vector*      fpVckm_mu = gsl_vector_alloc(Nexp);
  gsl_vector*  fpVckm_random = gsl_vector_alloc(Nexp);
  gsl_matrix* gsl_cov_fpVckm = gsl_matrix_alloc( Nexp, Nexp);

  for(int i = 0; i < Nexp; i++){
    gsl_vector_set( fpVckm_mu, i, fpVckm[i]);
  }// i

  for(int i = i_start; i <= i_end; i++){
    for(int j = i_start; j <= i_end; j++){
      gsl_matrix_set( gsl_cov_fpVckm, i-i_start, j-i_start, cov_exp_fpVckm[i*Ncov+j]/(clusterfile - 2.0));
    }// i
  }// j
	gsl_linalg_cholesky_decomp1( gsl_cov_fpVckm);

  //DICHIARO IL GENERATORE RANDOM
	const gsl_rng_type* ran_T;
	gsl_rng* ran_gen;

	gsl_rng_env_setup();

	ran_T = gsl_rng_default;
	ran_gen = gsl_rng_alloc(ran_T);

	gsl_rng_set( ran_gen, time(NULL));

  // CALCOLO Vckm
  double **fp_synt=(double**)malloc(sizeof(double)*(Nexp));
  for(int i = 0; i < Nexp; i++) fp_synt[i]=(double*)malloc(sizeof(double)*(Nev+1));

	params_ffs parameters( 1.0, 2.0, 3.0, 4.0, 5.0);

  for(int i = 0; i < Nexp; i++){
    fp_synt[i][0] = 0.0;
    Vckm[i*(Nev+1)+0] = 0.0;
  }// i


	for(int iev = 1; iev <= Nev; iev++){

		gsl_ran_multivariate_gaussian( ran_gen, fpVckm_mu, gsl_cov_fpVckm, fpVckm_random);

		parameters.res        = par_boot[iev][0];
		parameters.A_fzero    = par_boot[iev][1];
		parameters.A_fplus    = par_boot[iev][2];
		parameters.polo_fzero = par_boot[iev][3];
		parameters.polo_fplus = par_boot[iev][4];

		for(int i = 0; i < Nexp; i++){

      fp_synt[i][iev] = fplus_func( q2[i], &parameters);      
      fp_synt[i][0] += fp_synt[i][iev]/Nev;

      Vckm[i*(Nev+1)+iev] = gsl_vector_get( fpVckm_random, i)/fp_synt[i][iev];
      Vckm[i*(Nev+1)+0] += Vckm[i*(Nev+1)+iev]/Nev;
    }// i
  }// iev

  free(fp_synt);
  gsl_vector_free(fpVckm_mu);
  gsl_vector_free(fpVckm_random);
  gsl_matrix_free(gsl_cov_fpVckm);
	gsl_rng_free(ran_gen);

}// Vckm_from_fpVckm_boot_new


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


double Vckm_const_fit_with_cut_reverse( double* Vckm, double* Vckm_cov, double* q2_all, double* q2_cut, int Ntot, int Nev_tot, int i_start, int i_end, int iev, int ic){

  double num = 0.0, den = 0.0, *inv_Vckm_cov=(double*)malloc(sizeof(double)*(Ntot)*(Ntot));

	LU_invert( inv_Vckm_cov, Vckm_cov, Ntot);    

	for(int i = i_start; i < i_end; i++){
		for(int j = i_start; j < i_end; j++){

			if( q2_all[i] > q2_cut[ic] && q2_all[j] > q2_cut[ic]){
	  
				num += inv_Vckm_cov[i*Ntot+j]*Vckm[j*(Nev_tot+1)+iev];
				den += inv_Vckm_cov[i*Ntot+j];
	  
			}// if
		}// j
	}// i

  free(inv_Vckm_cov);

	return num/den;

}// Vckm_const_fit_with_cut_reverse

double Vckm_const_fit_with_cut_corr_reverse( double* Vckm, double* Vckm_cov, double* q2_all, double* q2_cut, int Ntot, int i_start, int i_end, int ic){

  double num = 0.0, den = 0.0, *inv_Vckm_cov=(double*)malloc(sizeof(double)*(Ntot)*(Ntot));

	LU_invert( inv_Vckm_cov, Vckm_cov, Ntot);    

	for(int i = i_start; i < i_end; i++){
		for(int j = i_start; j < i_end; j++){

			if( q2_all[i] > q2_cut[ic] && q2_all[j] > q2_cut[ic]){
	  
				num += inv_Vckm_cov[i*Ntot+j]*Vckm[j];
				den += inv_Vckm_cov[i*Ntot+j];
	  
			}// if
		}// j
	}// i

  free(inv_Vckm_cov);

	return num/den;

}// Vckm_const_fit_with_cut_corr_reverse






























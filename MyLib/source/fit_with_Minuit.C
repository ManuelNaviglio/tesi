#include <TMinuit.h>
#include <TF1.h>            // questo serve per la funzione cout

void fit_with_Minuit( double* mu_out, double* sigma, double* cov_out, double* step, double* par, double* min, double* max, std::string* cpar, double* fix_par, double epsilon, int npar, int nfix, int block_par, void (*chi2)( int&, double*, double&, double*, int), double (*chi2_num)(double*, int, int, int*), std::string check_string, int iev){

	double *chi2_f=(double*)calloc( 1, sizeof(double));
  int num_par = 0, dof;
  //int *pdof = ; 
 
	TMinuit minuit(npar);
  minuit.SetFCN(chi2);
  minuit.SetErrorDef(1.);
  
  for(int ipar = 0; ipar < npar; ipar++){
    minuit.DefineParameter( ipar, cpar[ipar].c_str(), par[ipar], step[ipar], min[ipar], max[ipar]);
  }// ipar
  
  for(int ipar = 0; ipar < npar; ipar++){
    for(int j = 0; j < nfix; j++){
      
      if( ipar == fix_par[j] ){
				minuit.DefineParameter(ipar, cpar[ipar].c_str(), 0, step[j], min[j], max[j]);
				minuit.FixParameter(ipar);	
      }// if( ipar == fix_par[j] )
      
    }// j
  }// ipar
  
  minuit.Migrad();
 
	printf("\n"); 
  for(int ipar = 0; ipar < npar; ipar++){
    minuit.GetParameter( ipar, mu_out[ipar], sigma[ipar]);
    printf("GG%dG\t%s=%f\tsigma=%f\n", iev, cpar[ipar].c_str(), mu_out[ipar], sigma[ipar]);
  }// ipar
  printf("\n");

  ////////////////// PRENDO LA MATRICE DI COVARIANZA DEI PARAMETRI
  minuit.mnemat( cov_out, npar);
  ////////////////// FINE PRENDO LA MATRICE DI COVARIANZA DEI PARAMETRI

  for(int ipar = 0; ipar < npar; ipar++){
    if( fabs(mu_out[ipar]) > epsilon){
      num_par += 1;
    }
  }// ipar
  
  *chi2_f = chi2_num(mu_out, num_par, block_par, &dof);

  //printf("%s chi2_ev[%d] = %f\t dof = %d\n", check_string.c_str(), iev, *chi2_f, *dof);
  printf("%s chi2_ev[%d] = %f\t dof = %d\n", check_string.c_str(), iev, *chi2_f, dof);
  
  //num_points = 0;
  //num_par = 0;
  //dof = 0;

	free(chi2_f);

}// fit_with_Minuit


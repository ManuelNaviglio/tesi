#include <TMinuit.h>
#include <TF1.h>
#include <fstream>
#include "fit_2pts_func.h"
#include "clean_str_cl.h"
#include "stat_analysis_func.h"
#include "const_fit.h"
#include "definitions.h"

using namespace std;

const int LEN_NAME = 1024;
#define PI 3.141592653589793

const int npar = 2;

int ibeta, imusea, im1, im2, ith1, ith2, ismear;
int t_min, t_max, t_dat;

const int n_dat = 50;

const int num_rows_half = 50;    // ATTENTO CHE QUESTO VALORE È FISSATO (CERTAMENTE DA MODIFICARE)

////// VARIABILI GLOBALI PER IL FIT
double x_fit[n_dat];
double y_fit[n_dat];
double ey_fit[n_dat];
double meff_check_glb, err_meff_check_glb;
double Z_check_glb, err_Z_check_glb;
///////////////////////////////////

int T_max[5] = {64, 48, 64, 48, 96}, Tmax;

int num_points = 0;

void Jack_Knife_Average( int, int, double *, double *);
void sigma_Jack_Knife_Average( int, int, double *, double *);
void effective_mass( int, int, double, double *, double *);
void m_eff_and_Z_check(int, int, double *, double *);
void Z_function_JK( double *, double *, double *, double *, int);
void corr_func_fit( int, int, int, int, double *, double *, double *, double *, int, int, int, int, int);
void output_Energy_and_Z_from_2pts(int, double *, double *, const string *, const string [Nbeta][Nmusea], const string *, const string *, const string *, const string *, int, int, int, int, int, int);
void chi2( int &, double *, double &, double *, int);
double chi_num( double *);
double corr_function( double, double *, int);
void out_eff_mass_graph_grace(int, int, int, int, double *, double *, double *, const string *, const string [Nbeta][Nmusea], const string *, const string *, const string *, const string *, int, int, int, int, int, int);


double sinh_DR( double, int);
double stdDR( double, int);
double my_arcsinh( double);

int main(){

  //////////////////// LEGGO IL FILE DI INPUT & TIME INTERVALS
  
  FILE *fi;

  string input_file[1] = {"Input_corr_2pts/time_intervals.out"};
  
  if ((fi = fopen("Input_corr_2pts/file_input_corr_2pts.out", "r")) == NULL ){
    printf("Error opening the input file!!\n");
    exit(EXIT_FAILURE);
  }
  fscanf(fi, "%d %d %d %d %d %d %d", &ibeta, &imusea, &im1, &im2, &ith1, &ith2, &ismear);
  fclose(fi);
  
  read_time_intervals( input_file, t_min, t_max, t_dat, im1, ibeta, im2);

  printf("\n\nTIME_INTERVAL: tmin=%d tmax=%d tdat=%d\n\n", t_min, t_max, t_dat);
  
  Tmax = T_max[ibeta];
  
  //////////////////// FINE LEGGO IL FILE DI INPUT & TIME INTERVALS
  
  ///// LEGGO IL CORRELATORE A 2 PUNTI

  double *corr_2pts=(double*)malloc(sizeof(double)*(clusterfile)*(Tmax));
  
  read_corr_2pts_from_2pts( corr_2pts, Tmax, clusterfile, beta_V, mu_sea_1, n_conf, mu_sea_2, m_, th_, sme_, ibeta, imusea, im1, im2, ith1, ismear);
  
  ///// FINE LEGGO IL CORRELATORE A 2 PUNTI
  
  /////////////////////////////////////////
  //                                     // 
  // START DOING THE JACK-KNIFE AVERAGE. //
  //                                     //
  /////////////////////////////////////////

  double *sigma_JK=(double*)malloc(sizeof(double)*(Tmax)), *corr_JK=(double*)malloc(sizeof(double)*(clusterfile)*(Tmax));
  
  Jack_Knife_Average(clusterfile, Tmax, corr_2pts, corr_JK); // funzione da cambiare

  sigma_Jack_Knife_Average(clusterfile, Tmax, corr_2pts, sigma_JK);
  
  free(corr_2pts);
  
  /////////////////////////////////////////////////////////////
  //                                                         //
  // START ANALYSIS OF THE EFFECTIVE MASSES FOR EACH CLUSTER //
  //                                                         //
  /////////////////////////////////////////////////////////////

  double prec = 0.000001;

  double *m_eff_array=(double*)malloc(sizeof(double)*(clusterfile)*(Tmax/2));
  
  effective_mass(clusterfile, Tmax, prec, corr_JK, m_eff_array);

  m_eff_and_Z_check(clusterfile, Tmax, corr_JK, m_eff_array); // questo forse si può spostare
  /*
  for(int ijk = 0; ijk < clusterfile; ijk++){
    printf("Cluster %d\n", ijk);

    for(int t = 0; t < Tmax/2; t++){    
      
      printf("%d %+016.16g\n", t, m_eff_array[ijk*(Tmax/2)+t]);

    }// t
  }// ijk
  */
  double mass_out[clusterfile], z_term[clusterfile], sigma_mass_out[1];
  
  corr_func_fit(clusterfile, Tmax, t_min, t_dat, corr_JK, sigma_JK, mass_out, z_term, ibeta, imusea, im1, im2, ith1);
  
  output_Energy_and_Z_from_2pts(clusterfile, mass_out, z_term, beta_V, mu_sea_1, mu_sea_2, m_, th_, sme_, ibeta, imusea, im1, im2, ith1, ismear);
  
  sigma_Jack_Knife_Average( clusterfile, 1, mass_out, sigma_mass_out);

  out_eff_mass_graph_grace(clusterfile, Tmax, t_min, t_max, m_eff_array, mass_out, sigma_mass_out, beta_V, mu_sea_1, mu_sea_2, m_, th_, sme_, ibeta, imusea, im1, im2, ith1, ismear);
  
  free(corr_JK);
  free(sigma_JK);
  free(m_eff_array);
  
  return 0;
  
}








void Jack_Knife_Average(int clusterfile, int Nt, double *array, double *corr_JK){
  
  int k;
  double sum = 0;
  
  double *average_JK=(double*)malloc(sizeof(double)*(Nt));
  
  //Next block writes the JK-average for each row in average_JK[]
  for (int t = 0; t < Nt; t++){
    for (int ijk = 0; ijk < clusterfile-1; ijk++){
      
      sum += array[ijk*Nt + t];

    }// ijk
    
    average_JK[t] = sum/(clusterfile-1);
    
    sum =  0;
  }// t

  ////////// QUESTA È LA PARTE CHE VA MODIFICATA
  k = 0;
  
  for(int j = 0; j < (clusterfile*Nt); j++){
    
    if( (j >= (clusterfile-1)*Nt) && (j < clusterfile*Nt) ){
      
      corr_JK[k] = average_JK[j - (clusterfile-1)*Nt];
      
    } else {
      
      for(int ijk = 0; ijk < clusterfile-1; ijk++){
	
	if( (j >= ijk*Nt) && (j < (ijk+1)*Nt) ){
	  
	  corr_JK[k] = array[j];
	  
	}// if
	
      }// ijk
      
    }// else
    
    k++;
  }// j

  free(average_JK);
  
}// Jack_Knife_Average


void sigma_Jack_Knife_Average(int clusterfile, int Nt, double *array, double *sigma_JK){

  double sum = 0, sumquad = 0, factorJK_1, factorJK_2;
  double *average_JK=(double*)malloc(sizeof(double)*(Nt));
  
  factorJK_1 = (( (double) clusterfile - 2.) / ( (double) clusterfile - 1.));
  factorJK_2 = (( (double) clusterfile - 2.) / pow(( (double) clusterfile - 1.), 2));
  
  //Next block writes the JK-average for each row in average_JK[] and writes the JK-uncertainty in sigma_JK[]
  for (int t = 0; t < Nt; t++){
    for (int ijk = 0; ijk < clusterfile-1; ijk++){
      
      sum += array[ijk*Nt + t];
      sumquad += pow( array[ijk*Nt + t], 2);
      
    }// ijk
    
    average_JK[t] = sum/(clusterfile-1);
    sigma_JK[t] = sqrt((factorJK_1 * sumquad) - (factorJK_2 * pow(sum, 2)) );
    
    sum = 0;
    sumquad = 0;
  }// t

  free(average_JK);
  
}// sigma_Jack_Knife_Average


void effective_mass(int clusterfile, int Tmax, double prec, double *corr_JK, double *m_eff_array){

  int j = 0;
  double temp, u, d;
  
  double **ratio=(double**)malloc(sizeof(double)*(clusterfile));
  for(int ijk = 0; ijk < clusterfile; ijk++) ratio[ijk]=(double*)malloc(sizeof(double)*(Tmax/2));
  
  double **m_eff=(double**)malloc(sizeof(double)*(clusterfile));
  for(int ijk = 0; ijk < clusterfile; ijk++) m_eff[ijk]=(double*)malloc(sizeof(double)*(Tmax/2));
  
  //Next block fills the double array ratio[cluster][row] with the values of the ratios C[t]/C[t+1]
  for (int ijk = 0; ijk < clusterfile; ijk++){
    for (int t = 0; t < Tmax/2; t++){ 
      
      ratio[ijk][t] = corr_JK[(ijk*Tmax) + t]/corr_JK[(ijk*Tmax) + (t+1)] ;
      m_eff[ijk][t] = log (  fabs( ratio[ijk][t] ) );
      
    }// t
  }// ijk
  
  //Iterative cicle to evaluate the effective mass for each JK and each row
  for (int ijk = 0; ijk < clusterfile; ijk++){
    for (int t = 0; t < Tmax/2; t++){ 
      
      if( (ratio[ijk][t] > 0) && (m_eff[ijk][t] > 0)  ){  // Questo è stato aggiunto perché per alcune masse (in genere alte), per alti valori del tempo i correlatori assumono valori insensati
	                                                  // e "ratio" ed "m_eff" possono diventare negative. Questo succede in genere per gli ultimi 2/3 valori del tempo.
	                                                  // Comunque questi valori sono esclusi dal fit, poiché stanno chiaramente fuori dal plateau.
	do{
	  
	  u = 1 + exp( -m_eff[ijk][t]*(Tmax - 2*(t + 1)) );
	  d = 1 + exp( -m_eff[ijk][t]*(Tmax - 2*t) );
	  temp = m_eff[ijk][t];
	  
	  m_eff[ijk][t] = log( ratio[ijk][t]*(u/d) );
	  
	} while ( fabs( (m_eff[ijk][t]-temp)/m_eff[ijk][t] ) > prec );
	
      } else {
	
	if( t_max >= t ){          // ABBASSO IL VALORE t_max ENTRO CUI FARE IL FIT
	  t_max = t-1;
	  t_dat = t_max-t_min+1;
	}
	printf("WARNING: JK=%d t=%d ratio=%f m_eff=%f\n", ijk, t, ratio[ijk][t], m_eff[ijk][t]);
		
      }// else
      
    }// t
  }// ijk
  
  for (int ijk = 0; ijk < clusterfile; ijk++){
    for (int t = 0; t < Tmax/2; t++){
      
      m_eff_array[j] = m_eff[ijk][t];
      
      j++;
      
    }// t
  }// ijk
  
  free(ratio);  
  free(m_eff);

}// effective_mass


void m_eff_and_Z_check(int clusterfile, int Tmax, double *corr_JK, double *m_eff_array){
  
  double *Z_array_check=(double*)malloc(sizeof(double)*(clusterfile)*(Tmax/2));
  
  double *meff_check=(double*)malloc(sizeof(double)*(clusterfile)), *err_JK=(double*)malloc(sizeof(double)*(Tmax/2));
  double *Z_check=(double*)malloc(sizeof(double)*(clusterfile)), *err_Z=(double*)malloc(sizeof(double)*(Tmax/2));
  
  double *err_meff_check=(double*)malloc(sizeof(double)), *err_Z_check=(double*)malloc(sizeof(double));
  
  sigma_Jack_Knife_Average(clusterfile, Tmax/2, m_eff_array, err_JK);
  
  Z_function_JK( corr_JK, m_eff_array, Z_array_check, err_Z, Tmax);
  
  constant_fit_scorrelated( meff_check, m_eff_array, err_JK, clusterfile, Tmax/2, t_min, t_max-1);   // questo non so se è giusto t_max-1
  constant_fit_scorrelated( Z_check, Z_array_check, err_Z, clusterfile, Tmax/2, t_min, t_max-1);     // questo non so se è giusto t_max-1
  
  sigma_Jack_Knife_Average(clusterfile, 1, meff_check, err_meff_check);
  sigma_Jack_Knife_Average(clusterfile, 1, Z_check, err_Z_check);
  
  meff_check_glb = meff_check[clusterfile-1];
  err_meff_check_glb = err_meff_check[0];
  Z_check_glb = Z_check[clusterfile-1];
  err_Z_check_glb = err_Z_check[0];
  
  free(Z_array_check);
  free(meff_check);
  free(Z_check);
  free(err_JK);
  free(err_Z);
  free(err_meff_check);
  free(err_Z_check);
  
}// m_eff_and_Z_check


void Z_function_JK( double *corr_array, double *mass_array, double *Z_array, double *err_Z_array, int Tmax){

  for (int ijk = 0; ijk < clusterfile; ijk++){
    for (int t = 0; t < Tmax/2; t++){

      Z_array[ijk*(Tmax/2) + t]  = sqrt( fabs( (mass_array[ijk*(Tmax/2) + t] * corr_array[(ijk*Tmax) + t])/( exp(-Tmax*( fabs(mass_array[ijk*(Tmax/2) + t])/2) )*cosh((Tmax/2 - t)*fabs(mass_array[ijk*(Tmax/2) + t]) ) ) ) );

    }// t
  }// ijk

  sigma_Jack_Knife_Average(clusterfile, Tmax/2, Z_array, err_Z_array);

}// Z_function_JK

void corr_func_fit(int clusterfile, int Nt, int t_min, int t_dat, double *corr_JK, double *sigma_JK, double *mass_out, double *z_term, int ibeta, int imusea, int im1, int im2, int ith1){
  
  if(ibeta==0 && imusea==2 && im1== 0 && im2==9 && ith1==2){
    Z_check_glb = Z_check_glb*sqrt(2);
  }

  double chi2_dof;
  
  double outpar[npar], err[npar], par_chi2[npar];
  double par[npar] = {pow(Z_check_glb,2), meff_check_glb};
  double step[npar] = {0.01, 0.01};
  double min[npar] = {0.00, 0.00};
  double max[npar] = {0.00, 0.00};
  string cpar[npar] = {"Coeff", "Energy"};
  //0.0005 pow(Z_check_glb,2)
  
  for(int t = 0; t < t_dat; t++){
    
    x_fit[t] = t + t_min;
    ey_fit[t] = sigma_JK[t + t_min];
    
  }// t
  
  for(int ijk = 0; ijk < clusterfile; ijk++){
    
    for(int t = t_min; t < t_min+t_dat; t++){
      y_fit[t - t_min] = corr_JK[ijk*Nt + t];
    }// t
    
    TMinuit minuit(npar);
    
    minuit.SetFCN(chi2);
    
    minuit.SetErrorDef(1.);
    
    for(int j = 0; j < npar; j++){
      minuit.DefineParameter(j, cpar[j].c_str(), par[j], step[j], min[j], max[j]);
    }// j
    
    minuit.Migrad();
    
    for(int j = 0; j < npar; j++){
      minuit.GetParameter(j, outpar[j], err[j]);
    }// j
    
    mass_out[ijk] = outpar[1];
    z_term[ijk] = sqrt( outpar[0] );
    
    par_chi2[0] = outpar[0];
    par_chi2[1] = outpar[1];

    chi2_dof = chi_num( par_chi2);
    printf("CCCC ibeta=%d\timusea=%d\tm_%d%d\tth_%d\tnum_points=%d\tchi2[%d]=%f\n", ibeta, imusea, im1, im2, ith1, num_points, ijk, chi2_dof);

    num_points = 0;
    
  }// ijk
  
}// corr_func_fit

double corr_function( double t, double *par, int Nt){

  double Tmax = (double) Nt;
  
  return (par[0]/par[1])*exp( -(Tmax*par[1])/2.)*cosh( ((Tmax/2.)-t)*par[1] );
  
}// corr_function

void chi2( int &npar, double *deriv, double &f, double *par, int iflag){
  
  f = 0;
  for(int t = 0; t < t_dat; t++){
    f += pow( y_fit[t] - corr_function( x_fit[t], par, Tmax), 2)/pow( ey_fit[t], 2);
  }// t
  
}// chi2

double chi_num( double *parameters){
  
  double f = 0;
  
  for(int t = 0; t < t_dat; t++){
    
    f += pow( y_fit[t] - corr_function( x_fit[t], parameters, Tmax), 2)/pow( ey_fit[t], 2);

    num_points = num_points + 1;
  }// t
  
  return f/(num_points - npar);
  
}// chi_num

void output_Energy_and_Z_from_2pts(int clusterfile, double *mass_out, double *z_term, const string *beta_V, const string mu_sea_1[Nbeta][Nmusea], const string *mu_sea_2, const string *m_, const string *th_, const string *sme_, int ibeta, int imusea, int im1, int im2, int ith1, int ismear){

  const int LEN_NAME = 1024;
  
  FILE *fout_E_fit, *fout_Z, *fw_M, *fout_E_std_DR, *fout_E_sinh_DR;

  double *E_from_std_DR=(double*)malloc(sizeof(double)*(clusterfile)), *err_E_from_std_DR=(double*)malloc(sizeof(double));
  double *E_from_sinh_DR=(double*)malloc(sizeof(double)*(clusterfile)), *err_E_from_sinh_DR=(double*)malloc(sizeof(double));
  double *err_z_term=(double*)malloc(sizeof(double)), *err_mass_out=(double*)malloc(sizeof(double));

  char *output_file_E_from_fit=(char*)malloc(sizeof(char)*(LEN_NAME)), *output_file_E_from_std_DR=(char*)malloc(sizeof(char)*(LEN_NAME)), *output_file_E_from_sinh_DR=(char*)malloc(sizeof(char)*(LEN_NAME));
  char *output_file_M=(char*)malloc(sizeof(char)*(LEN_NAME)), *output_file_Z=(char*)malloc(sizeof(char)*(LEN_NAME));
  
  if( ith1 == 3 ){
    
    for(int ith = 0; ith < 3; ith++){
      
      sprintf(output_file_E_from_std_DR, "OUTPUT_SMEAR/M/fit_data/%s/Fit_%s/fit_2pts_E_from_std_DR.%s.m1m2_%s%s.th_%s.sme_30_%s0.r1r2_Av_00_11.P5P5.out", beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith].c_str(), sme_[ismear].c_str());
      sprintf(output_file_E_from_sinh_DR, "OUTPUT_SMEAR/M/fit_data/%s/Fit_%s/fit_2pts_E_from_sinh_DR.%s.m1m2_%s%s.th_%s.sme_30_%s0.r1r2_Av_00_11.P5P5.out", beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith].c_str(), sme_[ismear].c_str());
      
      if ((fout_E_std_DR = fopen(output_file_E_from_std_DR, "w")) == NULL ){
        printf("Error opening the output file: %s\n", output_file_E_from_std_DR);
        exit(EXIT_FAILURE);
      }
      if ((fout_E_sinh_DR = fopen(output_file_E_from_sinh_DR, "w")) == NULL ){
        printf("Error opening the output file: %s\n", output_file_E_from_sinh_DR);
        exit(EXIT_FAILURE);
      }
      
      for(int ijk = 0; ijk < clusterfile; ijk++){
	
	E_from_std_DR[ijk] = stdDR(mass_out[ijk], ith);
	E_from_sinh_DR[ijk] = sinh_DR(mass_out[ijk], ith);
	
      }// ijk
      
      sigma_Jack_Knife_Average(clusterfile, 1, E_from_std_DR, err_E_from_std_DR);
      sigma_Jack_Knife_Average(clusterfile, 1, E_from_sinh_DR, err_E_from_sinh_DR);
      
      for(int ijk = 0; ijk < clusterfile; ijk++){
	
	fprintf(fout_E_std_DR, "%+016.16g\n", E_from_std_DR[ijk]);
	fprintf(fout_E_sinh_DR,"%+016.16g\n", E_from_sinh_DR[ijk]);
	
      }// ijk
      
      fprintf(fout_E_std_DR, "JK_AVERAGE: E_from_std_DR, err_JK_E_std\n%+016.16g\t%+016.16g\n", E_from_std_DR[15], err_E_from_std_DR[0]);
      fprintf(fout_E_sinh_DR, "JK_AVERAGE: E_from_sinh_DR, err_JK_E_sinh\n%+016.16g\t%+016.16g\n", E_from_sinh_DR[15], err_E_from_sinh_DR[0]);
      
      fclose(fout_E_std_DR);
      fclose(fout_E_sinh_DR);
      
    }// for( ith = 0; ith < 3; ith++)

    sprintf(output_file_M, "OUTPUT_SMEAR/M/fit_data/%s/Fit_%s/fit_2pts_M.%s.m1m2_%s%s.th_%s.sme_30_%s0.r1r2_Av_00_11.P5P5.out", beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith1].c_str(), sme_[ismear].c_str());
    sprintf(output_file_Z, "OUTPUT_SMEAR/M/fit_data/%s/Fit_%s/fit_2pts_Z.%s.m1m2_%s%s.th_%s.sme_30_%s0.r1r2_Av_00_11.P5P5.out", beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith1].c_str(), sme_[ismear].c_str());
    
    if ((fw_M = fopen(output_file_M, "w")) == NULL ){
      printf("Error opening the output file!!\n");
      exit(EXIT_FAILURE);
    }
    if ((fout_Z = fopen(output_file_Z, "w")) == NULL ){
      printf("Error opening the output file!!\n");
      exit(EXIT_FAILURE);
    }

    sigma_Jack_Knife_Average(clusterfile, 1, mass_out, err_mass_out);
    sigma_Jack_Knife_Average(clusterfile, 1, z_term, err_z_term);
    
    for(int ijk = 0; ijk < clusterfile; ijk++){
      
      fprintf(fw_M, "%+016.16g\n", mass_out[ijk]);
      fprintf(fout_Z, "%+016.16g\n", z_term[ijk]);
      
    }// ijk
    
    fprintf(fw_M, "JK_AVERAGE: MASS, err_JK_MASS, Z, err_JK_Z\n%+016.16g\t%+016.16g\t%+016.16g\t%+016.16g\n", mass_out[15], err_mass_out[0], z_term[15], err_z_term[0]);
    fprintf(fout_Z, "JK_AVERAGE: MASS, err_JK_MASS, Z, err_JK_Z\n%+016.16g\t%+016.16g\t%+016.16g\t%+016.16g\n", mass_out[15], err_mass_out[0], z_term[15], err_z_term[0]);
    
  }else{  // SE ith1 != 3

    sprintf(output_file_E_from_fit, "OUTPUT_SMEAR/M/fit_data/%s/Fit_%s/fit_2pts_E.%s.m1m2_%s%s.th_%s.sme_30_%s0.r1r2_Av_00_11.P5P5.out", beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith1].c_str(), sme_[ismear].c_str());
    sprintf(output_file_Z, "OUTPUT_SMEAR/M/fit_data/%s/Fit_%s/fit_2pts_Z.%s.m1m2_%s%s.th_%s.sme_30_%s0.r1r2_Av_00_11.P5P5.out", beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith1].c_str(), sme_[ismear].c_str());
    
    if ((fout_E_fit = fopen(output_file_E_from_fit, "w")) == NULL ){
      printf("Error opening the output file!!\n");
      exit(EXIT_FAILURE);
    }
    if ((fout_Z = fopen(output_file_Z, "w")) == NULL ){
      printf("Error opening the output file!!\n");
      exit(EXIT_FAILURE);
    }

    sigma_Jack_Knife_Average(clusterfile, 1, mass_out, err_mass_out);
    sigma_Jack_Knife_Average(clusterfile, 1, z_term, err_z_term);
    
    for(int ijk = 0; ijk < clusterfile; ijk++){
      
      fprintf(fout_E_fit, "%+016.16g\n", mass_out[ijk]);
      fprintf(fout_Z, "%+016.16g\n", z_term[ijk]);
      
    }// ijk
    
    fprintf(fout_E_fit, "JK_AVERAGE: E_from_fit, err_JK_E, Z, err_JK_Z\n%+016.16g\t%+016.16g\t%+016.16g\t%+016.16g\n", mass_out[15], err_mass_out[0], z_term[15], err_z_term[0]);
    fprintf(fout_Z, "JK_AVERAGE: E_from_fit, err_JK_E, Z, err_JK_Z\n%+016.16g\t%+016.16g\t%+016.16g\t%+016.16g\n", mass_out[15], err_mass_out[0], z_term[15], err_z_term[0]);
    
    fclose(fout_E_fit);
    fclose(fout_Z);
    
  }// else

  free(err_mass_out);
  free(err_z_term);
  free(E_from_std_DR);
  free(E_from_sinh_DR);
  free(err_E_from_std_DR);
  free(err_E_from_sinh_DR);
  free(output_file_E_from_fit);
  free(output_file_E_from_std_DR);
  free(output_file_E_from_sinh_DR);
  free(output_file_M);
  free(output_file_Z);
  
  // PRINT ON synthetic_log FOR CHECK
  FILE *fw_log;
  
  if ((fw_log = fopen("OUTPUT_SMEAR/M/LOG/synthetic_log.out", "a")) == NULL ){
    printf("Error opening the output file!!\n");
    exit(EXIT_FAILURE);
  }
  
  fprintf(fw_log, "EEEE  %s\t%d\t%d%d\t%d\t%+016.16g\t%+016.16g\n", beta_V[ibeta].c_str(), imusea, im1, im2, ith1, mass_out[15], err_mass_out[0]);     
  fprintf(fw_log, "CCCE  %f %f\t%f %f %f\n", mass_out[15], meff_check_glb, err_mass_out[0], err_meff_check_glb, err_mass_out[0]/err_meff_check_glb);
  fprintf(fw_log, "CCCZ  %f %f\t%f %f %f\n", z_term[15], Z_check_glb, err_z_term[0], err_Z_check_glb, err_z_term[0]/err_Z_check_glb);
  fprintf(fw_log, "ZZZZ  %s\t%d\t%d%d\t%d\t%+016.16g\t%+016.16g\n", beta_V[ibeta].c_str(), imusea, im1, im2, ith1, z_term[15], err_z_term[0]);
  
  fclose(fw_log);

  // PRINT ON STANDARD OUTPUT
  //printf("\n\nENERGY=%f  Z_CHECK=%f\n\n", meff_check_glb, pow(Z_check_glb,2));
  
  printf("RRRR ibeta=%d\timusea=%d\tm_%d%d\tth_%d\tE_check=%f\tZ_check=%f\n", ibeta, imusea, im1, im2, ith1, meff_check_glb, pow(Z_check_glb,2));
  
}// output_Energy_and_Z_from_2pts

void out_eff_mass_graph_grace(int clusterfile, int Tmax, int t_min, int t_max, double *m_eff_array, double *mass_out, double *sigma_mass_out, const string *beta_V, const string mu_sea_1[Nbeta][Nmusea], const string *mu_sea_2, const string *m_, const string *th_, const string *sme_, int ibeta, int imusea, int im1, int im2, int ith1, int ismear){
  
  FILE *fout_Graph;
  
  char *output_graph_grace=(char*)malloc(sizeof(char)*(LEN_NAME));
  
  double *sigma_m_eff_array=(double*)malloc(sizeof(double)*(Tmax/2));
  
  sigma_Jack_Knife_Average(clusterfile, Tmax/2, m_eff_array, sigma_m_eff_array);
  
  sprintf(output_graph_grace, "OUTPUT_SMEAR/M/Plateaux/%s/Gr_%s/plateau_grace_2pts.%s.m1m2_%s%s.th_%s.sme_30_%s0.r1r2_Av_00_11.P5P5.xmg", beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith1].c_str(), sme_[ismear].c_str());
  
  if ((fout_Graph = fopen(output_graph_grace, "w")) == NULL ){
    printf("Error opening the output file!!\n");
    exit(EXIT_FAILURE);
  }
  
  fprintf(fout_Graph,"@type xydy\n");
  
  for(int t = 0; t < Tmax/2; t++ ){
    fprintf(fout_Graph,"%d %f %f\n", t, m_eff_array[(clusterfile-1)*(Tmax/2) + t], sigma_m_eff_array[t]);
  }// t
  
  fprintf(fout_Graph, "@type xy\n");
  fprintf(fout_Graph, "%d %f\n", t_min, mass_out[clusterfile-1]);
  fprintf(fout_Graph, "%d %f\n", t_max, mass_out[clusterfile-1]);
  
  fprintf(fout_Graph, "@type xy\n");
  fprintf(fout_Graph, "%d %f\n", t_min, mass_out[clusterfile-1]+sigma_mass_out[0]);
  fprintf(fout_Graph, "%d %f\n", t_max, mass_out[clusterfile-1]+sigma_mass_out[0]);
  
  fprintf(fout_Graph, "@type xy\n");
  fprintf(fout_Graph, "%d %f\n", t_min, mass_out[clusterfile-1]-sigma_mass_out[0]);
  fprintf(fout_Graph, "%d %f\n", t_max, mass_out[clusterfile-1]-sigma_mass_out[0]);
  
  fclose(fout_Graph);
  
  free(sigma_m_eff_array);
  free(output_graph_grace);
  
}// out_eff_mass_graph_grace

double stdDR(double mass_out, int ith){

  double th_val, momentum, E;
  
  th_val = theta_value[ibeta][6-ith];
  momentum = (th_val*PI)/L[ibeta];

  E = sqrt( pow(mass_out, 2) + 3*pow(momentum, 2) );

  return E;
}

double sinh_DR(double mass_out, int ith){

  double th_val, momentum, E;
  
  th_val = theta_value[ibeta][6-ith];
  momentum = (th_val*PI)/L[ibeta];

  E = 2*my_arcsinh( sqrt( pow( sinh(mass_out/2) , 2) + 3*pow( sin(momentum/2) , 2) ) );

  return E;
}

double my_arcsinh(double x){
  return log( x + sqrt(pow(x,2) + 1)  );
}

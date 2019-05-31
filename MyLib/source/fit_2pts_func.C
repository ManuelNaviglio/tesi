#include <TF1.h>
#include <fstream>
#include "clean_str_cl.h"
#include "definitions.h"

using namespace std;

void read_time_intervals( string *input_file, int &t_min, int &t_max, int &t_dat, int im1, int ibeta, int im2){
  
  int ind;
  int *minimum_x=(int*)malloc(sizeof(int)*(Nmass_m1_time_intervals)*(Nbeta)*(Nmasses));
  int *maximum_x=(int*)malloc(sizeof(int)*(Nmass_m1_time_intervals)*(Nbeta)*(Nmasses));
  int *data=(int*)malloc(sizeof(int)*(Nmass_m1_time_intervals)*(Nbeta)*(Nmasses));
  
  string line;
  
  ifstream f_in (input_file[0].c_str());
  if (!f_in) cout<< "Cannot open "<< input_file[0].c_str() <<endl;
  
  for(int j = 0; j < Nmass_m1_time_intervals; j++){
    
    getline(f_in, line);  // #MASS m1=
    cout<<line<<endl;
    
    for(int i = 0; i < Nbeta; i++){
      
      getline(f_in, line);  // beta= L=
      cout<<line<<endl;
      
      for(int k = 0; k < Nmasses; k++){

	ind = (j*Nbeta + i) + (Nmass_m1_time_intervals*Nbeta)*k;
	
	f_in >> minimum_x[ind];
	f_in >> maximum_x[ind];
	f_in >> data[ind];

	printf("%d %d %d\n", minimum_x[ind], maximum_x[ind], data[ind]);
	
      }// k
      getline(f_in, line);  // newline
      
    }// i
    
  }// j
  f_in.close();

  ind = (im1*Nbeta + ibeta) + (Nmass_m1_time_intervals*Nbeta)*im2;
  
  *&t_min = minimum_x[ind];
  *&t_max = maximum_x[ind];
  *&t_dat = data[ind];
  /*
  printf("ibeta=%d\tim1=%d\tim2=%d\n", ibeta, im1, im2);

  printf("\n\nTIME_INTERVAL: tmin=%d tmax=%d tdat=%d\n\n", t_min, t_max, t_dat);
  */
  free(minimum_x);
  free(maximum_x);
  free(data);
  
}// read_time_intervals

void r_average(double *array1, double *array2, double *average, int Nrows_file){
  
  for(int i = 0; i < Nrows_file; i++){
    *(average + i) = (( *(array1 + i) + *(array2 + i) ) / 2 );
  }
}// r_average


void th_r_average(double *array1, double *array2, double *array3, double *array4, double *average, int Nrows_file){
  
  for(int i = 0; i < Nrows_file; i++){
    *(average + i) = (( *(array1 + i) + *(array2 + i) + *(array3 + i) + *(array4 + i) ) / 4 );
  }
}// th_r_average

void simmetrize(double *array, double *simmetrized, int Tmax, int clusterfile){

  for (int n = 0; n < clusterfile; n++){
    for (int j = 0; j < Tmax; j++){
      if(j == 0){
        simmetrized[(n * Tmax) + j] = array[(n * Tmax) + j];
      }else{
        simmetrized[(n * Tmax) + j] = (array[(n * Tmax) + j] + array[((n + 1) * Tmax) - j])/2;
      }
    }
  }
}// simmetrize


void read_corr_2pts_from_2pts( double *corr_2pts, int Nt, int clusterfile, const string *beta_V, const string mu_sea_1[Nbeta][Nmusea], const string n_conf[Nbeta][Nmusea], const string *mu_sea_2, const string *m_, const string *th_, const string *sme_, int ibeta, int imusea, int im1, int im2, int ith, int ismear){

  int LEN_NAME = 1024;
  
  FILE *fr_2pts_th_r1r2_00, *fr_2pts_th_r1r2_11;  
  char *open_2pts_th_r1r2_00=(char*)malloc(sizeof(char)*LEN_NAME), *open_2pts_th_r1r2_11=(char*)malloc(sizeof(char)*LEN_NAME);
  
  double *temp_2pts_th_r1r2_00=(double*)malloc(sizeof(double)*(clusterfile)*(Nt)), *temp_2pts_th_r1r2_11=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  double *temp_corr_2pts=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
  
  if(ith != 3){

    FILE *fr_2pts_th_rev_r1r2_00, *fr_2pts_th_rev_r1r2_11;
    char *open_2pts_th_rev_r1r2_00=(char*)malloc(sizeof(char)*LEN_NAME), *open_2pts_th_rev_r1r2_11=(char*)malloc(sizeof(char)*LEN_NAME);
    double *temp_2pts_th_rev_r1r2_00=(double*)malloc(sizeof(double)*(clusterfile)*(Nt)), *temp_2pts_th_rev_r1r2_11=(double*)malloc(sizeof(double)*(clusterfile)*(Nt));
    
    // CASO ith != 3
    sprintf(open_2pts_th_r1r2_00, "DATA_SMEAR/%s/%s/2pts.%s.%s.m1m2_%s%s.th_%s.sme_30_%s0.r1r2_00.P5P5.out", beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), n_conf[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith].c_str(), sme_[ismear].c_str());
    sprintf(open_2pts_th_r1r2_11, "DATA_SMEAR/%s/%s/2pts.%s.%s.m1m2_%s%s.th_%s.sme_30_%s0.r1r2_11.P5P5.out", beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), n_conf[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith].c_str(), sme_[ismear].c_str());
    sprintf(open_2pts_th_rev_r1r2_00, "DATA_SMEAR/%s/%s/2pts.%s.%s.m1m2_%s%s.th_%s.sme_30_%s0.r1r2_00.P5P5.out", beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), n_conf[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[6-ith].c_str(), sme_[ismear].c_str());
    sprintf(open_2pts_th_rev_r1r2_11, "DATA_SMEAR/%s/%s/2pts.%s.%s.m1m2_%s%s.th_%s.sme_30_%s0.r1r2_11.P5P5.out", beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), n_conf[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[6-ith].c_str(), sme_[ismear].c_str());
    

    if ((fr_2pts_th_r1r2_00 = fopen(open_2pts_th_r1r2_00, "r")) == NULL ){
      printf("Error opening the file to read: 2pts_0m1_th1_r1r2_00\n");
      exit(EXIT_FAILURE);
    }
    if ((fr_2pts_th_r1r2_11 = fopen(open_2pts_th_r1r2_11, "r")) == NULL ){
      printf("Error opening the file to read: 2pts_0m1_th1_r1r2_11\n");
      exit(EXIT_FAILURE);
    }
    if ((fr_2pts_th_rev_r1r2_00 = fopen(open_2pts_th_rev_r1r2_00, "r")) == NULL ){
      printf("Error opening the file to read: 2pts_0m1_th1_rev_r1r2_00\n");
      exit(EXIT_FAILURE);
    }
    if ((fr_2pts_th_rev_r1r2_11 = fopen(open_2pts_th_rev_r1r2_11, "r")) == NULL ){
      printf("Error opening the file to read: 2pts_0m1_th1_rev_r1r2_11\n");
      exit(EXIT_FAILURE);
    }
    free(open_2pts_th_rev_r1r2_00);
    free(open_2pts_th_rev_r1r2_11);
  
    clean_string_cluster(fr_2pts_th_r1r2_00, temp_2pts_th_r1r2_00, Nt, clusterfile);
    clean_string_cluster(fr_2pts_th_r1r2_11, temp_2pts_th_r1r2_11, Nt, clusterfile);
    clean_string_cluster(fr_2pts_th_rev_r1r2_00, temp_2pts_th_rev_r1r2_00, Nt, clusterfile);
    clean_string_cluster(fr_2pts_th_rev_r1r2_11, temp_2pts_th_rev_r1r2_11, Nt, clusterfile);
    
    fclose(fr_2pts_th_r1r2_00);
    fclose(fr_2pts_th_r1r2_11);
    fclose(fr_2pts_th_rev_r1r2_00);
    fclose(fr_2pts_th_rev_r1r2_11);
    
    th_r_average(temp_2pts_th_r1r2_00, temp_2pts_th_r1r2_11, temp_2pts_th_rev_r1r2_00, temp_2pts_th_rev_r1r2_11, temp_corr_2pts, clusterfile*Nt);
    
    simmetrize(temp_corr_2pts, corr_2pts, Nt, clusterfile);
    
    free(temp_2pts_th_rev_r1r2_00);
    free(temp_2pts_th_rev_r1r2_11);
    
  }// ith != 3

  if(ith == 3){

    // CASO ith1 == 3
    sprintf(open_2pts_th_r1r2_00, "DATA_SMEAR/%s/%s/2pts.%s.%s.m1m2_%s%s.th_%s.sme_30_%s0.r1r2_00.P5P5.out", beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), n_conf[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith].c_str(), sme_[ismear].c_str());
    sprintf(open_2pts_th_r1r2_11, "DATA_SMEAR/%s/%s/2pts.%s.%s.m1m2_%s%s.th_%s.sme_30_%s0.r1r2_11.P5P5.out", beta_V[ibeta].c_str(), mu_sea_1[ibeta][imusea].c_str(), n_conf[ibeta][imusea].c_str(), mu_sea_2[imusea].c_str(), m_[im1].c_str(), m_[im2].c_str(), th_[ith].c_str(), sme_[ismear].c_str());
    
    if ((fr_2pts_th_r1r2_00 = fopen(open_2pts_th_r1r2_00, "r")) == NULL ){
      printf("Error opening the file to read: 2pts_0m1_th1_r1r2_00\n");
      exit(EXIT_FAILURE);
    }
    if ((fr_2pts_th_r1r2_11 = fopen(open_2pts_th_r1r2_11, "r")) == NULL ){
      printf("Error opening the file to read: 2pts_0m1_th1_r1r2_11\n");
      exit(EXIT_FAILURE);
    }
    
    clean_string_cluster(fr_2pts_th_r1r2_00, temp_2pts_th_r1r2_00, Nt, clusterfile);
    clean_string_cluster(fr_2pts_th_r1r2_11, temp_2pts_th_r1r2_11, Nt, clusterfile);
    
    fclose(fr_2pts_th_r1r2_00);
    fclose(fr_2pts_th_r1r2_11);
    
    r_average(temp_2pts_th_r1r2_00, temp_2pts_th_r1r2_11, temp_corr_2pts, clusterfile*Nt);
    
    simmetrize(temp_corr_2pts, corr_2pts, Nt, clusterfile);
    
  }// ith == 3
  
  free(open_2pts_th_r1r2_00);
  free(open_2pts_th_r1r2_11);
  
  free(temp_2pts_th_r1r2_00);
  free(temp_2pts_th_r1r2_11);
  free(temp_corr_2pts);

} // read_corr_2pts_from_2pts

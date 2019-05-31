#include <TF1.h>            // questo serve per la funzione cout
#include <fstream>

using namespace std;

double sigma_std(double array[], int i_in, int i_fin){

  double sum = 0, sumquad = 0;
  double sigma_JK;

    for (int n = i_in; n <= i_fin; n++){

      sum += array[n]/(i_fin - i_in + 1);
      sumquad += pow( array[n], 2)/(i_fin - i_in + 1);

    }

    sigma_JK = sqrt( sumquad - pow( sum ,2) );

    return sigma_JK;
  
}


double sigma_JK(double array[], int clusterfile){

  double sum = 0, sumquad = 0, factorJK_1, factorJK_2;
  double sigma_JK;
  factorJK_1 = (( (double) clusterfile - 2.) / ( (double) clusterfile - 1.));
  factorJK_2 = (( (double) clusterfile - 2.) / pow(( (double) clusterfile - 1.), 2));

    for (int n = 0; n < clusterfile-1; n++){
      sum += *(array + n) ;
      sumquad += pow( *(array + n), 2);
    }

    sigma_JK = sqrt((factorJK_1 * sumquad) - (factorJK_2 * pow(sum, 2)) );

    return sigma_JK;
  
}// sigma_JK



double sigma_JK_modified_2(double array[], int i_in, int i_fin, int clusterfile){

  double sum = 0, sumquad = 0;
  double sigma_JK;

    for (int n = i_in; n <= i_fin; n++){

      sum += array[n]/(i_fin - i_in + 1);
      sumquad += pow( array[n], 2)/(i_fin - i_in + 1);

    }

    sigma_JK = sqrt( ( (double) clusterfile - 2. )*( sumquad - pow( sum ,2)  ) );

    return sigma_JK;
  
}



double sigma_bootstrap(double array[], int analysis_in, int analysis_fin, int Nev_an, int clusterfile){
  
  int num_analysis = (analysis_fin - analysis_in) + 1;
  
  int start_ev = (analysis_in - 1)*Nev_an + 1, end_ev = analysis_fin*Nev_an;
  
  double* sigma = (double*)calloc( num_analysis, sizeof(double));
  double* average = (double*)calloc( num_analysis, sizeof(double));
  
  int ian, Nev_real = num_analysis*Nev_an;
  
  double sigma_boot = 0, average_tot = 0;
  
  for(int iev = start_ev; iev <= end_ev; iev++){   // NOTA: L'ev 0 NON MI SERVE
    
    ian = ((iev-1)/Nev_an);
    
    average[ian-(analysis_in-1)] += array[iev]/Nev_an;

    average_tot += array[iev]/Nev_real;
    
  }
  
  
  for(int ianalysis = analysis_in-1; ianalysis <= analysis_fin-1; ianalysis++){
    
    sigma[ianalysis-(analysis_in-1)] = sigma_JK_modified_2( array, ianalysis*Nev_an + 1, ianalysis*Nev_an + Nev_an, clusterfile);
    
    sigma_boot += ( pow( sigma[ianalysis-(analysis_in-1)] ,2) + pow( average[ianalysis-(analysis_in-1)] - average_tot ,2) )/num_analysis;

  }
  free(sigma);
  free(average);

  return sqrt(sigma_boot);

}

double statistical_sigma_bootstrap(double array[], int analysis_in, int analysis_fin, int Nev_an, int clusterfile){

  int num_analysis = (analysis_fin - analysis_in) + 1;
  
  double* sigma = (double*)calloc( num_analysis, sizeof(double));

  double stat_sigma_boot = 0;

  for(int ianalysis = analysis_in-1; ianalysis <= analysis_fin-1; ianalysis++){

    sigma[ianalysis-(analysis_in-1)] = sigma_JK_modified_2( array, ianalysis*Nev_an + 1, ianalysis*Nev_an + Nev_an, clusterfile);

    stat_sigma_boot += pow( sigma[ianalysis-(analysis_in-1)] ,2)/num_analysis;

  }
  free(sigma);
  
  return sqrt(stat_sigma_boot);

}

double systematic_sigma_bootstrap(double array[], int analysis_in, int analysis_fin, int Nev_an){

  int num_analysis = (analysis_fin - analysis_in) + 1;

  int start_ev = (analysis_in - 1)*Nev_an + 1, end_ev = analysis_fin*Nev_an;
  
  double* average = (double*)calloc( num_analysis, sizeof(double));

  int ian, Nev_real = num_analysis*Nev_an;

  double syst_sigma_boot = 0, average_tot = 0;

  for(int iev = start_ev; iev <= end_ev; iev++){   // NOTA: L'ev 0 NON MI SERVE

      ian = ((iev-1)/Nev_an);

      average[ian-(analysis_in-1)] += array[iev]/Nev_an;

      average_tot += array[iev]/Nev_real;

  }

  for(int ianalysis = analysis_in-1; ianalysis <= analysis_fin-1; ianalysis++){

    syst_sigma_boot += pow( average[ianalysis-(analysis_in-1)] - average_tot ,2)/num_analysis;

  }
  free(average);
  
  return sqrt(syst_sigma_boot);

}


double check_sigma_bootstrap(double array[], int analysis_in, int analysis_fin, int Nev_an, int clusterfile){

  double stat_sigma_boot = 0, syst_sigma_boot = 0, sigma_boot;

  stat_sigma_boot = statistical_sigma_bootstrap( array, analysis_in, analysis_fin, Nev_an, clusterfile);

  syst_sigma_boot = systematic_sigma_bootstrap( array, analysis_in, analysis_fin, Nev_an);

  sigma_boot = sqrt( pow( stat_sigma_boot, 2) + pow( syst_sigma_boot, 2) );

  return sigma_boot;

}


double covariance(double array_1[], double array_2[], int i_in, int i_fin, int clusterfile){

  double av_1=0, av_2=0, av_12=0;
  double Nev=(i_fin - i_in + 1);
  double cov;
  
  for(int iev = i_in; iev <= i_fin; iev++){

    av_1 = av_1 + array_1[iev]/Nev;
    av_2 = av_2 + array_2[iev]/Nev;
    av_12 = av_12 + (array_1[iev]*array_2[iev])/Nev;
  }

  cov = ( (double) clusterfile - 2. )*( av_12-(av_1*av_2) );

  return cov;
  
}


double cov_bootstrap(double array_1[], double array_2[], int analysis_in, int analysis_fin, int Nev_an, int clusterfile){

  int num_analysis = (analysis_fin - analysis_in) + 1;

  double* av_1 = (double*)malloc(sizeof(double)*num_analysis);
  double* av_2 = (double*)malloc(sizeof(double)*num_analysis);
  double* cov_an = (double*)malloc(sizeof(double)*num_analysis);

  double av_cov=0, av_diff=0, temp_av_1=0, temp_av_2=0, av_1_tot=0, av_2_tot=0;
  
  double covariance_boot;


  for(int ianalysis = (analysis_in-1); ianalysis < num_analysis; ianalysis++){

    for(int iev = ianalysis*Nev_an + 1; iev <= (ianalysis+1)*Nev_an; iev++){

      temp_av_1 = temp_av_1 + array_1[iev]/Nev_an;
      temp_av_2 = temp_av_2 + array_2[iev]/Nev_an;

    } // iev

    av_1[ianalysis] = temp_av_1;
    av_2[ianalysis] = temp_av_2;
    cov_an[ianalysis] = covariance( array_1, array_2, ianalysis*Nev_an + 1, (ianalysis+1)*Nev_an, clusterfile);
    
    av_1_tot = av_1_tot + av_1[ianalysis]/num_analysis;
    av_2_tot = av_2_tot + av_2[ianalysis]/num_analysis;
    
    temp_av_1 = 0;
    temp_av_2 = 0;

  } // ianalysis


  for(int ianalysis = (analysis_in-1); ianalysis < num_analysis; ianalysis++){

    av_diff = av_diff + ( (av_1[ianalysis]-av_1_tot)*(av_2[ianalysis]-av_2_tot) )/num_analysis;

    av_cov = av_cov + cov_an[ianalysis]/num_analysis;
    
  }

  covariance_boot = av_cov + av_diff;


  free(av_1);
  free(av_2);
  free(cov_an);
  
  return covariance_boot;
  
}








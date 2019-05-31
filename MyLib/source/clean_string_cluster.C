#include <stdio.h>
#include <TF1.h>

void clean_string_cluster(FILE *fr, double *corr, int Tmax, int clusterfile){

  int i, j = 0, k, num_rows_input;
  int *time = (int*)malloc(sizeof(int)*(Tmax)*(clusterfile));
  char *res, C[200];
  double x;
  long int position_in_file;
  long int start_file;

  start_file = ftell(fr);
  num_rows_input = (Tmax + 1) * clusterfile;

  for (i = 1; i <= num_rows_input; i++){
    position_in_file = ftell(fr);
    res = fgets(C, 200, fr);
    if (*C == 'c' || *C == 'C' ) {
    } else {  
      fseek(fr, position_in_file, start_file);   
      fscanf(fr, "%d %lf\n", &k, &x);// Per questo \n ci ho perso una mattinata!!!
      *(time + j) = k;
      *(corr + j) = x;
      ++j;
    }
  }
  rewind(fr);

  free(time);
  
}

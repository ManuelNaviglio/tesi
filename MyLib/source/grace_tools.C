#include <TF1.h> 

void print_band_grace( double x_min, double x_max, int Nstep, std::string name_variable, double* array, double* sigma_array, int i_in){

	double *x_step=(double*)calloc( 1, sizeof(double));
  double *array_inf=(double*)malloc(sizeof(double)*(Nstep)), *array_sup=(double*)malloc(sizeof(double)*(Nstep));

  printf("#%s\n", name_variable.c_str());
  if(i_in == 0){
    printf("@type xy\n");
  }else if(i_in == 1){
    printf("&\n");
  }// if i_in

  for(int istep = 0; istep < Nstep; istep++){
		array_inf[istep] = array[istep] - sigma_array[istep];
		array_sup[istep] = array[istep] + sigma_array[istep];
  }// istep

  for(int istep = 0; istep < Nstep; istep++){
    *x_step = x_min + ( (x_max - x_min)/( (double) (Nstep - 1)))*istep;
    printf("%f %f\n", *x_step, array_sup[istep]);
  }// istep

  for(int istep = Nstep-1; istep >= 0; istep--){
    *x_step = x_min + ( (x_max - x_min)/( (double) (Nstep - 1)))*istep;
    printf("%f %f\n", *x_step, array_inf[istep]);
  }// istep

  free(x_step);
  free(array_inf);
  free(array_sup);

}// print_band_grace


void print_points_grace( double* x, double* y, int i_min, int i_max, std::string name_variable, int i_in){

  printf("#%s\n", name_variable.c_str());
  if(i_in == 0){
    printf("@type xy\n");
  }else if(i_in == 1){
    printf("&\n");
  }// if i_in

  for(int i = i_min; i <= i_max; i++){
    printf("%f %f\n", x[i], y[i]);
  }// i

}// print_points_grace


void print_points_with_err_grace( double* x, double* y, double* dy, int i_min, int i_max, std::string name_variable, int i_in){

  printf("#%s\n", name_variable.c_str());
  if(i_in == 0){
    printf("@type xydy\n");
  }else if(i_in == 1){
    printf("&\n");
  }// if i_in

  for(int i = i_min; i <= i_max; i++){
    printf("%f %f %f\n", x[i], y[i], dy[i]);
  }// i

}// print_points_with_err_grace


void print_vertical_line_grace( double x, double y_min, double y_max, std::string name_variable, int i_in){

  printf("#%s\n", name_variable.c_str());
  if(i_in == 0){
    printf("@type xy\n");
  }else if(i_in == 1){
    printf("&\n");
  }// if i_in

  printf("%f %f\n", x, y_min);
  printf("%f %f\n", x, y_max);

}// print_vertical_line_grace


void print_horizontal_line_grace( double x_min, double x_max, double y, std::string name_variable, int i_in){

  printf("#%s\n", name_variable.c_str());
  if(i_in == 0){
    printf("@type xy\n");
  }else if(i_in == 1){
    printf("&\n");
  }// if i_in

  printf("%f %f\n", x_min, y);
  printf("%f %f\n", x_max, y);

}// print_horizontal_line_grace


void print_horizontal_mean_with_err_grace( double x_min, double x_max, double y, double dy, std::string name_variable, int i_in){

	print_horizontal_line_grace( x_min, x_max, y+dy, name_variable, i_in);
	print_horizontal_line_grace( x_min, x_max, y-dy,            "", 1);
	print_horizontal_line_grace( x_min, x_max,    y,            "", 1);

}// print_horizontal_mean_with_err_grace


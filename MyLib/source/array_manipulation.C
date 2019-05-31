#include <TF1.h>            // questo serve per la funzione cout

using namespace std;

void append_two_arrays(double* array_tot, double* array_1, double* array_2, int NL1, int NL2, int Nev_arrays){

  double NL_tot = NL1+NL2;

	for(int i = 0; i < NL_tot; i++){
		array_tot[i*(Nev_arrays+1)+0] = 0.0;
	}// i

	if(Nev_arrays >= 1){
 
		for(int i = 0; i < NL_tot; i++){

			if(i < NL1){

				for(int iev = 1; iev <= Nev_arrays; iev++){
					array_tot[i*(Nev_arrays+1)+iev] = array_1[i*(Nev_arrays+1)+iev];
					array_tot[i*(Nev_arrays+1)+0] += array_tot[i*(Nev_arrays+1)+iev]/Nev_arrays;
				}// iev

			}else if(i >= NL1){

				for(int iev = 1; iev <= Nev_arrays; iev++){
					array_tot[i*(Nev_arrays+1)+iev] = array_2[(i-NL1)*(Nev_arrays+1)+iev];
					array_tot[i*(Nev_arrays+1)+0] += array_tot[i*(Nev_arrays+1)+iev]/Nev_arrays;
				}// iev

			}// if

		}// i

  }else if(Nev_arrays == 0){

		for(int i = 0; i < NL_tot; i++){
			if(i < NL1){
				array_tot[i] = array_1[i];
		  }else if(i >= NL1){
				array_tot[i] = array_2[i-NL1];
			}// if
		}// i

  }// if Nev_arrays

}// append_two_arrays


void append_three_arrays(double* array_tot, double* array_1, double* array_2, double* array_3, int NL1, int NL2, int NL3, int Nev_arrays){

  double *temp_two_arrays=(double*)malloc(sizeof(double)*(NL1+NL2)*(Nev_arrays+1));

	append_two_arrays( temp_two_arrays, array_1, array_2, NL1, NL2, Nev_arrays);

	append_two_arrays( array_tot, temp_two_arrays, array_3, NL1+NL2, NL3, Nev_arrays);

  free(temp_two_arrays);

}// append_three_arrays


void append_four_arrays(double* array_tot, double* array_1, double* array_2, double* array_3, double* array_4, int NL1, int NL2, int NL3, int NL4, int Nev_arrays){

  double *temp_three_arrays=(double*)malloc(sizeof(double)*(NL1+NL2+NL3)*(Nev_arrays+1));

	append_three_arrays( temp_three_arrays, array_1, array_2, array_3, NL1, NL2, NL3, Nev_arrays);

	append_two_arrays( array_tot, temp_three_arrays, array_4, NL1+NL2+NL3, NL4, Nev_arrays);

  free(temp_three_arrays);

}// append_four_arrays


void append_five_arrays(double* array_tot, double* array_1, double* array_2, double* array_3, double* array_4, double* array_5, int NL1, int NL2, int NL3, int NL4, int NL5, int Nev_arrays){

  double *temp_four_arrays=(double*)malloc(sizeof(double)*(NL1+NL2+NL3+NL4)*(Nev_arrays+1));

	append_four_arrays( temp_four_arrays, array_1, array_2, array_3, array_4, NL1, NL2, NL3, NL4, Nev_arrays);

	append_two_arrays( array_tot, temp_four_arrays, array_5, NL1+NL2+NL3+NL4, NL5, Nev_arrays);

  free(temp_four_arrays);

}// append_five_arrays


void append_six_arrays(double* array_tot, double* array_1, double* array_2, double* array_3, double* array_4, double* array_5, double* array_6, int NL1, int NL2, int NL3, int NL4, int NL5, int NL6, int Nev_arrays){

  double *temp_five_arrays=(double*)malloc(sizeof(double)*(NL1+NL2+NL3+NL4+NL5)*(Nev_arrays+1));

	append_five_arrays( temp_five_arrays, array_1, array_2, array_3, array_4, array_5, NL1, NL2, NL3, NL4, NL5, Nev_arrays);

	append_two_arrays( array_tot, temp_five_arrays, array_6, NL1+NL2+NL3+NL4+NL5, NL6, Nev_arrays);

  free(temp_five_arrays);

}// append_six_arrays




void sort_array(double* array_sort, double* array_in, int* conv_array, int Nq, int Nev_tot, int Nev_block, int Nev_an){

 	double Nblocks = Nev_tot/Nev_block;

  for(int iq = 0; iq < Nq; iq++){
    for(int iev = 1; iev <= Nev_block; iev++){
      for(int iblock_an = 0; iblock_an < Nblocks; iblock_an++){
	
 			  //array_sort[iq*(Nev_tot+1)+(iev+Nev_block*iblock_an)] = array_in[iq*(Nev_tot+1)+(iev+Nev_block*conv_array[iblock_an])];

 			  array_sort[iq*(Nev_tot+1)+(iev+Nev_block*conv_array[iblock_an])] = array_in[iq*(Nev_tot+1)+(iev+Nev_block*iblock_an)];
	
      }// iblock_an
    }// iev
  } // iq
  
}// sort_array


void sort_matrix( double* matrix, int Nrows, int Ncol, double* arr_ord_row, double* arr_ord_col){

  double *temp=(double*)malloc(sizeof(double)), *temp_arr=(double*)malloc(sizeof(double));  
	double *temp_array_rows=(double*)malloc(sizeof(double)*(Nrows)), *temp_array_col=(double*)malloc(sizeof(double)*(Ncol));


	for(int i = 0; i < Nrows; i++){
		temp_array_rows[i] = arr_ord_row[i];
  }// i

	for(int j = 0; j < Ncol; j++){
		temp_array_col[j] = arr_ord_col[j];
  }// j


	//// RIGHE
	for(int j = 0; j < Ncol; j++){  // scorro su tutte le colonne

		for(int k = 0; k < Nrows-1; k++){
			for(int i = Nrows-1; i > k; i--){
	  
				if( arr_ord_row[i-1] > arr_ord_row[i] ){

					*temp_arr = arr_ord_row[i-1];
					arr_ord_row[i-1] = arr_ord_row[i];
					arr_ord_row[i] = *temp_arr;

					*temp = matrix[(i-1)*Ncol+j];
					matrix[(i-1)*Ncol+j] = matrix[(i)*Ncol+j];
					matrix[(i)*Ncol+j] = *temp;

				}// if
			}// i
		}// k

  	for(int i = 0; i < Nrows; i++){
			arr_ord_row[i] = temp_array_rows[i];
		}// i

	}// j

  //// COLONNE
	for(int i = 0; i < Nrows; i++){  // scorro su tutte le righe
      
		for(int k = 0; k < Ncol-1; k++){
			for(int j = Ncol-1; j > k; j--){
	  
				if( arr_ord_col[j-1] > arr_ord_col[j] ){

					*temp_arr = arr_ord_col[j-1];
					arr_ord_col[j-1] = arr_ord_col[j];
					arr_ord_col[j] = *temp_arr;
	    
					*temp = matrix[i*Ncol+(j-1)];
					matrix[i*Ncol+(j-1)] = matrix[i*Ncol+(j)];
					matrix[i*Ncol+(j)] = *temp;
	    
				}// if
			}// j
		}// k

  	for(int j = 0; j < Ncol; j++){
			arr_ord_col[j] = temp_array_col[j];
		}// j

	}// i

  free(temp);
  free(temp_arr);
  free(temp_array_rows);
  free(temp_array_col);

}// sort_matrix


void sort_vector( double* vector, int N, double* arr_ord){

  double *temp=(double*)malloc(sizeof(double)), *temp_arr=(double*)malloc(sizeof(double));  
	double *temp_array=(double*)malloc(sizeof(double)*(N));

	for(int i = 0; i < N; i++){
		temp_array[i] = arr_ord[i];
  }// i

	for(int k = 0; k < N-1; k++){
		for(int i = N-1; i > k; i--){
	  
			if( arr_ord[i-1] > arr_ord[i] ){

				*temp_arr = arr_ord[i-1];
				arr_ord[i-1] = arr_ord[i];
				arr_ord[i] = *temp_arr;

				*temp = vector[i-1];
				vector[i-1] = vector[i];
				vector[i] = *temp;

			}// if
		}// i
	}// k

	for(int i = 0; i < N; i++){
		arr_ord[i] = temp_array[i];
	}// i

  free(temp);
  free(temp_arr);
  free(temp_array);

}// sort_vector























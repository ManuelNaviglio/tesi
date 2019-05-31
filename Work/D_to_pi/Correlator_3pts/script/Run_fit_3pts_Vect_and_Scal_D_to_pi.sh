#!/bin/bash

cd ../ # torno alla cartella principale

set -e

PATH=$PATH:/home/francesco/QCD/LAVORI/MANUEL/MyExe/

cp  Input_corr_3pts/null_file.out Input_corr_3pts/file_input_corr_3pts.out   # Copia preventiva di null_file.out su quello di input

Nkl=(2 3 3 0 2)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@       FIT DELLE TRE PUNTI            @@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#

#@@@ COMPILAZIONE fit_corr_3pts_Vect_and_Scal_D_to_pi.C  @@@#
optimization=3

dir_sourc=$(PathScript MySource)
header_path=$(PathScript MyInclude)

func1="${dir_sourc}/statitisical_analysis_functions.C"
func2="${dir_sourc}/fit_2pts_func.C"
func3="${dir_sourc}/clean_string_cluster.C"
func4="${dir_sourc}/fit_3pts_func.C"
func5="${dir_sourc}/const_fit_scorr.C"

g++ -O$optimization -o fit_corr_3pts_Vect_and_Scal_D_to_pi fit_corr_3pts_Vect_and_Scal_D_to_pi.C $func1 $func2 $func3 $func4 $func5 $(rootlib) -I$header_path

#@@@ CICL0 DI ESECUZIONE  @@@#
for ibeta in `seq 0 4`
do

  for imusea in `seq 0 ${Nkl[ibeta]}`
  do
    
	  for mu1 in 0
	  do
	    
	    for mu2 in `seq 4 6`
	    do
		
		    for th1 in `seq 0 3`
		    do
		    
		      for th2 in `seq 0 6`
		      do
			
			      for sme in 1
			      do
			    	    				
				      echo $th2 $th1 $mu1 $mu2 $imusea $ibeta
				      echo $ibeta $imusea $mu1 $mu2 $th1 $th2 $sme >> Input_corr_3pts/file_input_corr_3pts.out
				
				      ./fit_corr_3pts_Vect_and_Scal_D_to_pi 
				
				      cp  Input_corr_3pts/null_file.out  Input_corr_3pts/file_input_corr_3pts.out
			    
			      done  # sme
		      done  # th2
		    done  # th1
	    done  # mu2
	  done  # mu1
  done # imusea  
done #ibeta

#@@@ CANCELLO ESEGUIBILE  @@@#
#rm fit_corr_3pts_Vect_and_Scal_D_to_pi

#@@@ CANCELLO IL LOG  @@@#
#rm log_corr_3pts


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@       PLOT FATTORI DI FORMA          @@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#

cp  Input_corr_3pts/null_file.out Input_corr_3pts/file_input_corr_3pts.out   # Copia preventiva di null_file.out su quello di input

#@@@ COMPILAZIONE Vect_and_Scal_form_factor_vs_q2_D_to_pi.C  @@@#
optimization=3

dir_sourc=$(PathScript MySource)
header_path=$(PathScript MyInclude)

func1="${dir_sourc}/statitisical_analysis_functions.C"
func2="${dir_sourc}/fit_3pts_func.C"
func3="${dir_sourc}/clean_string_cluster.C"
func4="${dir_sourc}/fit_2pts_func.C"

g++ -O$optimization -o Vect_and_Scal_form_factor_vs_q2_D_to_pi Vect_and_Scal_form_factor_vs_q2_D_to_pi.C $func1 $func2 $func3 $func4 $(rootlib) -I$header_path

#@@@ CICL0 DI ESECUZIONE  @@@#
for ibeta in `seq 0 4`
do

  for imusea in `seq 0 ${Nkl[ibeta]}`
  do
    
	  for mu1 in 0
	  do
	    
	    for mu2 in `seq 4 6`
	    do
		
		    echo $mu1 $mu2 $imusea $ibeta
		    echo $ibeta $imusea $mu1 $mu2 >> Input_corr_3pts/file_input_corr_3pts.out
		
		    ./Vect_and_Scal_form_factor_vs_q2_D_to_pi
		
		    cp  Input_corr_3pts/null_file.out  Input_corr_3pts/file_input_corr_3pts.out
		
	    done  # mu2
	  done  # mu1
  done # imusea  
done # ibeta

#@@@ CANCELLO ESEGUIBILE  @@@#
rm Vect_and_Scal_form_factor_vs_q2_D_to_pi


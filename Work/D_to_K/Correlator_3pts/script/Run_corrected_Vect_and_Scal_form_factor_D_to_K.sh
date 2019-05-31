#!/bin/bash

set -e

PATH=$PATH:/home/francesco/QCD/LAVORI/MANUEL/MyExe/

cd ../ # torno alla cartella principale

cp  Input_corr_3pts/null_file.out Input_corr_3pts/file_input_corr_3pts.out   # Copia preventiva di null_file.out su quello di input

Nkl[0]=3
Nkl[1]=4
Nkl[2]=4
Nkl[3]=1
Nkl[4]=3

#@@@ SCEGLIERE SE APPLICARE O MENO LA CORREZIONE  @@@#

no_correction=1  # SCEGLIERE 0 PER APPLICARE LA CORREZIONE, 1 PER NON APPLICARLA

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@   corrected_Vect_and_Scal_form_factor_800   @@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#

#@@@ COMPILAZIONE corrected_Vect_and_Scal_form_factor_800_D_to_K.C  @@@#
optimization=3

dir_sourc=$(PathScript MySource)
header_path=$(PathScript MyInclude)

func1="${dir_sourc}/statitisical_analysis_functions.C"
func2="${dir_sourc}/lettura_ultimate_input.C"

g++ -O$optimization -g -o corrected_Vect_and_Scal_form_factor_800_D_to_K corrected_Vect_and_Scal_form_factor_800_D_to_K.C $func1 $func2 $(rootlib) -I$header_path


#@@@ CICL0 DI ESECUZIONE  @@@#
rm -fr log.out

for i in `seq 0 4`
do
    
    musea=0
    
    while [ $musea -lt $((Nkl[i])) ];
    do
	
			for th1 in `seq 0 3`
			do
	    
				for th2 in `seq 0 6`
				do
					
					echo $th2 $th1 $musea $i
					echo $i $musea $th1 $th2 $no_correction >> Input_corr_3pts/file_input_corr_3pts.out
		
					./corrected_Vect_and_Scal_form_factor_800_D_to_K
					
		
					cp  Input_corr_3pts/null_file.out  Input_corr_3pts/file_input_corr_3pts.out
		
				done  # th2
			done  # th1
	
		musea=$((musea+1))
		done # musea
    
done #i

#@@@ CANCELLO ESEGUIBILE  @@@#
#rm corrected_Vect_and_Scal_form_factor_800_D_to_K

#@@@ CANCELLO IL FILE LOG  @@@#
#rm log.out

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@   corrected_Vect_and_Scal_form_factor_800   @@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#

#@@@ COMPILAZIONE Vect_and_Scal_form_factor_vs_q2_correction_800_D_to_K.C  @@@#
optimization=3

dir_sourc=$(PathScript MySource)
header_path=$(PathScript MyInclude)

func1="${dir_sourc}/statitisical_analysis_functions.C"
func2="${dir_sourc}/lettura_ultimate_input.C"

g++ -O$optimization -o Vect_and_Scal_form_factor_vs_q2_correction_800_D_to_K Vect_and_Scal_form_factor_vs_q2_correction_800_D_to_K.C $func1 $func2 $(rootlib) -I$header_path

#@@@ CICL0 DI ESECUZIONE  @@@#
rm -fr log.out
for i in `seq 0 4`
do
    
    musea=0
    
    while [ $musea -lt $((Nkl[i])) ];
    do
	
			echo $musea $i
			echo $i $musea $no_correction >> Input_corr_3pts/file_input_corr_3pts.out
	
				./Vect_and_Scal_form_factor_vs_q2_correction_800_D_to_K |tee -a  log.out
	
				cp  Input_corr_3pts/null_file.out  Input_corr_3pts/file_input_corr_3pts.out
	
		musea=$((musea+1))
    done # musea
    
done #i

#@@@ CANCELLO L'ESEGUIBILE  @@@#
#rm Vect_and_Scal_form_factor_vs_q2_correction_800_D_to_K

#@@@ CANCELLO IL FILE LOG  @@@#
#rm log.out

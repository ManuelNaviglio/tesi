#!/bin/bash

cd ../ # torno alla cartella principale

PATH=$PATH:/home/francesco/QCD/LAVORI/MANUEL/MyExe/

cp  Input_corr_2pts/null_file.out Input_corr_2pts/file_input_corr_2pts.out   # Copia preventiva di null_file.out su quello di input

#cp Input_corr_2pts/null_file.out OUTPUT_SMEAR/M/LOG/sintetic_log.out   # Copia preventiva di null_file.out su quello sintetic_log perché altrimenti mette tutte le informazioni in append sotto quelle del lancio precedente

#cp Input_corr_2pts/null_file.out OUTPUT_SMEAR/sinh_M/LOG/sintetic_log.out   # Copia preventiva di null_file.out su quello sintetic_log perché altrimenti mette tutte le informazioni in append sotto quelle del lancio precedente

Nkl[0]=3
Nkl[1]=4
Nkl[2]=4
Nkl[3]=1
Nkl[4]=3

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@       FIT DELLE DUE PUNTI            @@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#

#@@@ COMPILAZIONE fit_corr_2pts.C  @@@#
optimization=3

dir_sourc=$(PathScript MySource)
header_path=$(PathScript MyInclude)

func1="${dir_sourc}/statitisical_analysis_functions.C"
func2="${dir_sourc}/fit_2pts_func.C"
func3="${dir_sourc}/clean_string_cluster.C"
func4="${dir_sourc}/const_fit_scorr.C"


g++ -O$optimization -o fit_corr_2pts fit_corr_2pts.C $func1 $func2 $func3 $func4 $(rootlib) -I$header_path

#@@@ CICL0 DI ESECUZIONE  @@@#
for i in `seq 0 4`
do
    
    musea=0
    
    while [ $musea -lt $((Nkl[i])) ];
    do
	
	for mu1 in 0
	do
	    
	    for mu2 in `seq 0 6`
	    do
		    
		for th1 in `seq 0 3`
		do
		    
		    if [ $th1 -eq 0 ] ; then th2=6 ; fi
		    if [ $th1 -eq 1 ] ; then th2=5 ; fi
		    if [ $th1 -eq 2 ] ; then th2=4 ; fi
		    if [ $th1 -eq 3 ] ; then th2=3 ; fi
		    
		    for sme in 1
		    do
			
			echo $th2 $th1 $mu1 $mu2 $musea $i
			echo $i $musea $mu1 $mu2 $th1 $th2 $sme >> Input_corr_2pts/file_input_corr_2pts.out
			
			./fit_corr_2pts > log_corr_2pts
			
			cp  Input_corr_2pts/null_file.out  Input_corr_2pts/file_input_corr_2pts.out
			
		    done  # sme
		done  # th1
		
	    done  # mu2
	done  # mu1
	
    musea=$((musea+1))
    done # musea
    
done #i

#@@@ CANCELLO ESEGUIBILE  @@@#
rm fit_corr_2pts

#@@@ CANCELLO IL LOG  @@@#
rm log_corr_2pts

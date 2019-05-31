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

MASS=1
MATRIX_EL_V_and_S=1

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@  SPLINE ELEMENTI DI MATRICE V_and_S  @@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#

if [ $MATRIX_EL_V_and_S -eq 1 ];
then

    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
    #@@@@@  SPLINE_V_S_DK_charm  @@@@@#
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
    
    #@@@ COMPILAZIONE spline_Vect_and_Scal.C  @@@#
    optimization="3 -g"

    dir_sourc=$(PathScript MySource)
    header_path=$(PathScript MyInclude)
    
    func1="${dir_sourc}/statitisical_analysis_functions.C"
    func2="${dir_sourc}/lettura_ultimate_input.C"
    func3="${dir_sourc}/spline_func.C"

    declared="-D SPLINE_V_S_DK -D SPLINE_V_S_DK_charm"
    
    g++ -O$optimization -o spline_Vect_and_Scal spline_Vect_and_Scal.C $func1 $func2 $func3 $(rootlib) -I$header_path $declared
    
    #@@@ CICL0 DI ESECUZIONE  @@@#
    for i in `seq 0 4`
    do
	
	musea=0
	
	while [ $musea -lt $((Nkl[i])) ];
	do
	    
	    for im in `seq 1 3`
	    do
		
		for th1 in `seq 0 3`
		do
		    
		    for th2 in `seq 0 6`
		    do
			
			echo $i $musea $im $th1 $th2
			echo $i $musea $im $th1 $th2 >> Input_corr_3pts/file_input_corr_3pts.out
			
			./spline_Vect_and_Scal
			
			cp  Input_corr_3pts/null_file.out  Input_corr_3pts/file_input_corr_3pts.out
			
		    done  # th2
		done  # th1
	    done # im
	    
	musea=$((musea+1))
	done # musea
	
    done #i
    
    #@@@ CANCELLO L'ESEGUIBILE  @@@#
    #rm spline_Vect_and_Scal
    
    #@@@ CANCELLO IL FILE LOG  @@@#
    #rm log_spline_V_and_S



    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
    #@@@@@  SPLINE_V_S_DK_strange  @@@@@#
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
    
    #@@@ COMPILAZIONE spline_Vect_and_Scal.C  @@@#
    optimization="3 -g"

    dir_sourc=$(PathScript MySource)
    header_path=$(PathScript MyInclude)

    func1="${dir_sourc}/statitisical_analysis_functions.C"
    func2="${dir_sourc}/lettura_ultimate_input.C"
    func3="${dir_sourc}/spline_func.C"

    declared="-D SPLINE_V_S_DK -D SPLINE_V_S_DK_strange"
    
    g++ -O$optimization -o spline_Vect_and_Scal spline_Vect_and_Scal.C $func1 $func2 $func3 $(rootlib) -I$header_path $declared
    
    #@@@ CICL0 DI ESECUZIONE  @@@#
    for i in `seq 0 4`
    do
	
	musea=0
	
	while [ $musea -lt $((Nkl[i])) ];
	do
	    
	    for im in 0 # questo non viene usato in questo caso (è l'indice del quark leggero nel caso spline_V_and_S_DK_charm)
	    do
		
		for th1 in `seq 0 3`
		do
		    
		    for th2 in `seq 0 6`
		    do
			
			echo $i $musea $im $th1 $th2
			echo $i $musea $im $th1 $th2 >> Input_corr_3pts/file_input_corr_3pts.out
			
			./spline_Vect_and_Scal 
			
			cp  Input_corr_3pts/null_file.out  Input_corr_3pts/file_input_corr_3pts.out
			
		    done  # th2
		done  # th1
	    done # im
	    
	musea=$((musea+1))
	done # musea
	
    done #i
    
    #@@@ CANCELLO L'ESEGUIBILE  @@@#
    #rm spline_Vect_and_Scal
    
    #@@@ CANCELLO IL FILE LOG  @@@#
    #rm log_spline_V_and_S
    
fi

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@    SPLINE MASSE     @@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#

if [ $MASS -eq 1 ];
then
    
    #@@@ COMPILAZIONE spline_Vect_and_Scal.C  @@@#
    optimization="3 -g"

    dir_sourc=$(PathScript MySource)
    header_path=$(PathScript MyInclude)

    func1="${dir_sourc}/statitisical_analysis_functions.C"
    func2="${dir_sourc}/lettura_ultimate_input.C"
    func3="${dir_sourc}/spline_func.C"

    declared="-D SPLINE_MD -D SPLINE_MPI -D SPLINE_MK"
    
    g++ -O$optimization -o spline_Vect_and_Scal spline_Vect_and_Scal.C $func1 $func2 $func3 $(rootlib) -I$header_path $declared
    
    #@@@ CICL0 DI ESECUZIONE  @@@#
    for i in `seq 0 4`
    do
	
	musea=0
	
	while [ $musea -lt $((Nkl[i])) ];
	do
	    
	    for im in 0 # questo non viene usato in questo caso (è l'indice del quark leggero nel caso spline_V_and_S_DK_charm)
	    do
		
		for th1 in 3            # fissati perché sto facendo la spline delle masse
		do
		    
		    for th2 in 3        # fissati perché sto facendo la spline delle masse
		    do
			
			echo $i $musea $im $th1 $th2
			echo $i $musea $im $th1 $th2 >> Input_corr_3pts/file_input_corr_3pts.out
			
			./spline_Vect_and_Scal 
			
			cp  Input_corr_3pts/null_file.out  Input_corr_3pts/file_input_corr_3pts.out
			
		    done  # th2
		done  # th1
	    done # im
	    
	musea=$((musea+1))
	done # musea
	
    done #i
    
    #@@@ CANCELLO L'ESEGUIBILE  @@@#
    #rm spline_Vect_and_Scal
    
    #@@@ CANCELLO IL FILE LOG  @@@#
    #rm log_spline_V_and_S
    
fi

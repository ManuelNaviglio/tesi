#!/bin/bash

cd ../ # torno alla cartella principale

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@     HFAG     @@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#

#@@@ COMPILAZIONE confronto_HFAG_D_to_K.C  @@@#
optimization=3

dir_sourc=$(PathScript MySource)
header_path=$(PathScript MyInclude)

func1="${dir_sourc}/statitisical_analysis_functions.C"
func2="${dir_sourc}/lettura_ultimate_input.C"
func3="${dir_sourc}/random_gen.C"
func4="${dir_sourc}/algebra_func.C"


g++ -O$optimization -o confronto_HFAG_D_to_K confronto_HFAG_D_to_K.C $func1 $func2 $func3 $func4 $(rootlib) -I$header_path

./confronto_HFAG_D_to_K > log_confronto_HFAG_D_to_K

#@@@ CANCELLO ESEGUIBILE  @@@#
rm confronto_HFAG_D_to_K

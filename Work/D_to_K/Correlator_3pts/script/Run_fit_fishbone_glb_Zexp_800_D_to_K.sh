#!/bin/bash

cd ../ # torno alla cartella principale

set -e

PATH=$PATH:/home/francesco/QCD/LAVORI/MANUEL/MyExe/

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@             FIT FISHBONE             @@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#

#@@@ COMPILAZIONE fit_fishbone_glb_Zexp_800_D_to_K.C  @@@#
optimization="3 -g"

dir_sourc=$(PathScript MySource)
header_path=$(PathScript MyInclude)

func1="${dir_sourc}/statitisical_analysis_functions.C"
func2="${dir_sourc}/lettura_ultimate_input.C"

declared="-D STD"
#declared="-D CUT"
#declared="-D A4"
#declared="-D A4_CUT"

#declared="-D NO_HYP"

g++ -O$optimization -o fit_fishbone_glb_Zexp_800_D_to_K fit_fishbone_glb_Zexp_800_D_to_K.C $func1 $func2 $(rootlib) -I$header_path $declared

./fit_fishbone_glb_Zexp_800_D_to_K

#@@@ CANCELLO ESEGUIBILE  @@@#
#rm fit_fishbone_glb_Zexp_800_D_to_K


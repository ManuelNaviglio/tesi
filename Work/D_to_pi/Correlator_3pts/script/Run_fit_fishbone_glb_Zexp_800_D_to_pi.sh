#!/bin/bash

set -e

PATH=$PATH:/home/francesco/QCD/LAVORI/MANUEL/MyExe/

cd ../ # torno alla cartella principale

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@             FIT FISHBONE             @@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#

#@@@ COMPILAZIONE fit_fishbone_glb_Zexp_800_D_to_pi.C  @@@#
optimization="0 -g"

dir_sourc=$(PathScript MySource)
header_path=$(PathScript MyInclude)

func1="${dir_sourc}/statitisical_analysis_functions.C"
func2="${dir_sourc}/lettura_ultimate_input.C"
declared="-D STD"

g++ -O$optimization -o fit_fishbone_glb_Zexp_800_D_to_pi fit_fishbone_glb_Zexp_800_D_to_pi.C $func1 $func2 $(rootlib) -I$header_path $declared

./fit_fishbone_glb_Zexp_800_D_to_pi #> log_fishbone

#@@@ CANCELLO ESEGUIBILE  @@@#
#rm fit_fishbone_glb_Zexp_800_D_to_pi


#!/bin/bash

set -e

PATH=$PATH:/home/francesco/QCD/LAVORI/MANUEL/MyExe/

cd ../ # torno alla cartella principale

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@     COVARIANCE Vcd AND Vcs             @@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#

#@@@ COMPILAZIONE matrice_covarianza_Vcd_Vcs.C @@@#
optimization="3 -g"

dir_sourc=$(PathScript MySource)
header_path=$(PathScript MyInclude)

func1="${dir_sourc}/statitisical_analysis_functions.C"

g++ -O$optimization -o matrice_covarianza_Vcd_Vcs matrice_covarianza_Vcd_Vcs.C $func1 $(rootlib) -I$header_path

./matrice_covarianza_Vcd_Vcs 

#@@@ CANCELLO ESEGUIBILE  @@@#
#rm matrice_covarianza_Vcd_Vcs

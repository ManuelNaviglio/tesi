#!/bin/bash

cd ../ # torno alla cartella principale

set -e

PATH=$PATH:/home/francesco/QCD/LAVORI/MANUEL/MyExe/

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@            FIT SYSTEMATICS             @@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#

#@@@ DICHIARAZIONI @@@#

declared_Vcs_cov_lat="VCS_COV_LAT_TOT"
#declared_Vcs_cov_lat="VCS_COV_LAT_BLOCK"
#declared_Vcs_cov_lat="VCS_COV_LAT_DIAG"
#declared_Vcs_cov_lat="NO_VCS_COV_LAT"

declared_Belle="BELLE"
#declared_Belle="NO_BELLE"

#@@@ COSTRUISCO IL TEMPLATE PER IL FILE DI ESTRAZIONE DI Vcs @@@#

cat Estrazione_Vcs/estrazione_Vcs_template.xmg > Estrazione_Vcs/Estrazione_Vcs_${declared_Vcs_cov_lat}_${declared_Belle}

#@@@ COMPILAZIONE fit_syst_D_to_K.C  @@@#
optimization="3 -g"



dir_sourc=$(PathScript MySource)
header_path=$(PathScript MyInclude)

func1="${dir_sourc}/statitisical_analysis_functions.C"
func2="${dir_sourc}/lettura_ultimate_input.C"
func3="${dir_sourc}/random_gen.C"
func4="${dir_sourc}/algebra_func.C"
func5="${dir_sourc}/array_manipulation.C"
func6="${dir_sourc}/synt_data_func.C"
func7="${dir_sourc}/grace_tools.C"
func8="${dir_sourc}/fit_with_Minuit.C"


g++ -O$optimization -o fit_syst_D_to_K fit_syst_D_to_K.C $func1 $func2 $func3 $func4 $func5 $func6 $func7 $func8 $(rootlib) $(gsllib) -I$header_path -D $declared_Vcs_cov_lat -D $declared_Belle

./fit_syst_D_to_K > log_syst_D_to_K_${declared_Vcs_cov_lat}_${declared_Belle}

#@@@ AGGIUNGO I DATI AL FILE DI ESTRAZIONE DI Vcs @@@#

grep -A 170 '##CUT_STEP' log_syst_D_to_K_${declared_Vcs_cov_lat}_${declared_Belle} >> Estrazione_Vcs/Estrazione_Vcs_${declared_Vcs_cov_lat}_${declared_Belle}

#@@@ CANCELLO ESEGUIBILE  @@@#
#rm fit_syst_D_to_K

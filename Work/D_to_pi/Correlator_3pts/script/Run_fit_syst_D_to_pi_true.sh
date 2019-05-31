#!/bin/bash

cd ../ # torno alla cartella principale

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@            FIT SYSTEMATICS             @@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#

#@@@ DICHIARAZIONI @@@#

declared_Vcd_cov_lat="VCD_COV_LAT_TOT"
#declared_Vcd_cov_lat="VCD_COV_LAT_BLOCK"
#declared_Vcd_cov_lat="VCD_COV_LAT_DIAG"
#declared_Vcd_cov_lat="NO_VCD_COV_LAT"

declared_Belle="BELLE"
#declared_Belle="NO_BELLE"

#@@@ COSTRUISCO IL TEMPLATE PER IL FILE DI ESTRAZIONE DI Vcd @@@#

cat Estrazione_Vcd/estrazione_Vcd_template.xmg > Estrazione_Vcd/Estrazione_Vcd_${declared_Vcd_cov_lat}_${declared_Belle}

#@@@ COMPILAZIONE fit_syst_D_to_pi.C  @@@#
optimization=3

dir_sourc="/home/giorgio/NUOVI/VECTOR_AND_SCALAR/sorgenti"

func1="${dir_sourc}/statitisical_analysis_functions.C"
func2="${dir_sourc}/lettura_ultimate_input.C"
func3="${dir_sourc}/random_gen.C"
func4="${dir_sourc}/algebra_func.C"
func5="${dir_sourc}/array_manipulation.C"
func6="${dir_sourc}/synt_data_func.C"
func7="${dir_sourc}/grace_tools.C"
func8="${dir_sourc}/fit_with_Minuit.C"
header_path="/home/giorgio/NUOVI/VECTOR_AND_SCALAR/header/"

g++ -O$optimization -o fit_syst_D_to_pi fit_syst_D_to_pi.C $func1 $func2 $func3 $func4 $func5 $func6 $func7 $func8 $(rootlib) -I$header_path -D $declared_Vcd_cov_lat -D $declared_Belle

./fit_syst_D_to_pi > log_syst_D_to_pi_${declared_Vcd_cov_lat}_${declared_Belle}

#@@@ AGGIUNGO I DATI AL FILE DI ESTRAZIONE DI Vcd @@@#

grep -A 128 '##CUT_STEP' log_syst_D_to_pi_${declared_Vcd_cov_lat}_${declared_Belle} >> Estrazione_Vcd/Estrazione_Vcd_${declared_Vcd_cov_lat}_${declared_Belle}

#@@@ CANCELLO ESEGUIBILE  @@@#
rm fit_syst_D_to_pi

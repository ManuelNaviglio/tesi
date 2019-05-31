#!/bin/bash

cd ../ # torno alla cartella principale

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@            FIT SYSTEMATICS             @@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#

#@@@ DICHIARAZIONI @@@#

#declared_Vcd_cov_lat="VCD_COV_LAT_TOT"
declared_Vcd_cov_lat="VCD_COV_LAT_BLOCK"
#declared_Vcd_cov_lat="VCD_COV_LAT_DIAG"
#declared_Vcd_cov_lat="NO_VCD_COV_LAT"

declared_Belle="BELLE"
#declared_Belle="NO_BELLE"

#@@@ COMPILAZIONE fit_syst_D_to_pi.C  @@@#
optimization=3

dir_sourc="/home/giorgio/MyLib/source"
header_path="/home/giorgio/MyLib/include"

func1="${dir_sourc}/statitisical_analysis_functions.C"
func2="${dir_sourc}/lettura_ultimate_input.C"
func3="${dir_sourc}/random_gen.C"
func4="${dir_sourc}/algebra_func.C"
func5="${dir_sourc}/array_manipulation.C"
func6="${dir_sourc}/synt_data_func.C"
func7="${dir_sourc}/grace_tools.C"
func8="${dir_sourc}/fit_with_Minuit.C"

g++ -O$optimization -o fit_syst_D_to_pi fit_syst_D_to_pi.C $func1 $func2 $func3 $func4 $func5 $func6 $func7 $func8 $(rootlib) -lgsl -lgslcblas -I$header_path -D $declared_Vcd_cov_lat -D $declared_Belle

scp fit_syst_D_to_pi giorgio7@entree:/lpti/giorgio7

#@@@ CANCELLO ESEGUIBILE  @@@#
rm fit_syst_D_to_pi

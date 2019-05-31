#!/bin/bash
for i in \
    Run_fit_3pts_Vect_and_Scal_D_to_K.sh \
	Run_spline_Vect_and_Scal_D_to_K.sh \
	Run_fit_fishbone_glb_Zexp_800_D_to_K.sh \
	Run_corrected_Vect_and_Scal_form_factor_D_to_K.sh
do
    bash $i
done
    

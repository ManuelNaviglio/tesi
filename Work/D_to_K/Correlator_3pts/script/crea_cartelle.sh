for ens in 1.90/24/0.0040/ 1.90/24/0.0060/ 1.90/24/0.0080/ 1.90/24/0.0100/ 1.90/32/0.0030/ 1.90/32/0.0040/ 1.90/32/0.0050/ 1.95/24/0.0085/ 1.95/32/0.0025/ 1.95/32/0.0035/ 1.95/32/0.0055/ 1.95/32/0.0075/ 2.10/48/0.0015/ 2.10/48/0.0020/ 2.10/48/0.0030/
do
    for i in {hyp_corrections,.}/{f0_S,fplus,fzero,fminus,Plateaux,S0,V0,Vi}/{Heavy_to_light,Hl_sym_av_S_ratio,Hl_sym,Hl_sym_av,light_to_Heavy} \
				{fzero,fplus,fminus}_corrected/{Heavy_to_light,Hl_sym_av_S_ratio,Hl_sym,Hl_sym_av,light_to_Heavy} \
				fzero_corrected/Hl_sym_av/fzero_vs_q2 \
				fplus_corrected/Hl_sym_av/fplus_vs_q2
    do
	mkdir -pv ../OUTPUT_SMEAR/with_S/E_from_stdDR/$i/$ens
    done
done

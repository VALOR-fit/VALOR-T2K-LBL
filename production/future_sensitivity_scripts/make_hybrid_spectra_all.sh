allocation=( FHC FHC_RHC )
POT=( 7.8E21 10E21 15E21)

for p in "${POT[@]}"; do
    for a in "${allocation[@]}"; do
        . ${T2KOSCANA}/production/future_sensitivity_scripts/make_data_from_spectra.sh spectra_${a}_${p}_excess_ data_${a}_${p}_excess.root
	hadd -f data_${a}_${p}_hybrid.root ${T2KOSCANA}/inputs/official/2018/data/SK_data_Run19d_valor.root data_${a}_${p}_excess.root
    done
done

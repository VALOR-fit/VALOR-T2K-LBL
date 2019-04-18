allocation=( FHC FHC_RHC )
POT=( 7.8E21 10E21 15E21)

for p in "${POT[@]}"; do
    for a in "${allocation[@]}"; do
        . ${T2KOSCANA}/production/future_sensitivity_scripts/make_data_from_spectra.sh spectra_${a}_${p}_total_ data_${a}_${p}_pure_AsimovA.root
    done
done

ok=1
disable_syst="false"


run="$1"
amount="$2"
config_type="$3"
config="analysis_configuration_oa2018_run1-9d_v1_postfit_future_${amount}_${run}_${config_type}_hybrid.xml"

if [ "$run" == "" ] || [ "$amount" == "" ] || [ "$config_type" == "" ]; then
    ok=0
fi
   
mixing12=0.307
mixing13=0.0212
mixing23=0.528
dmsq12=7.53E-5
dmsq23=2.509E-3
dcp=-1.601
octant="two"

if [ "$ok" == "1" ]; then
   sample=( 1Re 1Rmu MultiRe RHC1Re RHC1Rmu )
   for s in "${sample[@]}"; do
      if [[ "$s" == "1Rmu" || "$s" == "RHC1Rmu" ]]; then
         erec_high=30.0
      else
         erec_high=1.25
      fi
      command="${T2KOSCANA}/macros/spectra/test_sk_ereco_spectrum_nominal_mainmode_joint_2D.C"
      std_args="\"spectra_${run}_${amount}_${config_type}\",\"${s}\",\"3f_osc\",$mixing12,$mixing13,$mixing23,$dcp,$dmsq12,$dmsq23,-99,${disable_syst},\"${config}\",true,false,\"sinsqtheta\",\"normal\",\"auto\",\"$octant\",false,0.0,$erec_high,0,\"\",false"
      root2k -q -b "${command}(${std_args})"
      #echo "root2k -q -b ${command}(${std_args})"
   done
else
   echo -e "Usage: make_spectra.sh <run> <amount>"
   echo -e "   run: FHC (for extra POT = 100% FHC) | FHC_RHC (for extra POT = 50/50 FHC/RHC)"
   echo -e "   amount: 7.8E21 | 10E21 | 15E21"
   echo -e "   config_type: total (for config containing total future POT) | excess (for config containing future POT - current POT)"
fi


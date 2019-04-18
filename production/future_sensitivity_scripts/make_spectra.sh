ok=1
disable_syst="false"

if [[ "$1" == "fhc" ]]; then
   run="fhc"
   if [[ "$2" == "full" ]]; then
      amount="full"
      config="analysis_configuration_pac_201901_fhc_postfit.xml"
   elif [[ "$2" == "future" ]]; then
      amount="future"
      config="analysis_configuration_pac_201901_fhc_future_postfit.xml"
   elif [[ "$2" == "standard" ]]; then
      amount="standard"
      config="analysis_configuration_pac_201901_standard_postfit.xml"
   else
      ok=0
   fi
elif [[ "$1" == "rhc" ]]; then
   run="rhc"
   if [[ "$2" == "full" ]]; then
      amount="full"
      config="analysis_configuration_pac_201901_rhc_postfit.xml"
   elif [[ "$2" == "future" ]]; then
      amount="future"
      config="analysis_configuration_pac_201901_rhc_future_postfit.xml"
   elif [[ "$2" == "standard" ]]; then
      amount="standard"
      config="analysis_configuration_pac_201901_standard_postfit.xml"
   else
      ok=0
   fi
else
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
      std_args="\"spectra_${run}_${amount}\",\"${s}\",\"3f_osc\",$mixing12,$mixing13,$mixing23,$dcp,$dmsq12,$dmsq23,-99,${disable_syst},\"${config}\",true,false,\"sinsqtheta\",\"normal\",\"auto\",\"$octant\",false,0.0,$erec_high,0,\"\",false"
      bsub -q express -G t2kgrp -R "rusage[mem=3072,swap=3280]" -oo make_spectra_${run}_${amount}_${s}.log ${T2KOSCANA}/production/warwick_submit_command.sh root2k -q -b "${command}(${std_args})"
   done
else
   echo -e "Usage: make_spectra.sh <run> <amount>"
   echo -e "   run: fhc | rhc"
   echo -e "   amount: full | future"
fi


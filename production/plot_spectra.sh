ok=1
if [[ "$1" == "predicted" || "$1" == "observed" ]]; then
   type=$1
else
   ok=0
fi

if [[ "$2" != "" ]]; then
   osc=$2
else
   ok=0
fi

if [ "$ok" == "1" ]; then
   sample=(1Re 1Rmu MultiRe RHC1Re RHC1Rmu)
   command=${T2KOSCANA}/macros/tn_plotting_stuffs/plot_spectra.py
   for s in ${sample[@]}; do
      if [ "$type" == "observed" ]; then
         if [ "${osc}" == "asimov2017" ]; then
            data="-d $T2KOSCANA/inputs/official/2018/data/SK_data_Run18_valor.root"
         elif [ "${osc}" == "run9asimov2018" ]; then
            data="-d $T2KOSCANA/inputs/official/2018/data/SK_dataRun9only_valor.root"
         elif [ "${osc}" == "asimov_run1-9" ]; then
            data="-d $T2KOSCANA/inputs/official/2018/data/SK_data_Run19d_valor.root"
            data_old="-x $T2KOSCANA/inputs/official/2018/data/SK_data_Run19c_valor.root"
            data_new="-d $T2KOSCANA/inputs/official/2018/data/diff_9d.root"
         fi
         output="-o spectrum_${osc}_${s}_data"
      else
         data=""
         output="-o spectrum_${osc}_${s}"
      fi
      python ${command} -f spectra_${osc}_${s}.root ${output} ${data} -s ${s}
#      python ${command} -f spectra_${osc}_${s}.root ${output} ${data_new} ${data_old} -s ${s}
   done
else
   echo -e "Usage: plot_spectra,sh <type> <osc>"
   echo -e "   type: predicted | observed"
   echo -e "   osc: asimova | asimovb | unosc"
fi


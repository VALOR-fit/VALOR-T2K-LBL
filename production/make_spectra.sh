ok=1
if [[ "$1" != "" ]]; then
   osc=$1
else
   ok=0
fi

if [[ "$2" == "true" || "$2" == "false" ]]; then
   disable_syst="$2"
else
   ok=0
fi

if [ "$3" == "prefit" ]; then
   config="analysis_configuration_hkfluxsophie_10year_prefit.xml" #"analysis_configuration_oa2018_run1-9d_v1_prefit.xml"
   banff="_prefit"
else
   config="analysis_configuration_hkfluxsophie_10year_postfit.xml" #"analysis_configuration_oa2018_run1-9d_v1_postfit.xml"
   banff=""
fi

if [ "${osc}" == "unosc" ]; then
   mixing12=0.0
   mixing13=0.0
   mixing23=0.0
   dmsq12=7.53E-5
   dmsq23=2.509E-3
   dcp=0
   octant="two"
elif [ "${osc}" == "hkdr" ]; then
   mixing12=0.32
   mixing13=0.0256584
   mixing23=0.5
   dmsq12=7.6E-5
   dmsq23=2.4E-3
   dcp=0
   octant="two"
elif [ "${osc}" == "asimova" ]; then
   mixing12=0.304
   mixing13=0.0219
   mixing23=0.528
   dmsq12=7.53E-5
   dmsq23=2.509E-3
   dcp=-1.601
   octant="two"
elif [ "${osc}" == "asimovb" ]; then
   mixing12=0.304
   mixing13=0.0219
   mixing23=0.450
   dmsq12=7.53E-5
   dmsq23=2.509E-3
   dcp=0.0
   octant="one"
elif [ "${osc}" == "asimov2018" ]; then
   mixing12=0.304
   mixing13=0.02235
   mixing23=0.535
   dmsq12=7.53E-5
   dmsq23=2.440E-3
   dcp=-1.75929
   octant="two"
elif [ "${osc}" == "asimov_run1-9" ]; then
   mixing12=0.307
   mixing13=0.0215625
   mixing23=0.535               # From unsmeared fit
   dmsq12=7.53E-5
   dmsq23=2.452E-3              # From unsmeared fit
   dcp=-1.88496
   octant="two"
elif [ "${osc}" == "asimov2019" ]; then
   mixing12=0.304
   mixing13=0.02235
   mixing23=0.535
   dmsq12=7.53E-5
   dmsq23=2.440E-3
   dcp=-1.75929
   octant="two"
elif [ "${osc}" == "asimova_fake" ]; then
   mixing12=0.304
   mixing13=0.0219
   mixing23=0.528
   dmsq12=7.53E-5
   dmsq23=2.509E-3
   dcp=-1.601
   octant="two"
   config="fakedata/analysis_configuration_oa2018_run1-9c_asimov_postfit.xml"
elif [ "${osc}" == "asimovb_fake" ]; then
   mixing12=0.304
   mixing13=0.0219
   mixing23=0.450
   dmsq12=7.53E-5
   dmsq23=2.509E-3
   dcp=0.0
   octant="one"
   config="fakedata/analysis_configuration_oa2018_run1-9c_asimov_postfit.xml"
fi

if [ "$ok" == "1" ]; then
   sample=( 1Re 1Rmu MultiRe RHC1Re RHC1Rmu )
   for s in "${sample[@]}"; do
      if [[ "$s" == "1Rmu" || "$s" == "RHC1Rmu" ]]; then
         erec_high=30.0
      else
         erec_high=1.25
      fi
      command="${T2KOSCANA}/macros/spectra/test_sk_ereco_spectrum_nominal_mainmode_joint_2D.C"
      std_args="\"spectra_${osc}${banff}\",\"${s}\",\"3f_osc\",$mixing12,$mixing13,$mixing23,$dcp,$dmsq12,$dmsq23,-99,${disable_syst},\"${config}\",true,false,\"sinsqtheta\",\"normal\",\"auto\",\"$octant\",false,0.0,$erec_high,0,\"\",false"
      echo root2k -l -b -q $command'('$std_args')'
      root2k -l -b -q $command'('$std_args')'
      #bsub -q express -R "rusage[mem=3072,swap=3280]" -oo make_spectra_${osc}_${s}.log ${T2KOSCANA}/production/warwick_submit_command.sh root2k -q -b "${command}(${std_args})"
   done
else
   echo -e "Usage: make_spectra.sh <osc> <disable_systematics> [banff]"
   echo -e "   osc: asimova | asimovb | unosc | asimov2018 | asimova2019 | hkdr"
   echo -e "   disable_systematics: false | true"
   echo -e "   banff: prefit | postfit (default)"
fi


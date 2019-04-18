ok=1

# Defaults
mixing12="--mix12 0.307"
mixing13="--mix13 0.0212"
mixing23="--mix23 0.528"
dmsq12="--dmsq12 7.53E-5"
dmsq23="--dmsq23 2.509E-3"
deltacp="--deltacp -1.601"
octant="--mix-param-octant two"
osc_params="$mixing23 $dmsq23 $mixing13 $mixing12 $dmsq12 $deltacp $octant"

if [ "$1" == "asimov_fhc" ]; then
   datafit="--do-toys 0 --fakedata-root-file $T2KOSCANA/inputs/official/2018/data/SK_pac_fhc_asimov.root"
   config_file="--config-file analysis_configuration_pac_201901_fhc_postfit.xml"
elif [ "$1" == "asimov_rhc" ]; then
   datafit="--do-toys 0 --fakedata-root-file $T2KOSCANA/inputs/official/2018/data/SK_pac_rhc_asimov.root"
   config_file="--config-file analysis_configuration_pac_201901_rhc_postfit.xml"
elif [ "$1" == "asimov_standard" ]; then
   datafit="--do-toys 0 --fakedata-root-file $T2KOSCANA/inputs/official/2018/data/SK_pac_standard_asimov.root"
   config_file="--config-file analysis_configuration_pac_201901_standard_postfit.xml"
elif [ "$1" == "hybrid_fhc" ]; then
   datafit="--do-toys 0 --fakedata-root-file $T2KOSCANA/inputs/official/2018/data/SK_pac_fhc_hybrid.root"
   config_file="--config-file analysis_configuration_pac_201901_fhc_postfit.xml"
elif [ "$1" == "hybrid_rhc" ]; then
   datafit="--do-toys 0 --fakedata-root-file $T2KOSCANA/inputs/official/2018/data/SK_pac_rhc_hybrid.root"
   config_file="--config-file analysis_configuration_pac_201901_rhc_postfit.xml"
else
   ok=0
fi

if [ "$2" != "" ]
then
   dir="--file-margtoys /data/t2k/phrdqd/analysis/run1_9/toys/$2"
else
   ok=0
fi

if [[ "$3" == "normal" || "$3" == "inverted" ]]
then
   hierarchy="--hierarchy-marg $3"
else
   ok=0
fi

if [ "$4" == "reactor" ]; then
   penalty="--penalty-par0 1"
elif [ "$4" == "t2k" ]; then
   penalty="--penalty-par0 0"
else
   ok=0
fi

if [ "$5" == "" ]
then
   freepar0="--freepar0 deltacp"
else
   freepar0="--freepar0 $5"
fi

if [ "$6" == "" ]
then
   freepar1=""
else
   freepar1="--freepar1 $6"
fi

freepars="$freepar0 $freepar1"
constraint=$4

if [ "$ok" == "1" ]
then
   echo -e "Running fit - $freepars $hierarchy $constraint"
   resource="--resource warwick"
   $T2KOSCANA/production/submit_toy_fit_marginalization_exe $resource $datafit $hierarchy $penalty $osc_params $config_file $freepars $dir
else
   echo -e "Usage: fit.sh dataset toy_directory hierarchy constraint freepar0 [freepar1]"
fi

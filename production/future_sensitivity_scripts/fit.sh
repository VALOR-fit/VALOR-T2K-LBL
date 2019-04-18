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

study_type=$1
split=$2
POT=$3
marg_dir=$4
hier=$5
reac=$6
job_dir=$7
par0=$8
par1=$9

tag="future_sensitivity_Apr19_${hier}_${reac}_${split}_${POT}_${study_type}"
datafit="--do-toys 0 --fakedata-root-file /scratch/fbench/future_sensitivity_Apr2019/data/data_${split}_${POT}_${study_type}.root"
config_file="--config-file analysis_configuration_oa2018_run1-9d_v1_postfit_future_${POT}_${split}_total.xml"
if [ "$study_type" == "" ] || [ "$split" == "" ] || [ "$POT" == "" ]; then
    ok=0
fi

if [ "$marg_dir" != "" ]
then
   dir="--file-margtoys $marg_dir"
else
   ok=0
fi

if [[ "$hier" == "normal" || "$hier" == "inverted" ]]
then
   hierarchy="--hierarchy-marg $hier"
else
   ok=0
fi

if [ "$reac" == "reactor" ]; then
   penalty="--penalty-par0 1"
elif [ "$reac" == "t2k" ]; then
   penalty="--penalty-par0 0"
else
   ok=0
fi

if [ "$par0" == "" ]
then
   freepar0="--freepar0 deltacp"
else
   freepar0="--freepar0 $par0"
fi

if [ "$par1" == "" ]
then
   freepar1=""
else
   freepar1="--freepar1 $par1"
fi

freepars="$freepar0 $freepar1"
constraint=$reac

if [ "$ok" == "1" ]
then
   echo -e "Running fit - $tag"
   resource="--resource cedar"
   $T2KOSCANA/production/submit_toy_fit_marginalization_exe $resource $datafit $hierarchy $penalty $osc_params $config_file $freepars $dir --job-dir $job_dir --tag $tag --cpu-time 1:00:00
   echo
   echo
else
   echo -e "Usage: fit.sh study_type (pure_AsimovA | hybrid) split (FHC | FHC_RHC) POT (XE21) marg_toy_directory hierarchy constraint job_dir freepar0 [freepar1]"
fi

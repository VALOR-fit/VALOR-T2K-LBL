###################################################################################################
### Jobsub command to run multiple fits on the same node                                        ###
###                                                                                             ###
### To run submit the jobs, first setup your experiment                                         ###
### e.g source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh                      ###
### You might also need to "setup fife_utils"                                                   ###
###                                                                                             ###
### To submit the job run a command such as:                                                    ###
### "jobsub_submit -G $group --role=Analysis --memory=$memory --expected-lifetime=$cput   \     ###
###   --disk=$disk --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC,OFFSITE        \     ###
###   --OS=SL6 -N $numjobs --generate-email-summary                                       \     ###
###   file://$ENV{VALOR_EXPT}/prod/jobsub_submit_param_scan.sh                            \     ###
###       --tag $tag  --param_nullhyp1 $param_nullhyp1 --param_nullhyp2 $param_nullhyp2   \     ###
###       --physics-hypothesis $fitparamlist --config-file $config_file  --param1 $param1 \     ###
###       --min-value1 $min_value1 --max-value1 $max_value1 --log-param1 $log_param1      \     ###
###       --npoints1 $npoints1 --param2 $param2 --min-value2 $min_value2                  \     ###
###       --max-value2 $max_value2 --log-param2 $log_param2 --npoints2 $npoints2          \     ###
###       --discovery $discovery --hierarchy $mhchoice --store-spectra $store_spectra     \     ###
###       --fit-option $fitoption --load-prior-from-file $prior_file                      \     ###
###       --pot-override $pot_override --mesg-thresholds $mesg_thresholds                 \     ###
###       --param-values-option $param_values --tarball $tarball --jobdir $job_dir"       \     ###
###       --njobs $numjobs                                                                      ###
###                                                                                             ###
### The first section are jobsub specific parameters:                                           ###
### $group is your experiment e.g. sbnd or dune                                                 ###
### $memory is the memory allocation in MB e.g. 4000                                            ###
### $expected-lifetime e.g. 8h, 45m                                                             ###
### $disk is the disk allocation e.g. 40GB                                                      ###
### $numjobs is the number of jobs you wish to run                                              ###
### $load-prior-from-file default is - but make it none the scipt will change this back         ###
###                                                                                             ###
### The rest of the parameters belong to the script, sent to the node, and  are mainly          ###
### specific to valor. Exceptions, which are not used in run fit are as follows:                ###
###                                                                                             ###
### $tarball: is the valor tarball sent to the node. Please put your $VALOR-DUNE                ###
###           and $VALOR_FRMK in the same directory and tar the file using:                     ###
###           "tar -C /path/to/valor/ -czf valor.tar ."                                         ###
### $job_dir: This is where your output files go and it must be on pnfs.                        ###
### $numjobs: $numjobs should be the same as -N $numjobs in the jobsub section                  ###
###                                                                                             ###
### For Further help ask Dominic Barker                                                         ###
###################################################################################################


#!/bin/sh
#set -x
echo Start `date`
echo Site:${GLIDEIN_ResourceName}
echo "the worker node is " `hostname` "OS: " `uname -a`
echo "You are running as user `whoami`"
echo "The file you are executing is jobsub_submit_param_scan_multifits.sh"
voms-proxy-info --all
export IFDH_CP_MAXRETRIES=2

#Setup the sbnd experiment for ease -- Can change for DUNE GUYS
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
setup fife_utils;

TAG=""
PARAM1=0
NPOINTS1=0
MINVALUE1="none"
MAXVALUE1="none"
LOGPARAM1=0
PARAM2=0
MINVALUE2="none"
NPOINTS2=0
MAXVALUE2="none"
LOGPARAM2=0
DISCOVERY=0
FITPARAMLIST=""
HIERARCHY="" 
FITOPTION=""  
STORESPECTRA=0
SHOULDSTORE="--discard-spectra"
LOADPRORFROMFILE="-"
POTOVERIDE=0
CONFIGFILE=""
MESGTHRESHOLD=""
ONLY1DSCAN=false
SPACE=" "
PARAMVALUESOPTIONS=" "
PARAM1NULLHYP=0
PARAM2NULLHYP=0  
TARBALL=""
INTERACTIVE=0
JOBDIR=""
NJOBS=0

valorfailed=false

while [ $# -gt 0 ]; do
  case "$1" in

    --jobdir )
      if [ $# -gt 1 ]; then
        JOBDIR=$2
        shift
      fi
      ;;

    --tarball )
      if [ $# -gt 1 ]; then
        TARBALL=$2
        shift
      fi
      ;;

    --interactive )
      INTERACTIVE=1
      ;;

    --param_nullhyp1 )
      if [ $# -gt 1 ]; then
        PARAM1NULLHYP=$2
        shift
      fi
      ;;

    --param_nullhyp2 )
      if [ $# -gt 1 ]; then
        PARAM2NULLHYP=$2
        shift
      fi
      ;;


    --tag )
      if [ $# -gt 1 ]; then
        TAG=$2
        shift
      fi
      ;;

    --npoints1 )
      if [ $# -gt 1 ]; then
        NPOINTS1=$2
        shift
      fi
      ;;


    --param1 )
      if [ $# -gt 1 ]; then
        PARAM1=$2
        shift
      fi
      ;;

    --min-value1 )
      if [ $# -gt 1 ]; then
        MINVALUE1=$2
        shift
      fi
      ;;

    --max-value1 )
      if [ $# -gt 1 ]; then
        MAXVALUE1=$2
        shift
      fi
      ;;

    --log-param1 )
      if [ $# -gt 1 ]; then
        LOGPARAM1=$2
        shift
      fi
      ;;

    --npoints2 )
      if [ $# -gt 1 ]; then
        NPOINTS2=$2
        shift
      fi
      ;;


    --param2 )
      if [ $# -gt 1 ]; then
        PARAM2=$2
        shift
      fi
      ;;

    --min-value2 )
      if [ $# -gt 1 ]; then
        MINVALUE2=$2
        shift
      fi
      ;;

    --max-value2 )
      if [ $# -gt 1 ]; then
        MAXVALUE2=$2
        shift
      fi
      ;;

    --log-param2 )
      if [ $# -gt 1 ]; then
        echo "test $2"
        LOGPARAM2=$2
        shift
      fi
      ;;

    --discovery )
      if [ $# -gt 1 ]; then
        DISCOVERY=$2
        shift
      fi
      ;;

    --physics-hypothesis )
      if [ $# -gt 1 ]; then
        FITPARAMLIST=$2
        shift
      fi
      ;;

    --hierarchy )
      if [ $# -gt 1 ]; then
        HIERARCHY=$2
        shift
      fi
      ;;

    --fit-option ) 
      if [ $# -gt 1 ]; then
        FITOPTION=$2
        shift
      fi
      ;;

    --store-spectra )
      if [ $# -gt 1 ]; then
        STORESPECTRA=$2
        shift
      fi
      ;;

    --load-prior-from-file )
      if [ $# -gt 1 ]; then
        LOADPROIRFROMFILE=$2
        shift
      fi
      ;;

    --param-values-option )
      if [ $# -gt 1 ]; then
        PARAMVALUESOPTIONS=$2
        shift
      fi
      ;;


    --pot-override )
      if [ $# -gt 1 ]; then
        POTOVERIDE=$2
        shift
      fi
      ;;

    --config-file )
      if [ $# -gt 1 ]; then
        CONFIGFILE=$2
        shift
      fi
      ;;

    --mesg-thresholds )
      if [ $# -gt 1 ]; then
        MESGTHRESHOLD=$2
        shift
      fi
      ;;

    --njobs )
      if [ $# -gt 1 ]; then
        NJOBS=$2
        shift
      fi
      ;;

      # Other.
      * )
      echo "Invalid Setting $1 with value $2"
      exit 1;
  esac
  shift
done

# Do we want to store the spectra?
if [ $STORESPECTRA == 1 ]; then
  SHOULDSTORE="--store-spectra"
fi

echo "TAG = $TAG"
echo "NPOINTS1 = $NPOINTS1"
echo "PARAM1 = $PARAM1"
echo "MINVALUE1 = $MINVALUE1"
echo "MAXVALUE1 = $MAXVALUE1"
echo "LOGPARAM1 = $LOGPARAM1"
echo "NPOINTS2 = $NPOINTS2"
echo "PARAM2 = $PARAM2"
echo "MINVALUE2 = $MINVALUE2"
echo "MAXVALUE2 = $MAXVALUE2"
echo "LOGPARAM2 = $LOGPARAM2"
echo "DISCOVERY = $DISCOVERY"
echo "PHYSHYPOTHESIS = $FITPARAMLIST"
echo "HIERARCHY = $HIERARCHY" 
echo "FITOPTION = $FITOPTION"
echo "STORESPECTRA = $SHOULDSTORE"
echo "LOADPRORFROMFILE = $LOADPRORFROMFILE"
echo "POTOVERIDE = $POTOVERIDE"
echo "CONFIGFILE = $CONFIGFILE" 
echo "MESGTHRESHOLD $MESGTHRESHOLD" 
echo "TARBALL = $TARBALL"
echo "INTERACTIVE = $INTERACTIVE"
echo "NJOBS = $NJOBS" 

if [ $INTERACTIVE -eq 1 ]; then
  rm -rf $PWD/scratch
  mkdir $PWD/scratch
  _CONDOR_SCRATCH_DIR=$PWD/scratch
  PROCESS=0
  echo "In interactive mode"
fi

#Do Everything in SCRATCH
cd $_CONDOR_SCRATCH_DIR


#First Check that we have a config file 
if [ x$CONFIGFILE = x ]; then
  echo "No configuration was specified."
  exit 1
fi

#Check we have a directory to put things
if [ x$JOBDIR = x ]; then
  echo "No job directory"
  exit 1
fi

#Check we have a tarball
if [ x$TARBALL = x ]; then
  echo "No tarball path"
  exit 1
fi

JOBDIRCUT=${JOBDIR:0:5}

if [ "$JOBDIRCUT" != "/pnfs" ] && ! [ $INTERACTIVE ]; then 
  echo "Job directory is not on dCache. Please fix" 
  exit 1
fi 


TARBALLCUT=${TARBALL:0:5}

if [ "$TARBALLCUT" != "/pnfs" ] && ! [ $INTERACTIVE ]; then 
  echo "tarball  is not on dCache. Please fix" 
  exit 1
fi 

CONFIGCUT=${CONFIGFILE:0:5}

if [ "$CONFIGCUT" != "/pnfs" ] && ! [ $INTERACTIVE ]; then 
  echo "config file  is not on dCache. Please fix" 
  exit 1
fi 

ifdh lss $JOBDIR
if [ $? -ne 0 ]; then
  echo "job dir does not exist bailing."
  exit 73
fi 

#Make a subdirectory
JOBSUBDIR=$JOBDIR/$PROCESS/
echo "JOBSUBDIR = $JOBSUBDIR"

SUBDIRFILES=$(ifdh ls $JOBSUBDIR)
if [ $? -ne 0 ]; then
  echo "making sub dir";
  ifdh mkdir $JOBSUBDIR
else
  for filename in $SUBDIRFILES; do
    ifdh rm $filename
  done
fi

#Fix the proir file work around 
if [ $LOADPRORFROMFILE == "none" ]; then 
  $LOADPRORFROMFILE="-"
fi

## Only look for 2nd parameter if param2 != none
if [ $PARAM2 == "none" ]; then
  $ONLY1DSCAN=true;
fi

#In case someone wants to run the jobsub command not through perl script repeat the catches 
if [ ! $ONLY1DSCAN  ] && [[ $MINVALUE2 == "none" || $MAXVALUE2 == "none" ]]; then 
  echo "If there is a second parameter you must define min_value2 & max_value2, was given $MINVALUE2 and $MAXVALUE2 \n";
  exit 1;
fi

if [ $NPOINTS1 -le 1 ]; then
  echo "npoints1 must be 2 or greater, was given $NPOINTS1";
  exit 1;
fi

if [ "$ONLY1DSCAN" = false  ] && [ $NPOINTS2 -le 0 ]; then
  echo "npoints2 must be 2 or greater, was given $NPOINTS2";
  exit 1;
fi

#Check we are not going to be killed by the log 
if (( $(echo "${MINVALUE1} <= 0" |bc -l) )) && [ $LOGPARAM1 -eq 0 ]; then 
  echo "minval1 cannot be less than or equal to zero for log based calculations"
  exit 1;
fi

if (( $(echo "${MINVALUE2} <= 0" |bc -l) )) && [ $LOGPARAM2 -eq 0 ] && [ "$ONLY1DSCAN" = false  ]; then 
  echo "minval2 cannot be less than or equal to zero for log based calculations"
  exit 1;
fi

#Check we are not going to be killed by the log 
if (( $(echo "${MAXVALUE1} <= 0" |bc -l) )) && [ $LOGPARAM1 -eq 0 ]; then 
  echo "minval1 cannot be less than or equal to zero for log based calculations"
  exit 1;
fi

if (( $(echo "${MAXVALUE2} <= 0" |bc -l) )) && [ $LOGPARAM2 -eq 0 ] && [ "$ONLY1DSCAN" = false  ]; then 
  echo "minval2 cannot be less than or equal to zero for log based calculations"
  exit 1;
fi


if [ $NJOBS -eq 0 ]; then 
  echo "Jobs number not set properly" 
  exit 1 
fi

if [ "$CONFIGCUT" == "/pnfs" ]; then 
  echo "Copying config file"
  timeout 20s ifdh cp -D $CONFIGFILE $_CONDOR_SCRATCH_DIR
  if [ $? -ne 0 ]; then
    echo "Error $? when copying the config file"
    exit 73
  fi

  CONFIGNAME="${CONFIGFILE##*/}"
  CONFIGFILE="config.xml"
  mv $CONFIGNAME config.xml
fi 


#Copy the tarball over 
echo "Copying  tar file"
echo $_CONDOR_SCRATCH_DIR
echo $TARBALL
cd $_CONDOR_SCRATCH_DIR
timeout 200s ifdh cp $TARBALL local.tar
if [ $? -ne 0 ]; then
  echo "Error $? when copying tar file from dCache"
  exit 73
fi

#untar the ball 
tar -xf local.tar

#Setup valor 
cd VALOR-SDK
source setup/setup_at_cvmfs_root6.sh; 
perl configure; 
make; 
cd $_CONDOR_SCRATCH_DIR
cd VALOR-DUNE-LBL
source setup/setup_at_cvmfs_root6.sh;
perl configure;
make;

cd $_CONDOR_SCRATCH_DIR


#The run number to start from for this job 
NMIN=$(perl -E "use POSIX;say floor($PROCESS*$NPOINTS1*$NPOINTS2/$NJOBS)")
echo "NMIN: $NMIN"

#The run number to stop at in this job
NMAX=$(perl -E "use POSIX;say floor(($PROCESS+1)*$NPOINTS1*$NPOINTS2/$NJOBS)")
echo "NMAX: $NMAX" 

STEP1=$(perl -E "say (($MAXVALUE1-$MINVALUE1)/($NPOINTS1-1))")
STEP2=$(perl -E "say (($MAXVALUE2-$MINVALUE2)/($NPOINTS2-1))") 

echo "STEP1 is $STEP1"
echo "STEP2 is $STEP2"

if [ $LOGPARAM1 ];  then
  STEP1=$(perl -E "say ((log($MAXVALUE1/$MINVALUE1)/log(10))/($NPOINTS1-1))")

fi
if [ $LOGPARAM2 ]; then 
  STEP2=$(perl -E "say ((log($MAXVALUE2/$MINVALUE2)/log(10))/($NPOINTS2-1))")
fi



for ((i=$NMIN;i<$NMAX;i++)); 
do 

  #Calculate what the values of the paramers should be for this job. PROCESS is the job number 

  echo "MINVALUE1:$MINVALUE1, i:$i, NPOINTS2:$NPOINTS2, STEP1:$STEP1" 
  VAL1=$(perl -E "use POSIX;say ($MINVALUE1 + floor($i/$NPOINTS2)*$STEP1)")  
  if [ $LOGPARAM1 ];  then
    VAL1=$(perl -E "use POSIX;say (10**($STEP1*(floor($i/$NPOINTS2)) + log($MINVALUE1)/log(10)))")
  fi
  VAL2=0
  if [ $ONLY1DSCAN = false ] ; then 
    VAL2=$(perl -E "say ($MINVALUE2 + ($i%$NPOINTS2)*$STEP2)")
    if [ $LOGPARAM2 ]; then
      VAL2=$(perl -E "say (10**($STEP2*($i%$NPOINTS2) + log($MINVALUE2)/log(10)))")
    fi

    echo "VAL1: $VAL1"
    echo "VAl2: $VAL2"
    echo " " 
  fi

  #
  # Perform all fits
  # Fixing the parameter each value in the scan

  #Make the valor command.
  standard_options=" --run $PROCESS --ntoymc 1  --physics-hypothesis $FITPARAMLIST \
    --param-values-option $PARAMVALUESOPTIONS --fitter $FITOPTION \
    --config $CONFIGFILE --mesg-thresholds $MESGTHRESHOLD --without-stat-fluct \
    $SHOULDSTORE --$HIERARCHY --load-prior-from-file $LOADPRORFROMFILE \
    --pot-override $POTOVERIDE ";

  if [ $DISCOVERY == 0 ]; then 
    arguments="--tag $TAG-$PARAM1-$VAL1-$PARAM2-$VAL2 --first-octant  --fix-at-$PARAM1 $VAL1 --fix-$PARAM1 ";
    arguments=$arguments$standard_options

    if [ "$ONLY1DSCAN" = false ]; then 
      argparm2=" --fix-at-$PARAM2 $VAL2 --fix-$PARAM2 ";
      arguments=$arguments$argparm2
    fi
  else 
    arguments="--tag $TAG-$PARAM1-$VAL1-$PARAM2-$VAL2 --first-octant "
    arguments=$arguments$standard_options
    argparam1=" --input-$PARAM1 $VAL1 \
      --fix-at-$PARAM1 $PARAM1NULLHYP \
      --fix-$PARAM1 ";
    arguments=$arguments$argparam1
    if  [ "$ONLY1DSCAN" = false  ]; then
      arguments=$arguments
      argparam2=" --input-$PARAM2 $VAL2 \
        --fix-at-$PARAM2 $PARAM2NULLHYP \
        --fix-$PARAM2 ";
      arguments=$arguments$argparam2
    fi
  fi

  #Might be able to make this generic? Leave for now
  program="valor_dune-run_fit";
  valorcommand=$program$SPACE$arguments;

  #Run Valor
  echo "Running: $valorcommand";
  eval $valorcommand >> valor.out 
  if [ $? -ne 0 ]; then
    valorfailed=true
    echo "Err $?. Something went wrong with valor please look at the text file"
  fi
done 

#Clean the directory before moving things. I think the scratch will be moved anyways.
timeout 200s ifdh rm  local.tar


#Copy over scratch directory 
for filename in $_CONDOR_SCRATCH_DIR/*; do
  if [ -f "$filename" ]; then
    echo "Copying $filename"
    timeout 200s ifdh cp -D $filename $JOBSUBDIR
  fi
done

if [ $? -ne 0 ]; then
  echo "Error $? when copying output to dCache"
  exit 73
fi

if [ $valorfailed = true ]; then
  exit 1
fi

exit 0





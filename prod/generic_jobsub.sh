#!/bin/sh
#set -x
echo Start `date`
echo Site:${GLIDEIN_ResourceName}
echo "the worker node is " `hostname` "OS: " `uname -a`
echo "You are running as user `whoami`"
voms-proxy-info --all
export IFDH_CP_MAXRETRIES=2

#Setup the sbnd experiment for ease -- Can change for DUNE GUYS
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
setup fife_utils;

SHELLSCRIPT=""
JOBDIR=""
TARBALLSDK=""
TARBALLDUNE=""
INTERACTIVE=0

#Collect command params if set 
while [ $# -gt 0 ]; do
  case "$1" in

      --shellscript )
	  if [ $# -gt 1 ]; then
              SHELLSCRIPT=$2
              shift
	  fi
	  ;;

      --jobdir )
	  if [ $# -gt 1 ]; then
              JOBDIR=$2
              shift
	  fi
	  ;;

      --tarballsdk )
	  if [ $# -gt 1 ]; then
              TARBALLSDK=$2
              shift
	  fi
	  ;;

      --tarballdune )
	  if [ $# -gt 1 ]; then
              TARBALLDUNE=$2
              shift
	  fi
	  ;;

      --interactive )
	  INTERACTIVE=1
	  ;;

      # Other.
      * )
	  echo "Invalid Setting $1 with value $2"
	  exit 1;
  esac
  shift
done

echo "SHELLSCRIPT=$SHELLSCRIPT"
echo "JOBDIR=$JOBDIR"
echo "TARBALLSDK=$TARBALLSDK"
echo "TARBALLDUNE=$TARBALLDUNE"

if [ $INTERACTIVE -eq 1 ]; then
    rm -rf $PWD/scratch
    mkdir $PWD/scratch
    _CONDOR_SCRATCH_DIR=$PWD/scratch
    PROCESS=1
    echo "In interactive mode"
fi

#Make the subdir
SHELLSCRIPTNAME="${SHELLSCRIPT##*/}"
SUBDIR=${SHELLSCRIPTNAME%.*}
JOBSUBDIR=$JOBDIR/$SUBDIR
if [ ! -d $JOBSUBDIR ]; then
    echo "making sub dir";
    ifdh mkdir $JOBSUBDIR
fi

#Check everything is filled
if [ x$TARBALLSDK = x ]; then 
    echo "SDK tarball path empty"
    exit 2
fi

#Check everything is filled
if [ x$TARBALLDUNE = x ]; then 
    echo "DUNE tarball path empty"
    exit 2
fi

#Check the jobdirctory.
if [ x$JOBDIR = x ]; then
    echo "No job directory"
    exit 2
fi

JOBDIRCUT=${JOBDIR:0:5}

if [ "$JOBDIRCUT" != "/pnfs" ] && ! [ $INTERACTIVE ]; then
    echo "Job directory is not on dCache. Please fix"
    exit 2
fi


#Check the shell script location
if [ x$SHELLSCRIPT = x ]; then 
    echo "shell script path empty"
    exit 2
fi 

timeout 200s ifdh cp $SHELLSCRIPT $_CONDOR_SCRATCH_DIR/shellscript.sh

if [ $? -ne 0 ]; then
  echo "Error $? when copying shell script from dCache"
  exit 73
fi

chmod -x shellscript.sh

#Copy the balls over 
echo "Copying  tar files"
echo $_CONDOR_SCRATCH_DIR
echo $TARBALLDUNE
echo $TARBALLSDK
cd $_CONDOR_SCRATCH_DIR
timeout 200s ifdh cp $TARBALLDUNE localDUNE.tar
timeout 200s ifdh cp $TARBALLSDK localSDK.tar

if [ $? -ne 0 ]; then
  echo "Error $? when copying tar files from dCache"
  exit 73
fi

#Untar the balls!
mkdir VALOR-DUNE-LBL
cd ./VALOR-DUNE-LBL
mv ../localDUNE.tar ./
tar -xf localDUNE.tar

cd $_CONDOR_SCRATCH_DIR/
mkdir VALOR-SDK
cd ./VALOR-SDK
mv ../localSDK.tar ./
tar -xf localSDK.tar 

cd $_CONDOR_SCRATCH_DIR/
#Setup valor 
cd VALOR-SDK
source setup/setup_at_cvmfs_root6.sh; 
perl configure; 
make; 
cd $_CONDOR_SCRATCH_DIR/
cd VALOR-DUNE-LBL
source setup/setup_at_cvmfs_root6.sh;
perl configure;
make;

valorfailed=false;
cd $_CONDOR_SCRATCH_DIR
#Run Valor
echo "Running shellscript";
sh shellscript.sh > valor.out 
if [ $? -ne 0 ]; then
    echo "Err $?. Something went wrong with valor please look at the text file"
    valorfailed=true;
fi

#Copy over scratch directory 
for filename in $_CONDOR_SCRATCH_DIR/*; do
    echo "Copying $filename" 
    if [ -f "$filename" ]; then
	timeout 200s ifdh cp -D $filename $JOBSUBDIR
    fi
done 

if [ $? -ne 0 ]; then
  echo "Error $? when copying output to dCache"
  exit 73
fi

if [ $valorfailed = true ]; then
    exit 99
fi

exit 0
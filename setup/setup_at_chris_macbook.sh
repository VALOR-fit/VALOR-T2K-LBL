#!/bin/bash

#
# T2KANA oscillation analysis code setup script for Steve Dennis' laptop
#
# Input argument: oscillation analysis code top-level directory
# If no argument given, just use the directory containing the script
#

if [ -z ${VALOR_FWMK} ] ; then
  if [ $USER == chrisbarry ] ; then 
    #source "/home/steve/analysis/fmwk/setup/setup_at_steve_laptop.sh"
    export VALOR_FMWK="/Users/chrisbarry/software/VaLOR/fmwk/trunk"
    source ${VALOR_FMWK}/setup/setup_at_chris_macbook.sh ${VALOR_FMWK}
  fi
  if [ $USER == phrhbq ] ; then 
    source "/storage/epp2/phrhbq/analysis/fmwk/setup/setup_at_steve_laptop.sh"
  fi
  if [ $USER == jat98257 ] ; then
    source "/opt/ppd/t2k/users/steve/valor/fmwk/trunk//setup/setup_at_cvmfs.sh"
  fi
fi

if [ -z $VALOR_EXPT ]; then
  source="$( greadlink -m "${BASH_SOURCE[0]}" )"
  export VALOR_EXPT="$( dirname "$( dirname "${source}" )" )/" # Directory above the one containing this file
fi

export ROOTSYS=/Users/chrisbarry/software/root_v5.34.34_osx
export PYTHONPATH=${ROOTSYS}/lib:${ROOTSYS}/bin:${ROOTSYS}/bindings/pyroot

echo -e "VALOR_FMWK: \n $VALOR_FMWK"
echo -e "VALOR_EXPT: \n $VALOR_EXPT"

export DYLD_LIBRARY_PATH="${VALOR_EXPT}/lib:${LD_LIBRARY_PATH}"
export PATH="${VALOR_EXPT}/bin:${VALOR_EXPT}/prod:${PATH}"


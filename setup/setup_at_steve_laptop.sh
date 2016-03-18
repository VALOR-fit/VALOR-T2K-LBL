#!/bin/bash

#
# T2KANA oscillation analysis code setup script for Steve Dennis' laptop
#
# Input argument: oscillation analysis code top-level directory
# If no argument given, just use the directory containing the script
#

if [ -z ${VALOR_FWMK} ] ; then
  if [ $USER == steve ] ; then 
    source "/home/steve/analysis/fmwk/setup/setup_at_steve_laptop.sh"
  fi
  if [ $USER == phrhbq ] ; then 
    source "/storage/epp2/phrhbq/analysis/fmwk/setup/setup_at_steve_laptop.sh"
  fi
  if [ $USER == jat98257 ] ; then
    source "/opt/ppd/t2k/users/steve/valor/fmwk/trunk//setup/setup_at_cvmfs.sh"
  fi
fi

export T2KOSCANA=$1

if [ -z $T2KOSCANA ]; then
  source="$( readlink -m "${BASH_SOURCE[0]}" )"
  export T2KOSCANA="$( dirname "$( dirname "${source}" )" )/" # Directory above the one containing this file
fi

export LD_LIBRARY_PATH="${T2KOSCANA}/lib:${LD_LIBRARY_PATH}"
export PATH="${T2KOSCANA}/bin:${T2KOSCANA}/production:${PATH}"


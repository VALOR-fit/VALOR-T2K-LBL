#!/bin/bash

#
# T2KANA oscillation analysis code setup script for Steve Dennis' laptop
#
# Input argument: oscillation analysis code top-level directory
# If no argument given, just use the directory containing the script
#

if [ -z ${VALOR_FWMK} ] ; then
  if [ ! -z "${2}" ] ; then
    source "${2}/setup/setup_at_steve_laptop.sh" "${2}"
  fi
fi

if [ -z ${VALOR_FMWK} ] ; then
  if [ $USER == steve ] ; then 
    source "/home/steve/steve/analysis/github/VALOR-SDK/setup/setup_at_steve_laptop.sh"
  fi
  if [ $USER == steved ] ; then 
    source "/hepstore/steved/analysis/VALOR-SDK/setup/setup_at_cvmfs_root6.sh"
  fi
  if [ $USER == phrhbq ] ; then 
    source "/storage/epp2/phrhbq/analysis/fmwk/setup/setup_at_steve_laptop.sh"
  fi
  if [ $USER == jat98257 ] ; then
    source "/opt/ppd/t2k/users/steve/valor/fmwk/trunk//setup/setup_at_cvmfs.sh"
  fi
fi

if [ ! -z ${1} ] ; then
  export VALOR_EXPT="$(readlink -m ${1})"
fi

if [ -z $VALOR_EXPT ]; then
  source="$( readlink -m "${BASH_SOURCE[0]}" )"
  export VALOR_EXPT="$( dirname "$( dirname "${source}" )" )/" # Directory above the one containing this file
fi

export LD_LIBRARY_PATH="${VALOR_EXPT}/lib:${LD_LIBRARY_PATH}"
export PATH="${VALOR_EXPT}/bin:${VALOR_EXPT}/prod:${PATH}"


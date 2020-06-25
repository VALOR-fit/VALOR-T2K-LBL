#!/bin/sh

#
# VALOR_EXPT oscillation analysis code setup script for Steve Dennis' laptop
#
# Input arguments: experiment code, framework code
#
# Neither is necessary - if the framework has already been setup, it won't use that argument
# If experiment code isn't given, uses the folder containing the folder containing this setup file.
#

if [ -z "${VALOR_FMWK}" ] ; then
  if [ -z "${2}" ] ; then
    echo "VALOR framework must already be set up or specified in argument \$2 ($2)"
    return
  fi
  export VALOR_FMWK="$(readlink -m $2)" # Requires GNU readlink - should work on SL builds (but not Macs - try greadlink)
  
  . "${VALOR_FMWK}/setup/setup_at_cvmfs.sh" "${VALOR_FMWK}"
fi

export VALOR_EXPT="$(readlink -m "$1")"

if [ -z "${VALOR_EXPT}" ]; then
  sdir="$( readlink -m "${BASH_SOURCE[0]}" )"
  export VALOR_EXPT="$( dirname "$( dirname "${sdir}" )" )/" # Directory above the one containing this file
fi

export LD_LIBRARY_PATH="${VALOR_EXPT}/lib:${LD_LIBRARY_PATH}"
export PATH="${VALOR_EXPT}/bin:${VALOR_EXPT}/prod:${PATH}"

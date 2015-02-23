#!/bin/bash

#
# T2KANA oscillation analysis code setup script for Steve Dennis' laptop
#
# Input argument: oscillation analysis code top-level directory
# If no argument given, just use the directory containing the script
#

if [ -z ${VALOR_FWMK} ] ; then
  if [ $USER == steve ] ; then 
    source /home/steve/analysis/fmwk/setup/setup_at_steve_laptop.sh 
  fi
  if [ $USER == phrhbq ] ; then 
    source /home/epp/phrhbq/analysis/fmwk/setup/setup_at_steve_laptop.sh
  fi
fi

export T2KOSCANA=$1

if [ -z $T2KOSCANA ]; then
  source="$( readlink -m "${BASH_SOURCE[0]}" )"
  export T2KOSCANA="$( dirname "$( dirname "${source}" )" )/" # Directory above the one containing this file
fi

export LD_LIBRARY_PATH=${T2KOSCANA}/fitter:${LD_LIBRARY_PATH}
export PATH=${T2KOSCANA}/root:${T2KOSCANA}/production:${PATH}

alias root2k="root2k -l"



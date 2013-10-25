#!/bin/bash

#
# T2KANA oscillation analysis code setup script for Steve Dennis' laptop
#
# Input argument: oscillation analysis code top-level directory
# If no argument given, just use the directory containing the script
#

export T2KOSCANA=$1

if [ -z $T2KOSCANA ]; then
  export T2KOSCANA="$( dirname `cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd `)/" # Directory above the one containing this file
fi

source /home/steve/packages/root/root-v5-34/bin/thisroot.sh

export LD_LIBRARY_PATH=${T2KOSCANA}/fitter:${LD_LIBRARY_PATH}
export PATH=${T2KOSCANA}/root:${T2KOSCANA}/production:${PATH}

alias root2k="root2k -l"


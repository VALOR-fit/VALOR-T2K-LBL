#!/bin/bash

#
# T2KANA oscillation analysis code setup script for Liverpool
#
# Input argument: oscillation analysis code top-level directory
#

export T2KOSCANA=$1
export VALOR_FMWK=$2
export ROOTSYS=/user/cbarry/software/root_v5.34.23
export LOG4CPP=/user/cbarry/software/log4cpp-install
export LD_LIBRARY_PATH=${T2KOSCANA}/lib:${VALOR_FMWK}/lib:${ROOTSYS}/lib:${LOG4CPP}/lib:${LD_LIBRARY_PATH}
export PATH=${T2KOSCANA}/bin:${T2KOSCANA}/root:${T2KOSCANA}/production:${ROOTSYS}/bin:${PATH}
export PYTHONPATH="${ROOTSYS}/lib:${PYTHONPATH}"

echo "T2KOSCANA set as:"
echo $T2KOSCANA
echo "VALOR_FMWK set as:"
echo $VALOR_FMWK


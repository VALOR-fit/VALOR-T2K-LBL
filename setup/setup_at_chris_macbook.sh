#!/bin/bash

#
# T2KANA oscillation analysis code setup script for Chris' MacBookPro
#
# Input arguments: oscillation analysis code top-level directory
#                  framework code top-level directory
#

export T2KOSCANA=$1
export VALOR_FMWK=$2
export ROOTSYS=/Users/chrisbarry/software/root_v5.34.34_osx
export PYTHONPATH=${ROOTSYS}/lib:${ROOTSYS}/bin:${ROOTSYS}/bindings/pyroot

echo -e "VALOR_FMWK: \n $VALOR_FMWK"
echo -e "T2KOSCANA: \n $T2KOSCANA"

export DYLD_LIBRARY_PATH=${VALOR_FMWK}/lib:${DYLD_LIBRARY_PATH}
export PATH=${VALOR_FMWK}/root:${PATH}

export DYLD_LIBRARY_PATH=${T2KOSCANA}/lib:${ROOTSYS}/lib:${DYLD_LIBRARY_PATH}
#export DYLD_LIBRARY_PATH=${T2KOSCANA}/fitter:${ROOTSYS}/lib:${DYLD_LIBRARY_PATH}
export PATH=${T2KOSCANA}/bin:${T2KOSCANA}/root:${T2KOSCANA}/production:${ROOTSYS}/bin:${PATH}
#export PATH=${T2KOSCANA}/root:${T2KOSCANA}/production:${ROOTSYS}/bin:${PATH}


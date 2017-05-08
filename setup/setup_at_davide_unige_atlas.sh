#!/bin/bash

#
# T2KANA oscillation analysis code setup script for Costas' MacBookPro
#
# Input argument: oscillation analysis code top-level directory
#

export T2KOSCANA=$1

export LD_LIBRARY_PATH=${T2KOSCANA}/fitter:${ROOTSYS}/fitter:${ROOTSYS}/lib:$T2KOSCANA/lib:${LD_LIBRARY_PATH}:/usr/lib64/
export DYLD_LIBRARY_PATH=${T2KOSCANA}/fitter:${ROOTSYS}/lib:$T2KOSCANA/lib:${DYLD_LIBRARY_PATH}
export PATH=${T2KOSCANA}/root:${T2KOSCANA}/production:${ROOTSYS}/bin:$T2KOSCANA/bin:${PATH}

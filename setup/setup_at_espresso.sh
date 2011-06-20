#!/bin/bash

#
# T2KANA oscillation analysis code setup script for Costas' MacBookPro
#
# Input argument: oscillation analysis code top-level directory
#

export T2KOSCANA=$1
#export ROOTSYS=/usr/local/hepsoft/root/v5.26.00b/
export ROOTSYS=/usr/local/hepsoft/root/v5.28.00a/

export DYLD_LIBRARY_PATH=${T2KOSCANA}/fitter:${ROOTSYS}/lib:${DYLD_LIBRARY_PATH}
export PATH=${T2KOSCANA}/root:${T2KOSCANA}/production:${ROOTSYS}/bin:${PATH}

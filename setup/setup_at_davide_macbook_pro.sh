#!/bin/bash

#
# T2KANA oscillation analysis code setup script for Costas' MacBookPro
#
# Input argument: oscillation analysis code top-level directory
#

#ROOTwork_nd280

export T2KOSCANA=$1
#export ROOTSYS=/Users/davidesgalaberna/root/

export DYLD_LIBRARY_PATH=${T2KOSCANA}/fitter:${ROOTSYS}/lib:${DYLD_LIBRARY_PATH}
export PATH=${T2KOSCANA}/root:${T2KOSCANA}/production:${ROOTSYS}/bin:${PATH}:

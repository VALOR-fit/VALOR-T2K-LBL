#!/bin/bash

#
# T2KANA oscillation analysis code setup script for Toms' MacBookPro
#
# Input argument: top-level directory
#

export T2KOSCANA=$1
export ROOTSYS=/Users/thomasdealtry/T2K/software/ROOT/root-v5-34-00/

export DYLD_LIBRARY_PATH=${T2KOSCANA}/fitter:${ROOTSYS}/lib:${DYLD_LIBRARY_PATH}
export PATH=${T2KOSCANA}/root:${T2KOSCANA}/production:${ROOTSYS}/bin:${PATH}

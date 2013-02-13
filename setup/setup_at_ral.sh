#!/bin/bash

#
# T2KANA oscillation analysis code setup script for RAL/PPD Tier2
#
# Input argument: oscillation analysis code top-level directory
#

export T2KOSCANA=$1
export ROOTSYS=/opt/ppd/t2k/softw/ROOT/builds/SL5_64bit/v5.34.03/
export LD_LIBRARY_PATH=${T2KOSCANA}/fitter:${ROOTSYS}/lib:${LD_LIBRARY_PATH}
export PATH=${T2KOSCANA}/root:${T2KOSCANA}/production:${ROOTSYS}/bin:${PATH}


#!/bin/bash

#
# T2KANA oscillation analysis code setup script for RAL/PPD Tier2
#
# Input argument: oscillation analysis code top-level directory
#

export T2KOSCANA=$1
export VALOR_FMWK=$2
export LOG4CPP=/opt/ppd/t2k/softw/LOG4CPP/builds/SL6.x86_64/v1.1/src/.libs/
export ROOTSYS=/opt/ppd/t2k/softw/ROOT/builds/SL6.x86_64/v5.34.30/
export LD_LIBRARY_PATH=${T2KOSCANA}/fitter:${VALOR_FMWK}/lib:${ROOTSYS}/lib:${LOG4CPP}/src/.libs/:${LD_LIBRARY_PATH}
export PATH=${T2KOSCANA}/root:${T2KOSCANA}/production:${ROOTSYS}/bin:${PATH}
export PYTHONPATH="${ROOTSYS}/lib:${PYTHONPATH}"


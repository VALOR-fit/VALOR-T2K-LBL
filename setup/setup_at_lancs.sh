#!/bin/bash

#
# T2KANA oscillation analysis code setup script for RAL/PPD Tier2
#
# Input argument: oscillation analysis code top-level directory
#

export T2KOSCANA=$1
export VALOR_FMWK=$2
export LOG4CPP=/t2k/Software/GENIE/support64/log4cpp
export ROOTSYS=/home/t2k/ng/root_5_34_08_minuit2/
export LD_LIBRARY_PATH=${T2KOSCANA}/fitter:${VALOR_FMWK}/lib:${ROOTSYS}/lib:${LOG4CPP}/lib/:${LD_LIBRARY_PATH}
export PATH=${T2KOSCANA}/root:${T2KOSCANA}/production:${ROOTSYS}/bin:${PATH}
export PYTHONPATH="${ROOTSYS}/lib:${PYTHONPATH}"


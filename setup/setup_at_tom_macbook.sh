#!/bin/bash

#
# VALOR_EXPT oscillation analysis code setup script for Costas' MacBookPro
#
# Input arguments: 
#   1) oscillation analysis code top-level directory
#   2) valor framework code top-level directory
#

export VALOR_EXPT=$1
export VALOR_FMWK=${VALOR_FMWK}
export ROOTSYS=/Users/thomasdealtry/T2K/software/ROOT/root-v5-34-00/

export DYLD_LIBRARY_PATH=${VALOR_EXPT}/lib:${VALOR_FMWK}/lib:${ROOTSYS}/lib:${DYLD_LIBRARY_PATH}
export PATH=${VALOR_EXPT}/root:${VALOR_EXPT}/production:${ROOTSYS}/bin:${PATH}

#!/bin/bash

#
# VALOR_EXPT oscillation analysis code setup script for for Costas' MacBookPro2
#
# Make sure VALOR_FMWK is setup first. This also determines the version of ROOT used in the build.
# If GENIE is needed, make sure it setup before building the VALOR/T2K code.
# Make sure that GENIE and VALOR and built against the same version of ROOT
#   1) oscillation analysis code top-level directory
#   2) valor framework code top-level directory
#


export VALOR_EXPT="${1}"
if [ -z "$VALOR_EXPT" ]; then
    source="$( greadlink -m "${BASH_SOURCE[0]}" )"
    export VALOR_EXPT="$( dirname "$( dirname "${source}" )" )/" # Directory above the one containing this file
fi

#export VALOR_EXPT=$1
#export VALOR_FMWK=$2
#export ROOTSYS=/opt/ppd/t2k/softw/ROOT/builds/SL6.x86_64/v5.34.30/
#export LD_LIBRARY_PATH=${VALOR_EXPT}/lib:${VALOR_FMWK}/lib:${ROOTSYS}/lib:${LD_LIBRARY_PATH}:/opt/ppd/t2k/softw/LOG4CPP/builds/SL5_64bit/v1.0/src/.libs/:/usr/lib64/
#export PATH=${VALOR_EXPT}/bin:${VALOR_EXPT}/root:${VALOR_EXPT}/production:${ROOTSYS}/bin:${PATH}

export LD_LIBRARY_PATH=${VALOR_EXPT}/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=/opt/ppd/t2k/softw/LOG4CPP/builds/SL6.x86_64/v1.1/src/.libs/:/usr/lib64/:${LD_LIBRARY_PATH}
export PATH=${VALOR_EXPT}/bin:${VALOR_EXPT}/root:${VALOR_EXPT}/production:${PATH}


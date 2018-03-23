#!/bin/bash

#
# T2KANA oscillation analysis code setup script for Compute Canada's cedar cluster (uses setup_at_cvmfs.sh)
#
# Input arguments: analysis code directory, framework directory
#

export T2KOSCANA="$(readlink -m $1)"
export VALOR_FMWK="$(readlink -m $2)"

source "${VALOR_FMWK}/setup/setup_at_cvmfs.sh" "${VALOR_FMWK}"
export LD_LIBRARY_PATH="${T2KOSCANA}/lib:${LD_LIBRARY_PATH}:${T2KOSCANA}/setup/lib_cedar"
export PATH="${T2KOSCANA}/bin:${T2KOSCANA}/production:${PATH}"

#export PYTHONHOME="/cvmfs/fermilab.opensciencegrid.org/products/larsoft/python/v2_7_10/Linux64bit+3.10-2.17/:/cvmfs/fermilab.opensciencegrid.org/products/larsoft/python/v2_7_10/Linux64bit+3.10-2.17/"
#module load gsl

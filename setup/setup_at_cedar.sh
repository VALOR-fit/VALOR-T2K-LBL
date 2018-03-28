#!/bin/bash

#
# T2KANA oscillation analysis code setup script for Compute Canada's cedar cluster
#
# Input arguments: analysis code directory, framework directory
#

export T2KOSCANA="$(readlink -m $1)"
export VALOR_FMWK="$(readlink -m $2)"

source "${VALOR_FMWK}/setup/setup_at_cedar_no_CVMFS.sh" "${VALOR_FMWK}"
export LD_LIBRARY_PATH="/home/fbench/project/fbench/log4cpp/build/lib:${T2KOSCANA}/lib:${LD_LIBRARY_PATH}"
export PATH="${T2KOSCANA}/bin:${T2KOSCANA}/production:${PATH}"
export OPT_WITH_CXX_FLAGS=-fWno-deprecated

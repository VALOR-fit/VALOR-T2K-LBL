#!/bin/bash

#
# DUNE oscillation analysis code setup script for Compute Canada's cedar cluster
#
# Input arguments: analysis code directory, framework directory
#

export VALOR_EXPT="$(readlink -m $1)"
export VALOR_FMWK="$(readlink -m $2)"

log4cpp_lib_path=$(find $HOME -name "liblog4cpp.so" -print -quit)
log4cpp_lib_path=$(dirname "$log4cpp_lib_path")

source "${VALOR_FMWK}/setup/setup_at_cedar_root6.sh" "${VALOR_FMWK}"
export LD_LIBRARY_PATH="$log4cpp_lib_path:${VALOR_EXPT}/lib:${LD_LIBRARY_PATH}"
export PATH="${VALOR_EXPT}/bin:${VALOR_EXPT}/prod:${PATH}"
export OPT_WITH_CXX_FLAGS=-fWno-deprecated

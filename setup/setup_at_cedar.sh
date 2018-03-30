#!/bin/bash

#
# T2KANA oscillation analysis code setup script for Compute Canada's cedar cluster
#
# Input arguments: analysis code directory, framework directory
#

export T2KOSCANA="$(readlink -m $1)"
export VALOR_FMWK="$(readlink -m $2)"

log4cpp_lib_path=$(find $HOME -name "liblog4cpp.so" -print -quit)
log4cpp_lib_path=$(dirname "$log4cpp_lib_path")

source "${VALOR_FMWK}/setup/setup_at_cedar.sh" "${VALOR_FMWK}"
export LD_LIBRARY_PATH="$log4cpp_lib_path:${T2KOSCANA}/lib:${LD_LIBRARY_PATH}"
export PATH="${T2KOSCANA}/bin:${T2KOSCANA}/production:${PATH}"
export OPT_WITH_CXX_FLAGS=-fWno-deprecated

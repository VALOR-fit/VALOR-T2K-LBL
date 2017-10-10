#!/bin/bash

#
# T2KANA oscillation analysis code setup script for any grid site
#
# Input arguments: analysis code directory, framework directory
#

echo "Setting up $T2KOSCANA at IN2P3"

export T2KOSCANA="$(readlink -m $1)"
export VALOR_FMWK="$(readlink -m $2)"

source "${VALOR_FMWK}/setup/setup_at_cvmfs_in2p3.sh" "${VALOR_FMWK}"
export LD_LIBRARY_PATH="${T2KOSCANA}/lib:${LD_LIBRARY_PATH}"
export PATH="${T2KOSCANA}/bin:${T2KOSCANA}/production:${PATH}"


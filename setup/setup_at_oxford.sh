#!/bin/bash

#
# T2KANA oscillation analysis code setup script for Oxford batch farm
#
# Input argument: oscillation analysis code top-level directory
#



#source /data/t2k/raj/software/root/bin/thisroot.sh
#source /data/t2k/dealtry/software/v10r11p9/ROOT/v5r30p02n01/amd64_linux26/bin/thisroot.sh
#source /system/SL5/cern/root/x86_64/OxfordSetup-5.32.01.sh


export T2KOSCANA="$(readlink -m $1)"
export VALOR_FMWK="$(readlink -m $2)"

source "${VALOR_FMWK}/setup/setup_at_oxford.sh" "${VALOR_FMWK}"
export LD_LIBRARY_PATH="${T2KOSCANA}/fitter:${T2KOSCANA}/lib:${LD_LIBRARY_PATH}"
export PATH="${T2KOSCANA}/root:${T2KOSCANA}/production:${PATH}"



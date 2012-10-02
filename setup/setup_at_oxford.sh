#!/bin/bash

#
# T2KANA oscillation analysis code setup script for Oxford batch farm
#
# Input argument: oscillation analysis code top-level directory
#

export T2KOSCANA=$1

source /data/t2k/dealtry/software/v10r11p9/ROOT/v5r30p02n01/amd64_linux26/bin/thisroot.sh
#source /system/SL5/cern/root/x86_64/OxfordSetup-5.32.01.sh

export LD_LIBRARY_PATH=${T2KOSCANA}/fitter:${LD_LIBRARY_PATH}
export PATH=${T2KOSCANA}/root:${T2KOSCANA}/production:${PATH}



#!/bin/bash

#
# T2KANA oscillation analysis code setup script for Oxford batch farm
#
# Input argument: oscillation analysis code top-level directory
#

export T2KOSCANA=$1

source /system/SL5/cern/root/x86_64/OxfordSetup-current-pro.sh

export LD_LIBRARY_PATH=${T2KOSCANA}/fitter:${LD_LIBRARY_PATH}
export PATH=${T2KOSCANA}/root:${T2KOSCANA}/production:${PATH}



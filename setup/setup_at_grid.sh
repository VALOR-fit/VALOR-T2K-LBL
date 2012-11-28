#!/bin/bash 

#
# T2KANA oscillation analysis code setup script for RAL/PPD Tier2
#
# Input argument: oscillation analysis code top-level directory
#

export T2KOSCANA=$PWD/trunk_GRID/
echo ${T2KOSCANA}

#export ROOTSYS=$VO_T2K_ORG_SW_DIR/nd280v10r11p17/ROOT/v5r30p02n01/
export ROOTSYS=$VO_T2K_ORG_SW_DIR/nd280v9r9p1/ROOT/v5r28p00n01/
echo ${ROOTSYS}

export LD_LIBRARY_PATH=${T2KOSCANA}/fitter:${ROOTSYS}/lib:${LD_LIBRARY_PATH}

export PATH=${T2KOSCANA}/root:${T2KOSCANA}/production:${ROOTSYS}/bin:${PATH}

#!/bin/bash 

#
# VALOR_EXPT oscillation analysis code setup script for RAL/PPD Tier2
#
# Input argument: oscillation analysis code top-level directory
#

export VALOR_EXPT=$PWD/trunk_GRID/
echo ${VALOR_EXPT}

#export ROOTSYS=$VO_T2K_ORG_SW_DIR/nd280v10r11p17/ROOT/v5r30p02n01/amd64_linux26/
#export ROOTSYS=$VO_T2K_ORG_SW_DIR/nd280v9r9p1/ROOT/v5r28p00n01/
unset ROOTSYS
source ${VO_T2K_ORG_SW_DIR}/nd280v10r11p17/ROOT/v5r30p02n01/cmt/setup.sh
if [ ! -n "$ROOTSYS" ]
    then
      #queen mary hack
      source ${VO_T2K_ORG_SW_DIR}/../../sl5/t2k.org/nd280v10r11p17/ROOT/v5r30p02n01/cmt/setup.sh
fi
if [ ! -n "$ROOTSYS" ]
    then
      source ${VO_T2K_ORG_SW_DIR}/nd280v10r11p21/ROOT/v5r30p02n01/cmt/setup.sh;
fi
if [ ! -n "$ROOTSYS" ]
    then
      source ${VO_T2K_ORG_SW_DIR}/nd280v10r11p25/ROOT/v5r30p02n01/cmt/setup.sh;
fi
if [ ! -n "$ROOTSYS" ]
    then
      source ${VO_T2K_ORG_SW_DIR}/nd280v10r11p31/ROOT/v5r30p02n01/cmt/setup.sh;
fi
if [ ! -n "$ROOTSYS" ]
    then
      source ${VO_T2K_ORG_SW_DIR}/nd280v11r11/ROOT/v5r30p02n01/cmt/setup.sh;
fi
if [ ! -n "$ROOTSYS" ]
    then
      echo "ROOT not found";
fi


echo ${ROOTSYS}

export LD_LIBRARY_PATH=${VALOR_EXPT}/lib:${ROOTSYS}/lib:${LD_LIBRARY_PATH}

export PATH=${VALOR_EXPT}/root:${VALOR_EXPT}/production:${ROOTSYS}/bin:${PATH}

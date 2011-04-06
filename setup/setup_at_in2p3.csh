#!/bin/tcsh

#
# T2KANA oscillation analysis code setup script for IN2P3 farm
#
# Input argument: oscillation analysis code top-level directory
#

setenv T2KOSCANA $1
setenv ROOTSYS /usr/local/root/root_v5.26.00b/root
setenv LD_LIBRARY_PATH ${T2KOSCANA}/fitter:${ROOTSYS}/lib:${LD_LIBRARY_PATH}
setenv PATH ${T2KOSCANA}/root:${T2KOSCANA}/production:${ROOTSYS}/bin:${PATH}

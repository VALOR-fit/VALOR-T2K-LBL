#!/bin/bash

export T2KOSCANA=${PWD}
export ROOTSYS=/opt/ppd/t2k/GENIE/builds/SL5_64bit/ext/root-v5.28.00a/

export LD_LIBRARY_PATH=${T2KOSCANA}/fitter:${ROOTSYS}/lib:${LD_LIBRARY_PATH}
export PATH=${T2KOSCANA}/root:${T2KOSCANA}/production:${ROOTSYS}/bin:${PATH}


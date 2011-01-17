#!/bin/tcsh

setenv ROOTSYS /usr/local/root/root_v5.26.00b/root
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${ROOTSYS}/lib:.
setenv PATH ${PATH}:${ROOTSYS}/bin

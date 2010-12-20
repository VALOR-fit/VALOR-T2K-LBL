#!/bin/bash

source /opt/ppd/t2k/GENIE/builds/SL5_64bit/v2.7.1-setup

export T2KANA=${PWD}/../../

export PATH=${T2KANA}/bin:$PATH
export LD_LIBRARY_PATH=${T2KANA}/lib:${LD_LIBRARY_PATH}

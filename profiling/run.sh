#!/bin/bash

export GPERF=/Users/thomasdealtry/T2K/ralrepo/t2k/t2kana/nuosc/trunk/gperftools-2.1

g++ -m64 -pipe -W -Wall -fsigned-char -fno-common -flat_namespace -Wall -fPIC -O2   -Wno-strict-aliasing -ffriend-injection -c test_sk_ereco_spectrum_twk_sys_throw.C -I$T2KOSCANA/fitter/ -I$ROOTSYS/include/ -I$GPERF/include/install_dir/include/gperftools/

g++ -g test_sk_ereco_spectrum_twk_sys_throw.o -m64 -pipe -W -Wall -fsigned-char -fno-common -flat_namespace -L$T2KOSCANA/fitter/ -lT2KNuOsc -L$ROOTSYS/lib/ -lGui -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMinuit -lMinuit2 -lXMLIO -lXMLParser -L$GPERF/install_dir/lib/ -lprofiler -o throw.exe 

CPUPROFILE=exe.prof ./throw.exe

$GPERF/install_dir/bin/pprof --pdf throw.exe exe.prof > test.pdf


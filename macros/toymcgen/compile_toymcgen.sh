# Compiles generate_toymc_ensemble.C macro
# Works on Linux/g++ with log4cpp in library path.
# Other configurations may need some changes, but hopefully OK

command="g++ -pipe -W -Wall -fsigned-char -fno-common -Wall -fPIC -O2   -Wno-strict-aliasing -ffriend-injection -c $T2KOSCANA/macros/toymcgen/generate_toymc_ensemble_exe.C -I$T2KOSCANA/fitter/ -I${VALOR_FMWK}/src/ `root-config --cflags`"

echo $command
$command

# Order of lirbaries is important here. VALOR libs must come before ROOT/logcpp
command="g++ -W -Wall -fsigned-char  -O2   -Wno-strict-aliasing  generate_toymc_ensemble_exe.o -L$T2KOSCANA/fitter/ -L${VALOR_FMWK}/lib/  -lValorConv -lValorFitCore -lValorMesg -lValorParse -lValorUtil -lT2KNuOsc -o generate_toymc_ensemble.exe  -llog4cpp `root-config --evelibs --ldflags` -lMinuit -lMinuit2 -lHistPainter -lXMLParser -lXMLIO "
echo $command
$command

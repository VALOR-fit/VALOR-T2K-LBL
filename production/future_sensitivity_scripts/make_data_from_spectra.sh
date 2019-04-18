script="$T2KOSCANA/macros/spectra/make_data_from_spectra.py"
fake="$1"
outname="$2"

if [[ "$1" != "" && "$2" != "" ]]; then
   python $script -f $fake -o $outname
else
   echo -e "Usage: make_data_from_spectra.py <spectra_prefix> <output_filename>"
   echo -e "   spectra_prefix: The spectra filename prefix excluding the sample"
   echo -e "   output_filename: The output data filename"
fi

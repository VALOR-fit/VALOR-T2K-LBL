dir=$1

deltacpval="-1.601"
mixing13val="0.021"
mixing23val_nh="0.528"
mixing23val_ih=$mixing23val_nh
dmsq23val="0.003"
if [ "$2" == "asimov_fhc" ]; then
   asimov="asimov_fhc"
elif [ "$2" == "asimov_rhc" ]; then
   asimov="asimov_rhc"
elif [ "$2" == "asimov_standard" ]; then
   asimov="asimov_standard"
elif [ "$2" == "hybrid_fhc" ]; then
   asimov="hybrid_fhc"
elif [ "$2" == "hybrid_rhc" ]; then
   asimov="hybrid_rhc"
fi

if [ "$3" == "t2k" ]; then
   constraint="t2k"
elif [ "$3" == "reactor" ]; then
   constraint="reactor"
fi

par0="$4"
par1="$5"

if [ "$par0" == "deltacp" ]; then
   par0val_nh=$deltacpval
   par0val_ih=$deltacpval
fi

if [ "$par1" == "deltacp" ]; then
   par1val=$deltacpval
fi

dir="~/analysis/pac/results/fits/${dir}/${asimov}"
path="$T2KOSCANA/macros/misc/draw_contour"

if [ "$par1" == "" ]; then
   nhfilepre="${dir}/fit_oa2018_normal_${constraint}_${par0}_nexp_1_${par0val_nh}_-999.000-job0_-_0_0"
   ihfilepre="${dir}/fit_oa2018_inverted_${constraint}_${par0}_nexp_1_${par0val_ih}_-999.000-job0_-_0_0"
   root2k -l -q -b "${path}/get_contour.C(\"${nhfilepre}.root\")"
   root2k -l -q -b "${path}/get_contour.C(\"${ihfilepre}.root\")"
   draw_reactor=0
   global=1
   outfile="${par0}_${asimov}_${constraint}_global"
   script="CompareContour1D.C"
   root2k -l -q -b "${path}/${script}(\"${nhfilepre}_contour.root\",\"${ihfilepre}_contour.root\",\"${outfile}\",\"Normal\",\"Inverted\",${global},1,1,0,1,0,0,0,${draw_reactor},1)"
fi


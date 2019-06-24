#!/bin/bash

[ "$#" -eq 2 ] || (echo "2 arguments required, $# provided. e.g. produce_osc3inputs.sh 190606 <process|finish|checktex|checkroot|checkrootcopy>" && exit)

today=$1


if [ "$2" = "checkrootcopy" ]; then
    rm -f $T2KOSCANA/runs-lbl/inputsforosc3pp/$today/nodiff.txt
fi
for type in nominal 8 11a 11b flux_lukas; do
    if [ "$type" = "nominal" ]; then
	config=analysis_configuration_hkfluxsophie_10year_postfit.xml
	export systtag=0
    elif [ "$type" = "8" ]; then
	config=analysis_configuration_hkfluxsophie_10year_postfit_morenueflux.xml
	export systtag=8
    elif [ "$type" = "11a" ]; then
	config=analysis_configuration_hkfluxsophie_10year_postfit_escale29.xml
	export systtag=11a
    elif [ "$type" = "11b" ]; then
	config=analysis_configuration_hkfluxsophie_10year_postfit_escale19.xml
	export systtag=11b
    elif [ "$type" = "flux_lukas" ]; then
	config=analysis_configuration_hkfluxsophie_10year_postfit_lukasflux.xml
	export systtag=flux_lukas
    else
	echo Unknown option $type
	continue
    fi
    dir=$T2KOSCANA/runs-lbl/inputsforosc3pp/$today/$systtag
    echo $dir
    if [ ! -d "$dir" ]; then
	mkdir -p $dir
    fi
    cd $dir
    if [ "$2" = "process" ]; then
	for sample in 1Rmu 1Re RHC1Re RHC1Rmu; do
	    for x in true false; do
		if [ "$x" = "true"]; then
		    fijtag = "origsyst"
		else
		    fijtag = ""
		fi
		run_on_condor.py -e $T2KOSCANA/root/root2k -a '-b -q $T2KOSCANA/macros/spectra/test_systematic_error_propagation_all_2015_joint_2dtempl.C+(0,"'$sample'","3f_osc",0.32,0.0256584,0.5,0,7.6E-5,2.4E-3,false,"'$config'","nuosc3f+nuisance","sinsqtheta","normal",-999,"f_tamesu",false,true,"postbanff_osc_",true,'$x','$fijtag')' -l valorsyst_${systtag}_${sample}_${fijtag} --escape-quotes;
	done
    elif [ "$2" = "finish" ]; then
	mkdir -p ~/public_html/oa/$today/$systtag
	python $T2KOSCANA/macros/inputs/find_cov_matrices_from_configfile.py -c $config -d $today
	cp *pdf *.fij.t2k.root ~/public_html/oa/$today/$systtag
	for f in *tex; do
	    newfile=~/public_html/oa/$today/$systtag/$f
	    echo \\documentclass[10pt]{article} >> $newfile
	    echo \\usepackage[table,x11names,dvipsnames,table]{xcolor} >> $newfile
	    echo \\usepackage{siunitx} >> $newfile
	    echo \\usepackage{booktabs} >> $newfile
	    echo \\begin{document} >> $newfile
	    cat $f >> ~/public_html/oa/$today/$systtag/$f
	    echo \\end{document} >> ~/public_html/oa/$today/$systtag/$f
	done
    elif [ "$2" = "checkroot" ]; then
	if [ "$type" = "nominal" ]; then
	    continue;
	fi
	for f in ../0/*fij.t2k{,.true,.mode}.root; do
	    echo $f
	    if [ -f "$(basename $f)" ]; then
		for checker in syserre_ syserretrue_ truespectrum_ recospectrum_; do
		    if [[ $f == *"fijorigsyst.t2k.root"* ]]; then
			if [ "$checker" = "syserretrue_" ]; then
			    continue
			fi
		    fi
		    if [[ $checker = "syserre"* ]]; then
			subtract=true
		    else
			subtract=false
		    fi
		    if [[ "$(basename $f)" = *"1Re"* ]]; then
			maxx=1.25
		    else
			maxx=29.99
		    fi
		    run_on_condor.py -e $T2KOSCANA/root/root2k -a '-b -q $T2KOSCANA/macros/misc/compare/compare_1dhist_all.C("'$f'","'$(basename $f)'","'$(basename $f)_$checker'","'$checker'","XXXXXXXX",1.0,1E-6,false,false,-99999999,'$maxx',true,false,'$subtract')' -l systdiff_$(basename $f)_$checker --escape-quotes;
		done
	    fi
	done
    elif [ "$2" = "checkrootcopy" ]; then
	if [ "$type" = "nominal" ]; then
	    continue;
	fi
	for f in comparisons*.pdf; do
	    filesize=`stat -c%s $f`
	    if [ "$filesize" -le "10000" ]; then
		echo $type/$(basename $f)/$checker $filesize >> ../nodiff.txt
	    else
		\cp $f ~/public_html/oa/$today/$type
	    fi
	done
    elif [ "$2" = "checktex" ]; then
	if [ "$type" = "nominal" ]; then
	    continue;
	fi
	cd $dir/..
	for f in 0/*.tex; do
	    echo $(basename $f);
	    diff 0/$(basename $f) $type/$(basename $f) | colordiff.pl;
	done
    fi
done

if [ "$2" = "checkrootcopy" ]; then
    if [ -f ../nodiff.txt ]; then
	echo "No differences in these files"
	cat $T2KOSCANA/runs-lbl/inputsforosc3pp/$today/nodiff.txt
    fi
fi
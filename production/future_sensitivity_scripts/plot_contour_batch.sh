file_pattern=$1
tag=$2

shares=(FHC FHC_RHC)
POTs=(7.8E21 10E21 15E21)
types=(hybrid pure_AsimovA)
hierarchy=(normal inverted)

for t in ${types[@]}; do
	for s in ${shares[@]}; do
		for p in ${POTs[@]}; do
			for h in ${hierarchy[@]}; do
				filenames=""
				for b in `seq 0 3`; do
					f=${file_pattern/\[BATCH\]/$b}
					f=${f/\[SHARE\]/$s}
					f=${f/\[POT\]/$p}
					f=${f/\[TYPE\]/$t}
					f=${f/\[HIER\]/$h}
					filenames="$filenames $f"
				done
				labels="batch0 batch1 batch2 batch3"
				python $T2KOSCANA/macros/contour/plot_contours_TN.py --filename ${filenames} --label ${labels} --t2k --output contour_${2}_${h}_${s}_${p}_${t}.pdf	
			done
		done
	done
done






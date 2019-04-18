script=$T2KOSCANA/macros/contour/plot_contours.py

constraint=$1
if [ "$constraint" == "t2k" ]; then
   t2k_flag="--t2k"
else
   t2k_flag=""
fi
std_file=fit_oa2018_normal_${constraint}_deltacp_nexp_1_-1.601_-999.000-job0_-_0_0_contour.root

# Asimovs

file1="~/analysis/pac/results/fits/deltacp/asimov_standard/$std_file"
file2="~/analysis/pac/results/fits/deltacp/asimov_fhc/$std_file"
file3="~/analysis/pac/results/fits/deltacp/asimov_rhc/$std_file"

input_files="$file1 $file2 $file3 $file4"
labels="Run1-9 Run1-10FHC Run1-10RHC"

# Specify normal hierarchy by default, adding the global tag will produce inverted as wel
python $script -f $input_files -l $labels -o pac_asimov_${constraint}_normal.pdf --global $t2k_flag

# Hybrids

file1="/data/neutrino/t2k/phrdqd/analysis/run1_9/fits/deltacp/data/$std_file"
file2="~/analysis/pac/results/fits/deltacp/hybrid_fhc/$std_file"
file3="~/analysis/pac/results/fits/deltacp/hybrid_rhc/$std_file"

input_files="$file1 $file2 $file3"
labels="Run1-9 Run1-10FHC Run1-10RHC"

# Specify normal hierarchy by default, adding the global tag will produce inverted as wel
python $script -f $input_files -l $labels -o pac_hybrid_${constraint}_normal.pdf --global $t2k_flag


all_study_type=(pure_AsimovA hybrid)
all_allocation=( FHC FHC_RHC )
all_POT=( 7.8E21 10E21 15E21)

all_marg_dir=$1
all_hier=$2
all_reac=$3
all_job_dir=$4
all_par0=$5
all_par1=$6

for all_s in "${all_study_type[@]}"; do
    for all_a in "${all_allocation[@]}"; do
        for all_p in "${all_POT[@]}"; do
            . ${T2KOSCANA}/production/future_sensitivity_scripts/fit.sh ${all_s} ${all_a} ${all_p} ${all_marg_dir} ${all_hier} ${all_reac} ${all_job_dir} ${all_par0} ${all_par1}
        done
    done
done

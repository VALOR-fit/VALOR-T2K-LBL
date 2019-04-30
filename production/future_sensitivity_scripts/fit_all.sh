all_study_type=(pure_AsimovA hybrid)
all_allocation=( FHC FHC_RHC )
all_POT=( 7.8E21 10E21 15E21)
all_hier=(normal inverted)
#all_study_type=(hybrid)
#all_allocation=(FHC_RHC)
#all_POT=(7.8E21)

all_marg_dir=$1
all_ntoys=$2
all_cput=$3
all_mem=$4
all_reac=$5
all_job_dir=$6
all_par0=$7
all_par1=$8

for all_h in "${all_hier[@]}"; do
    for all_s in "${all_study_type[@]}"; do
        for all_a in "${all_allocation[@]}"; do
            for all_p in "${all_POT[@]}"; do
                . ${T2KOSCANA}/production/future_sensitivity_scripts/fit.sh ${all_s} ${all_a} ${all_p} ${all_marg_dir} ${all_ntoys} ${all_cput} ${all_mem} ${all_h} ${all_reac} ${all_job_dir} ${all_par0} ${all_par1}
	    done
        done
    done
done

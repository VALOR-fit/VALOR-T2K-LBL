mix23=0.528
dcp=-1.601
mix12=0.304
mix13=0.0219
dmsq12=7.53E-5
dmsq23=2.509E-3

command="${T2KOSCANA}/macros/tn_plotting_stuffs/make_table_systematic_effects.C"
std_args="\"3f_osc\",${mix12},${mix13},${mix23},${dcp},${dmsq12},${dmsq23},\"sinsqtheta\",\"normal\""
outfile="sk_error_table"
queue="long -M 5000"
queuex="long -M 8000"
#queue="short"

samples=(1Re 1Rmu MultiRe RHC1Re RHC1Rmu)
groups=(all allsyst osc banff_constrained xsec_nue xsec_ncoth xsec_eb xsec_ncgamma skdetfsi skdet)
for sample in ${samples[@]}; do
    for group in ${groups[@]}; do
	if [ "${group}" == "skdet" ]; then
	    config="systematic_configuration_run1-9d_skdet"
	elif [ "${group}" == "skfsi" ]; then
	    config="systematic_configuration_run1-9d_skfsi"
	else
	    config="analysis_configuration_oa2018_run1-9d_v1"
	fi
	echo root2k -b -q $command'('"\"${sample}\"",${std_args},"\"${config}\"","\"${group}\"","\"${outfile}_${sample}\""')'
	exec="${T2KOSCANA}/bin/root2k"
	args="-b -q $command'(\""${sample}"\","${std_args}",\""${config}"\","\"${group}"\",\""${outfile}"_"${sample}"\")'"
	python ${T2KOSCANA}/production/run_on_condor.py -e $exec -a "${args}" -l ${outfile}_${sample}_${group} --escape-quotes
	#root2k -b -q $command'('"\"${sample}\"",${std_args},"\"${config}\"","\"${group}\"","\"${outfile}_${sample}\""')'
    done
   #bsub -q ${queue} -oo log_syst_effects_${sample}_skdet.log ${T2KOSCANA}/production/warwick_submit_command.sh root2k -q -b "${command}(\"${sample}\",${std_args},\"systematic_configuration_run1-9d_skdet\",\"skdet\",\"${outfile}_${sample}\")"
   #bsub -q ${queue} -oo log_syst_effects_${sample}_skfsi.log ${T2KOSCANA}/production/warwick_submit_command.sh root2k -q -b "${command}(\"${sample}\",${std_args},\"systematic_configuration_run1-9d_skfsi\",\"skfsi\",\"${outfile}_${sample}\")"
   #bsub -q ${queue} -oo log_syst_effects_${sample}_banff_constrained.log ${T2KOSCANA}/production/warwick_submit_command.sh root2k -q -b "${command}(\"${sample}\",${std_args},\"analysis_configuration_oa2018_run1-9d_v1\",\"banff_constrained\",\"${outfile}_${sample}\")"
   #bsub -q ${queue} -oo log_syst_effects_${sample}_xsec_nue.log ${T2KOSCANA}/production/warwick_submit_command.sh root2k -q -b "${command}(\"${sample}\",${std_args},\"analysis_configuration_oa2018_run1-9d_v1\",\"xsec_nue\",\"${outfile}_${sample}\")"
   #bsub -q ${queue} -oo log_syst_effects_${sample}_xsec_ncoth.log ${T2KOSCANA}/production/warwick_submit_command.sh root2k -q -b "${command}(\"${sample}\",${std_args},\"analysis_configuration_oa2018_run1-9d_v1\",\"xsec_ncoth\",\"${outfile}_${sample}\")"
   #bsub -q ${queue} -oo log_syst_effects_${sample}_xsec_ncgamma.log ${T2KOSCANA}/production/warwick_submit_command.sh root2k -q -b "${command}(\"${sample}\",${std_args},\"analysis_configuration_oa2018_run1-9d_v1\",\"xsec_ncgamma\",\"${outfile}_${sample}\")"
   #bsub -q ${queue} -oo log_syst_effects_${sample}_xsec_eb.log ${T2KOSCANA}/production/warwick_submit_command.sh root2k -q -b "${command}(\"${sample}\",${std_args},\"analysis_configuration_oa2018_run1-9d_v1\",\"xsec_eb\",\"${outfile}_${sample}\")"
   #bsub -q ${queue} -oo log_syst_effects_${sample}_osc.log ${T2KOSCANA}/production/warwick_submit_command.sh root2k -q -b "${command}(\"${sample}\",${std_args},\"analysis_configuration_oa2018_run1-9d_v1\",\"osc\",\"${outfile}_${sample}\")"
   #bsub -q ${queuex} -oo log_syst_effects_${sample}_allsyst.log ${T2KOSCANA}/production/warwick_submit_command.sh root2k -q -b "${command}(\"${sample}\",${std_args},\"analysis_configuration_oa2018_run1-9d_v1\",\"allsyst\",\"${outfile}_${sample}\")"
   #bsub -q ${queuex} -oo log_syst_effects_${sample}_all.log ${T2KOSCANA}/production/warwick_submit_command.sh root2k -q -b "${command}(\"${sample}\",${std_args},\"analysis_configuration_oa2018_run1-9d_v1\",\"all\",\"${outfile}_${sample}\")"
done


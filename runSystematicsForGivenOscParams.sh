#\!/bin/bash
for syst in `cat syst.list`
do
jobID=0
      while [ $jobID -le 20 ]
      do
      echo $syst-$jobID
      echo "#\!/bin/bash">$syst-$jobID.job
      echo "cd ${PWD}" >>$syst-$jobID.job
      echo "source setup_at_ral_batch_farm.sh">>$syst-$jobID.job
      echo "root -l -b -q 'run_systematic_study_for_given_osc_params.C(\"$syst\", $jobID)' ">>$syst-$jobID.job
      chmod +x $syst-$jobID.job
      qsub -q prod -l cput=6:00:00 -l walltime=10:00:00 -j oe $syst-$jobID.job

      sleep 1
      let jobID=jobID+1    
done  
done

#\!/bin/bash
name="contour_newchisq"
jobID=0
while [ $jobID -le 250 ]
      do
      echo $name-$jobID
      echo "#\!/bin/bash">$name-$jobID.job
      echo "cd ${PWD}" >>$name-$jobID.job
      echo "source ../../../dev-setup">>$name-$jobID.job
      echo "root -l -b -q 'run_oscillation_contour_generator.C($jobID,\"$name\")' ">>$name-$jobID.job
      chmod +x $name-$jobID.job
      qsub -q prod -l cput=6:00:00 -l walltime=6:00:00 -j oe $name-$jobID.job

      sleep 2
      let jobID=jobID+1    
done

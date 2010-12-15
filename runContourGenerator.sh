#\!/bin/bash
name="contour_newchisq"
jobID=0
while [ $jobID -le 250 ]
      do
      echo $name-$jobID
      echo "#\!/bin/bash">$name-$jobID.job
      echo "cd ${PWD}" >>$name-$jobID.job
      echo "source setup-genie.sh">>$name-$jobID.job
      echo "source ../../../dev-setup">>$name-$jobID.job
      echo "root -l -b -q 'run_syst_contour_calc.C($jobID,\"$name\")' ">>$name-$jobID.job
      chmod +x $name-$jobID.job
      qsub -q prod -l cput=1:00:00 -l walltime=1:00:00 -j oe $name-$jobID.job

      sleep 2
      let jobID=jobID+1    
done

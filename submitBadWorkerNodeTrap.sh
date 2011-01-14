#\!/bin/bash
irun=0
while [ $irun -le 30 ]
do
   echo trap-$irun
   echo "#\!/bin/bash">trap-$irun.job
   echo "cd ${PWD}" >>trap-$irun.job
   echo "perl trap_bad_worker_nodes ">>trap-$irun.job
   chmod +x trap-$irun.job
   qsub -q prod -r y -l cput=12:00:00 -j oe trap-$irun.job
   let irun=irun+1    
done

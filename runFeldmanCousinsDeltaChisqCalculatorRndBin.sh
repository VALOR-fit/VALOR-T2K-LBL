#\!/bin/bash
name="geniefc_rndbiniter0"
irun=0
while [ $irun -le 500 ]
do
   echo $name-$irun
   echo "#\!/bin/bash">$name-$irun.job
   echo "cd ${PWD}" >>$name-$irun.job
   echo "source setup_at_ral_batch_farm.sh">>$name-$irun.job
   echo "root -l -b -q 'run_feldman_cousins_deltachisq_calculator_rndbin.C($irun,\"$name\")' ">>$name-$irun.job
   chmod +x $name-$irun.job
   qsub -q prod -r y -l cput=12:00:00 -j oe $name-$irun.job
   let irun=irun+1    
done

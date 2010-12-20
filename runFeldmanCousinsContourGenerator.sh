#\!/bin/bash
name="geniefc"
idmsqbin=1
while [ $idmsqbin -le 10 ]
do
	isinsq2thetabin=1
	while [ $isinsq2thetabin -le 10 ]
	do
		irun=0
		while [ $irun -le 10 ]
		do
		      echo $name-$irun-$idmsqbin-$isinsq2thetabin
		      echo "#\!/bin/bash">$name-$irun-$idmsqbin-$isinsq2thetabin.job
		      echo "cd ${PWD}" >>$name-$irun-$idmsqbin-$isinsq2thetabin.job
		      echo "source ../../../dev-setup">>$name-$irun-$idmsqbin-$isinsq2thetabin.job
		      echo "root -l -b -q 'run_oscillation_contour_generator.C($irun,\"$name\")' ">>$name-$irun-$idmsqbin-$isinsq2thetabin.job
		      chmod +x $name-$irun-$idmsqbin-$isinsq2thetabin.job
		      qsub -q prod -l cput=6:00:00 -l walltime=6:00:00 -j oe $name-$irun-$idmsqbin-$isinsq2thetabin.job
		      sleep 1
		      let irun=irun+1    
		done
		let isinsq2thetabin=isinsq2thetabin+1
	done
	let idmsqbin=idmsqbin+1
done

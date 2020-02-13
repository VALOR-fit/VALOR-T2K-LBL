nmargtoys=50
param1nbins=32
maxbinx=7
param2nbins=22
maxbiny=5

for ((i=0; i<=maxbinx; i++)); do
   stepx=$(echo "2 * 3.141592654 / ${maxbinx}" | bc -l)
   cpv=$(bc <<< "-3.141592654+${i}*${stepx}")
   for ((j=0; j<=maxbiny; j++)); do
      stepy=$(echo "0.4 / ${maxbiny}" | bc -l)
      mixing23=$(bc <<< "0.3+${j}*${stepy}")
      valor_dune-run_fit_marg \
         --config $VALOR_EXPT/config/test/marginalisation/analysis.xml \
         --seed-toys-syst 1 \
         --seed-toys-osc 1 \
         --seed-other 1 \
         --nmargtoy $nmargtoys \
         --nbatches 1 \
         --fit-param1 cpv \
         --fit-param1-value $cpv \
         --fit-param1-bin 1 \
         --fit-param1-min-value -3.141592654 \
         --fit-param1-max-value +3.141592654 \
         --fit-param1-nbins $param1nbins \
         --fit-param2 mixing23 \
         --fit-param2-value $mixing23 \
         --fit-param2-bin 1 \
         --fit-param2-min-value 0.3 \
         --fit-param2-max-value 0.7 \
         --fit-param2-nbins $param2nbins \
         --tag fits/fit_${i}_${j}_cpv_mixing23 \
         --margtoy-tag toys/2d/margtoy_cpv_mixing23
   done
done


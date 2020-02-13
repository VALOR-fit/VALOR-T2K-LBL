maxbinx=31
maxbiny=21
for ((i=0; i<=maxbinx; i++)); do
   stepx=$(echo "2 * 3.141592654 / ${maxbinx}" | bc -l)
   cpv=$(bc <<< "-3.141592654+${i}*${stepx}")
   for ((j=0; j<=maxbiny; j++)); do
      >&2 echo "Building $i $j"
      stepy=$(echo "0.4 / ${maxbiny}" | bc -l)
      theta23=$(bc <<< "0.3+${j}*${stepy}")
      valor_dune-run_build_margtoys \
         --config $VALOR_EXPT/config/test/marginalisation/analysis.xml \
         --seed-toys-syst 1 \
         --seed-toys-osc 1 \
         --seed-other 1 \
         --nmargtoy 50 \
         --batch 0 \
         --fit-param1 cpv \
         --fit-param1-value $cpv \
         --fit-param1-bin $i \
         --fit-param2 mixing23 \
         --fit-param2-value $theta23 \
         --fit-param2-bin $j \
         --tag toys/2d/margtoy_cpv_mixing23
   done
done

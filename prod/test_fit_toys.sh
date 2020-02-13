nmargtoys=50
param1nbins=101
maxbin=100
for ((i=0; i<=maxbin; i++)); do
   step=$(echo "2 * 3.141592654 / ${maxbin}" | bc -l)
   cpv=$(bc <<< "-3.141592654+${i}*${step}")
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
      --tag fits/fit_${i}_cpv \
      --margtoy-tag toys/1d/margtoy_cpv
done


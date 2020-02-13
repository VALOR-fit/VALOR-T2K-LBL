maxbin=100
for ((i=0; i<=maxbin; i++)); do
   echo $i
   step=$(echo "2 * 3.141592654 / ${maxbin}" | bc -l)
   cpv=$(bc <<< "-3.141592654+${i}*${step}")
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
      --tag toys/1d/margtoy_cpv
done

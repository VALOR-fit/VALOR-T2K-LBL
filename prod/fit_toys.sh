nmargtoys=50
param1nbins=101
cpv=-1.570796327
valor_dune-run_fit_marg \
   --config $VALOR_EXPT/config/test/marginalisation/analysis.xml \
   --seed-toys-syst 1 \
   --seed-toys-osc 1 \
   --seed-other 1 \
   --nmargtoy $nmargtoys \
   --nbatches 1 \
   --fit-param1 cpv \
   --fit-param1-value $cpv \
   --fix-at-cpv $cpv \
   --fit-param1-bin 1 \
   --fit-param1-min-value -3.141592654 \
   --fit-param1-max-value +3.141592654 \
   --fit-param1-nbins $param1nbins \
   --tag fits/fit_cpv \
   --margtoy-tag toys/1d/margtoy_cpv


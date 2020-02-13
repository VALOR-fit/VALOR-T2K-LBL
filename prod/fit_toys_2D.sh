nmargtoys=50
param1nbins=32
param2nbins=22
cpv=-1.570796327
mixing23=0.53
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
   --tag fits/fit_cpv_mixing23 \
   --margtoy-tag toys/2d/margtoy_cpv_mixing23


for file in $VALOR_DATA/official/2021/superk/mc/Minituples_21b_031221_EMiss_withcuts/t2ksk19b.fqv4r0.21bv2.*.root; do
       root -n -l -b -q 'validate_mtuples_OA2021_hard_cuts.C(" '${file}' ")'
done

for file in $VALOR_DATA/official/2021/superk/mc/Minituples_21b_031221_EMiss_withcuts/t2ksk19b.fqv4r0.21bv2.*selec.root; do
       root -n -l -b -q "ApplyCutsToMinituples.C(\"${file}\", \"${file/.root}_cut.root\")"
done

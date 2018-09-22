#!/bin/bash
str="$@"
tokens=($str)
newarray=()

for (( i=0; i<${#tokens[@]}; i++ )); do
   if [ "${tokens[i]}" == "--file-margtoys-nh" ]; then
      nh="${tokens[i + 1]//\"/}"
      x="${str%%$nh*}"
      val=$(echo $nh | sed 's/.*\///')
      newdir="${TMPDIR}/$val"
      newarray+="${tokens[i]} "
      newarray+="$newdir "
      i=$((i+1))
   elif [ "${tokens[i]}" == "--file-margtoys-ih" ]; then
      ih="${tokens[i + 1]//\"/}"
      x="${str%%$ih*}"
      val=$(echo $ih | sed 's/.*\///')
      newdir="${TMPDIR}/$val"
      newarray+="${tokens[i]} "
      newarray+="$newdir "
      i=$((i+1))
   elif [ "${tokens[i]}" == "--file-tag" ]; then
      newdir="${TMPDIR}/${tokens[i + 1]}"
      out="$newdir"
      newarray+="${tokens[i]} "
      newarray+="$newdir "
      i=$((i+1))
   elif [[ ${tokens[i]} == \-\-* ]]; then
      if [[ ${tokens[i + 1]} == \-\-* ]]; then
         newarray+="${tokens[i]} "
         newarray+="\"\" "
      else
         newarray+="${tokens[i]} "
      fi
   else
      newarray+="${tokens[i]} "
   fi
done

newarray+="--top-level-data-dir "
newarray+="${TMPDIR}/"

cp $nh*.bin ${TMPDIR}
cp $ih*.bin ${TMPDIR}
#rsync -aq $ih*.bin ${TMPDIR}

mkdir -p ${TMPDIR}/official/2018/banff
mkdir -p ${TMPDIR}/official/2018/flux
mkdir -p ${TMPDIR}/official/2018/skdetfsi
mkdir -p ${TMPDIR}/official/2018/superk/mc/mtuples_v2/valor_format
mkdir -p ${TMPDIR}/official/2018/xsec/splines2018v1/valor_format

origdatadir=$T2KOSCANA/inputs/official/2018

cp $origdatadir/banff/*valor_order_all_plus_scc.root ${TMPDIR}/official/2018/banff/
cp -r $origdatadir/flux/tuned13av3 ${TMPDIR}/official/2018/flux/
cp -r $origdatadir/skdetfsi/*VALOR_order.root ${TMPDIR}/official/2018/skdetfsi/
cp -r $origdatadir/superk/mc/mtuples_v2/valor_format/*.root ${TMPDIR}/official/2018/superk/mc/mtuples_v2/valor_format/
cp -r $origdatadir/xsec/splines2018v1/valor_format/*.root ${TMPDIR}/official/2018/xsec/splines2018v1/valor_format/

cmd="${newarray[*]}"

eval "$cmd"

echo -e "The command $cmd"

cp $out*.root .

rm -rf ${TMPDIR}/official
rm -f ${TMPDIR}/*


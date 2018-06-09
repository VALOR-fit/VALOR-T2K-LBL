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
   elif [ "${tokens[i]}" == "--file-margtoys-ih" ]; then
      ih="${tokens[i + 1]//\"/}"
      x="${str%%$ih*}"
      val=$(echo $ih | sed 's/.*\///')
      newdir="${TMPDIR}/$val"
      newarray+="${tokens[i]} "
      newarray+="$newdir "
   elif [ "${tokens[i]}" == "--file-tag" ]; then
      newdir="${TMPDIR}/${tokens[i + 1]}"
      out="$newdir"
      newarray+="${tokens[i]} "
      newarray+="$newdir "
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

rsync -aq $nh*.bin ${TMPDIR}
rsync -aq $ih*.bin ${TMPDIR}

cmd="${newarray[*]}"

eval "$cmd"

echo -e $cmd

cp $out*.root .

rm -f ${TMPDIR}/*.bin ${TMPDIR}/*.root


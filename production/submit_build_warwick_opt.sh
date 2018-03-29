#!/bin/bash

str="$@"
IFS=',' read -r -a tokens <<< "$@"
newarray=()

for (( i=0; i<${#tokens[@]}; i++ )); do
   if [ $i == 6 ]; then
      # file output tag, prefix node directory (need to strip and readd quotes)
      tag=`eval echo -e ${tokens[i]}`
      newarray+="\"${TMPDIR}/${tag}\","
   elif [[ ${tokens[i]} == *[\)] ]]; then
      # Closing bracket, no comma afterwards
      newarray+="${tokens[i]}"
   else 
      # Standard argument, add a comma afterwards
      newarray+="${tokens[i]},"
   fi
done

cmd="${newarray[*]}"
eval "$cmd"
rsync -aq ${TMPDIR}/* .
rm -f ${TMPDIR}/*


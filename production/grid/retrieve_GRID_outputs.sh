#!/bin/bash

GRIDPATH=/grid/t2k.org/nd280/contrib/dealtry/2012/;
LOCALPATH=./;

voms-proxy-init -voms t2k.org;
export LFC_HOST=lfc.gridpp.rl.ac.uk;
lfc-ls $GRIDPATH > out;

FILENAME="out";
while read name;
do 
    if [ ! -f $LOCALPATH/$name ]
    then
	lcg-cp -v --vo t2k.org lfn:$GRIDPATH/$name $LOCALPATH/$name;
    fi
    #if [ -f $LOCALPATH/$name ]
    #then
	#lcg-del  -v --vo t2k.org -a lfn:$GRIDPATH/$name;
    #fi
done < $FILENAME;
rm -fr $FILENAME;

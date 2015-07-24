voms-proxy-init -voms t2k.org:/Role=production
export LFC_HOST=lfc.gridpp.rl.ac.uk
lfc-ls /grid/t2k.org/contrib/lorena/outputs/sensitivity > out

FILENAME="out"
while read name
do 
    lcg-cp -v --vo t2k.org lfn:/grid/t2k.org/contrib/lorena/outputs/sensitivity/$name GRID_outputs/sensitivity/$name
    lcg-del  -v --vo t2k.org -a lfn:/grid/t2k.org/contrib/lorena/outputs/sensitivity/$name
done < $FILENAME
rm -fr $FILENAME

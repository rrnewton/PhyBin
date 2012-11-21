#!/bin/bash
set -e 
if [ "$1" == "" ]; then echo "ERR: Two args required!"; exit 1; fi
if [ "$2" == "" ]; then echo "ERR: Two args required!"; exit 1; fi

taxa=$1
dataset=$2
branchlens="0.0001 0.001 0.01 0.05 0.1 0.5 1.0 5.0"

echo Wiping "$dataset/phybin_output/" ...
rm -rf "$dataset/phybin_output/"

echo "[parallel] Starting jobs..."

for bl in $branchlens; do
    echo "Running with branch len $bl"
    echo "=============================="
    PHYBIN="./phybin.exe -b $bl -g -n $taxa"
    outdir="$dataset/phybin_output/brlen$bl/"
    LOG="$dataset/phybin_output/brlen$bl.log"

    mkdir -p $outdir
    
    # ASSUMPTION1! We need the renaming table at a particular spot:
    # ASSUMPTION2! We likewise expect the inputs in a certain spot.
    RENAMER="$dataset/renaming_table.txt"
    echo "  RUNNING:  $PHYBIN -m $RENAMER -s '_'  -o $outdir"
    $PHYBIN -m $RENAMER -s '_'  -o $outdir $dataset/final_trees/*BranchLab*.out | tee $LOG
done

echo "[parallel] Waiting for outstanding jobs"
FAIL=0
for job in `jobs -p`
do
echo $job
    wait $job || let "FAIL+=1"
done
if [ "$FAIL" == "0" ];
then
echo "All jobs completed successfully"
else
echo "ERROR some jobs ($FAIL) failed!"
exit 1 
fi

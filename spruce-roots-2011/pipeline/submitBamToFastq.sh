#!/bin/bash

proj="umea"
in=/mnt/picea/projects/spruce/14_SpruceRoots_Project/STAR
out=/mnt/picea/projects/spruce/14_SpruceRoots_Project/trimmomatic
mail="nicolas.delhomme@umu.se"

## stop on error
set -ex

## check for UPSCb
if [ -z $UPSCb ]; then
    echo "The UPSCb env. var. needs to be set to your UPSCb Git checkout directory"
    exit 1;
fi

## run
for f in `find $in -name "*.bam" -type f`
do 
    fnam=`basename ${f//.bam/}`
    sbatch -A $proj --mail-user $mail -e $out/$fnam.err -o $out/$fnam.out -J b2f-$fnam $UPSCb/pipeline/runBedToolsBamToFastq.sh $f $out/${fnam}_1.fq $out/${fnam}_2.fq
done


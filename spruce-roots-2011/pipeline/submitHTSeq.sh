#!/bin/bash

proj=b2011227
tmp=/glob/delhomme/tmp
in=/proj/$proj/spruceRoots/STAR
out=/proj/$proj/spruceRoots/HTSeq
gff3=/proj/b2011227/reference/gff3/Eugene.gff3
mail="nicolas.delhomme@umu.se"

## stop on error
set -e

## check for UPSCb
if [ -z $UPSCb ]; then
    echo "The UPSCb env. var. needs to be set to your UPSCb Git checkout directory"
    exit 1;
fi

## create out dir
if [ ! -d $out ]; then
    mkdir $out
fi

## run
for f in `find $in -name "*.bam" -type f`
do 
    fnam=`basename ${f//.bam/}`
    sbatch -A $proj --mail-user $mail -C "usage_mail" -e $out/$fnam.err -o $out/$fnam.out -J htseq-$fnam $UPSCb/pipeline/runHTSeq.sh $out $tmp $f $gff3
done


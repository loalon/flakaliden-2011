#!/bin/bash

## global vars
proj=b2011227
mail="nicolas.delhomme@umu.se"
in=/proj/$proj/nobackup/spruceRoots/raw
out=/proj/$proj/nobackup/spruceRoots/fastQValidator

## raw
for f in `find $in -name "*.fastq.gz"`; do
fnam=`basename $f`
sbatch -A $proj --mail-user $mail -e $out/${fnam//.fastq.gz/.err} -o $out/${fnam//.fastq.gz/.out} -J fqv-$fnam ../../../pipeline/runFastQValidator.sh $f
done

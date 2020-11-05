#!/bin/bash

## args
mail="nicolas.delhomme@umu.se"
proj=b2011227
in=/proj/$proj/nobackup/spruceRoots/raw
out=/proj/$proj/nobackup/spruceRoots/trimmomatic

echo "Only for re-running the 1_120626_BC0YAHACXX_229_* samples"

for f in `find $in -name "1_120626_BC0YAHACXX_229_*"`; do echo `basename ${f//_[1,2].fastq.gz*/}` ; done | sort | uniq | while read line;
do sbatch -A b2010042 --mail-user $mail -e $out/$line.err -o $out/$line.out -J Trim-$line ../../../pipeline/runTrimmomatic.sh $in/${line}_1.fastq.gz $in/${line}_2.fastq.gz $out LEADING:20 SLIDINGWINDOW:5:20 MINLEN:50
done


#!/bin/bash

## args
mail="nicolas.delhomme@umu.se"
proj=b2011227
in=/proj/$proj/nobackup/spruceRoots/raw
out=/proj/$proj/nobackup/spruceRoots/trimmomatic

for f in `ls $in`; do echo `basename ${f//_[1,2].fastq.gz*/}` ; done | sort | uniq | while read line;
do sbatch -A $proj --mail-user $mail -e $out/$line.err -o $out/$line.out -J Trim-$line ../../../pipeline/runTrimmomatic.sh $in/${line}_1.fastq.gz $in/${line}_2.fastq.gz $out SLIDINGWINDOW:5:20 MINLEN:50
done


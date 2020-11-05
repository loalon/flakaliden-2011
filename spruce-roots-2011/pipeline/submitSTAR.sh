#!/bin/bash -l

## global vars
proj=b2010042
mail="nicolas.delhomme@umu.se"
in=/proj/$proj/nobackup/est/seq_data/illumina/Z4006_BGI/sortmerna
out=/proj/$proj/nobackup/est/seq_data/illumina/Z4006_BGI/STAR
intron=70000
gff=/proj/$proj/annotation/GenePrediction/ASSEMBLYLOCK_2012_NOV/PREDICTIONS_2013_JAN/Sprucev02_combined.gtf
genome=/proj/$proj/assembly/ASSEMBLYLOCK_2012_NOV/indices/STAR

## raw
for f in `find $in -name "*.fq.gz" -type f`; do echo `basename ${f//_[1,2].fq.gz/}` ; done | sort | uniq | while read line;
do
sbatch -A $proj --mail-user $mail -e $out/$line.err -o $out/$line.out -J STAR-$line ../../../pipeline/runSTAR.sh -o $out -m $intron $in/${line}_1.fq.gz $in/${line}_2.fq.gz $genome $gff -- --outQSconversionAdd -31 --outReadsUnmapped Fastx
done

#!/bin/bash

proj=b2011227
tmp=/glob/delhomme/tmp
mail="nicolas.delhomme@umu.se"
in=/proj/$proj/nobackup/spruceRoots/trimmomatic
out=/proj/$proj/spruceRoots/sortmerna

## check that the needed data is in the user glob
if [ ! -e ../../../data/sortmerna ]; then
    if [ ! -d /glob/$USER/data/sortmerna ]; then
	mkdir -p /glob/$USER/data/sortmerna
	echo -n "Enter your username at picea.plantphys.umu.se: "
	read picea
	scp -P 922 -r $picea@picea.plantphys.umu.se:/mnt/picea/storage/reference/rRNA/sortmerna /glob/$USER/data/sortmerna
    fi
    ln -sf /glob/$USER/data/sortmerna -T ../../../data/sortmerna
else
    if [ ! -L ../../../data/sortmerna ]; then
	echo "Aborting as we do not want to overwrite data." `pwd`"/../../../data/sortmerna should be a symbolic link."
	exit 1
    fi
fi

## then start the script
for f in `find $in -name "*trimmomatic_[1,2].fq.gz" -type f`; do echo "${f//_trimmomatic_[1,2].fq.gz/}" ; done | sort | uniq | while read line;
do 
fnam=`basename $line`
sbatch -A $proj --mail-user $mail -e $out/$fnam.err -o $out/$fnam.out -J smr-$fnam ../../../pipeline/runSortmerna.sh $out $tmp ${line}_trimmomatic_1.fq.gz ${line}_trimmomatic_2.fq.gz
done


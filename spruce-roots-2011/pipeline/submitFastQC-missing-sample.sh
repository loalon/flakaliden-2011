#!/bin/bash

## vars
proj=b2011227
mail="nicolas.delhomme@umu.se"

## check the arguments
if [ $# != 1 ]; then
   echo "This script takes one argument: the tool name"
   exit 1
fi

tool=
nam="*.fq.gz"
case "$1" in
    raw) 
	tool=raw
	nam="*.fastq.gz";;
    trimmomatic) 
	tool=trimmomatic;;
    sortmerna) 
	tool=sortmerna;;
esac

if [ -z $tool ]; then
    echo "The second argument should be one of raw, trimmomatic or sortmerna"
    exit 1
fi  

## check $in
in=/proj/$proj/nobackup/spruceRoots/$tool/
if [ ! -d $in ]; then
    echo "The tool directory $1 does not exist. Make sure that the tool was run first."
    exit 1
fi

## create the dir
out=/proj/$proj/nobackup/spruceRoots/FastQC/$tool
mkdir -p $out

## a warning
echo "Only for running the 2_120823_AD16FAACXX_229_197_index27_* samples"
nam=2_120823_AD16FAACXX_229_197_index27_$nam

for file in `find $in -name $nam`; do 
fnam=`basename $file`
fnam=${fnam//.f*q.gz/}
sbatch -A b2010042 --mail-user $mail -e $out/$fnam.err -o $out/$fnam.out -J FastQC-$fnam ../../../pipeline/runFastQC.sh $out $file
done
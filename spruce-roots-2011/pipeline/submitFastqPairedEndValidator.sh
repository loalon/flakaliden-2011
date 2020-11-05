#!/bin/bash -l

## global vars
proj=b2011227
mail="nicolas.delhomme@umu.se"

## error
set -e

## check the arguments
if [ $# != 1 ]; then
   echo "This script takes one argument: the tool name"
   exit 1
fi

tool=
nam="*.fq.gz"
ext="fq"
case "$1" in
    raw) 
	tool=raw
	nam="*.fastq.gz"
	ext="fastq";;
    trimmomatic) 
	tool=trimmomatic;;
    sortmerna) 
	tool=sortmerna;;
esac

if [ -z $tool ]; then
    echo "The first argument should be one of raw, trimmomatic or sortmerna"
    exit 1
fi  

## check $in
in=/proj/$proj/nobackup/spruceRoots/$tool/
if [ ! -d $in ]; then
    echo "The tool directory $1 does not exist. Make sure that the tool was run first."
    exit 1
fi

## create the dir
out=/proj/$proj/nobackup/spruceRoots/FastqPairedEndValidator/$tool
mkdir -p $out

for f in `find $in -name $nam`; do echo `basename ${f//_[1,2].f*q.gz/}` ; done | sort | uniq | while read line;
do
sbatch -A $proj -p core -n 1 -t 01:00:00 --mail-type=ALL --mail-user $mail -e $out/$line.err -o $out/$line.out -J FPEV-$line ../../../src/perl/FastqPairedEndValidator.pl $in/${line}_1.$ext.gz $in/${line}_2.$ext.gz
done

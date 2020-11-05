#!/bin/bash

## vars
proj=b2011227
mail="nicolas.delhomme@umu.se"

## check input
in=/proj/$proj/spruceRoots/sortmerna
out=/proj/$proj/nobackup/spruceRoots/sortmerna
tmp=/glob/delhomme/tmp

## find the files
cd $in
for f in `find -type f`
do
    if [ ! -d $out/`dirname $f` ]; then
	mkdir -p $out/`dirname $f`
    fi
    sbatch -A $proj -e $tmp/$f.err -o $tmp/$f.out -J trf-$f $UPSCb/pipeline/runTransfer.sh $in/$f $out/$f
done

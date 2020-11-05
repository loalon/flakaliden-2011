#!/bin/bash -l

proj=b2011227
in=/proj/$proj/spruceRoots/STAR
tmp=/glob/delhomme/tmp

for f in `find $in -name "*mate[1,2]"`
do
fnam=`basename $f`
sbatch -A $proj -e $tmp/$fnam.err -o $tmp/$fnam.out  -J gzip-$fnam $UPSCb/pipeline/runGzip.sh $f
done

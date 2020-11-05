#!/bin/bash -l

## define a function
usage () {
    echo "This function take one argument as parameter; one of 'umea','uppmax'"
    echo "The UPSCb env. var. needs to be set to your Git UPSCb checkout directory."
    exit 1
}

## args number
if [ $# != 1 ]; then
    usage
fi

## process the argument
runLoc=
case "$1" in
    umea)
	runLoc=$1
	in=/mnt/picea/projects/spruce/14_SpruceRoots_Project/STAR
	out=/mnt/picea/projects/spruce/14_SpruceRoots_Project/diginorm/raw/others
	cfgNum=1
	cfg=/mnt/picea/projects/spruce/14_SpruceRoots_Project/cfg
	;;
    uppmax) 
	echo Not implemented yet.
	exit 1
	;;
    *) usage;;
esac

## check vars
if [ -z $UPSCb ]; then
    echo "The UPSCb var needs to be set."
    usage
fi

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## clean
if [ -f $cfg/$cfgNum.in ]; then
 rm $cfg/$cfgNum.in
fi

## process
if [ $runLoc == "umea" ]; then
    for f in `find $in -type f -name "*mate1.gz"`; do
	echo bash $UPSCb/pipeline/runRePair.sh $f ${f//mate1/mate2} >> $cfg/$cfgNum.in
    done
    cd $cfg
    runParallel.pl -t 32 -l $cfgNum -m $cfgNum > log-$cfgNum.out 2> log-$cfgNum.err
else
    echo some sbatch
fi


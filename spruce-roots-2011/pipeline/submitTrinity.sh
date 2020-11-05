#!/bin/bash
set -ex

## usage
usage () {
    echo "This function takes one argument: either 'unmapped' or 'all'"
    echo "The UPSCb env. var. needs to be set to your Git UPSCb checkout directory."
    echo "The TRINITY_RNASEQ_ROOT env. var. needs to be set to the root directory of your trinity installation"
    exit 1
}

#Check arg number
if [ $# -ne 1 ];then
    usage
fi

case "$1" in
    unmapped)
	in=/mnt/picea/projects/spruce/14_SpruceRoots_Project/trinity/unmapped
	out=/mnt/picea/projects/spruce/14_SpruceRoots_Project/trinity/unmapped-k1
	;;
    all)
	in=/mnt/picea/projects/spruce/14_SpruceRoots_Project/trinity/all
	out=/mnt/picea/storage/projects/14_SpruceRoots_Project/trinity/all-k1
	;;
    *) usage;;
esac

# error checks
if [ ! -d $in ];then
    usage
fi

# check UPSCb 
if [ -z $UPSCb ]; then
    echo "The UPSCb env. var. needs to be set."
    usage
fi

# create out dir
if [ ! -d $out ];then
    mkdir -p $out
fi

#Generate comma separated lists of forward and reverse reads
export TRINITY_RNASEQ_ROOT=~/opt/trinityrnaseq_r20140413

case "$1" in
    unmapped)
	bash $UPSCb/pipeline/runTrinity.sh -b -n -p 64 -m 400G $out `find $in -name '*.mate1.gz' -printf %p\, | sed "s/,$//g"` `find $in -name '*.mate2.gz' -printf %p\, | sed "s/,$//g"`
	;;
    all)
	bash $UPSCb/pipeline/runTrinity.sh -b -n -p 64 -m 400G $out \
	$(find $in -name "*sortmerna_1.fq.gz" | sort | paste -s --delimiters=,) \
	$(find $in -name "*sortmerna_2.fq.gz" | sort | paste -s --delimiters=,) \
	$(find $in -name "*sunpaired_[1,2].fq.gz" | sort | paste -s --delimiters=,)
	;;
esac



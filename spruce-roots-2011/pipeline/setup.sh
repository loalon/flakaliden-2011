#!/bin/sh

## create the dirs
cd /proj/b2011127/nobackup
mkdir -p spruceRoots/raw
mkdir -p spruceRoots/trimmomatic
mkdir -p spruceRoots/sortmerna
mkdir -p spruceRoots/STAR
mkdir -p spruceRoots/FastQC/raw
mkdir -p spruceRoots/FastQC/trimmomatic
mkdir -p spruceRoots/FastQC/sortmerna
mkdir -p spruceRoots/fastQValidator

## link the raw data
cd spruceRoots/raw
ln -s /proj/b2011227/sequence_data/raw/spruceRoots/N.Street_12_01/120626_BC0YAHACXX/*.fastq.gz .
ln -s /proj/b2011227/sequence_data/raw/spruceRoots/N.Street_12_01/120704_BD14F1ACXX/*.fastq.gz .
ln -s /proj/b2011227/sequence_data/raw/spruceRoots/N.Street_11_13/120425_sN188_0274_BD0TRAACXX/*.fastq.gz .

## or cat them
rm 2_120425_BD0TRAACXX_P191_10_index10_*
zcat /proj/b2011227/sequence_data/raw/spruceRoots/N.Street_11_13/120425_sN188_0274_BD0TRAACXX/2_120425_BD0TRAACXX_P191_10_index10_1.fastq.gz /proj/b2011227/sequence_data/raw/spruceRoots/N.Street_11_13/120602_AC0RC2ACXX/2_120602_AC0RC2ACXX_P191_10_index10_1.fastq.gz > P191_10_index10_1.fastq.gz

zcat /proj/b2011227/sequence_data/raw/spruceRoots/N.Street_11_13/120425_sN188_0274_BD0TRAACXX/2_120425_BD0TRAACXX_P191_10_index10_2.fastq.gz /proj/b2011227/sequence_data/raw/spruceRoots/N.Street_11_13/120602_AC0RC2ACXX/2_120602_AC0RC2ACXX_P191_10_index10_2.fastq.gz > P191_10_index10_2.fastq.gz


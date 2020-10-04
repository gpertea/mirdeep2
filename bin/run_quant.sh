#!/bin/bash
dir='/mnt/compute_data/smallRNA/smallRNA_platesALL/cutadapt_trimmed_cleaned'
USAGE="Usage: \n\
 run_mirDeep.sh <sample_prefix>\n\
 \n\
 Using input base directory: $dir\n\
 Creates a <sample_prefix> directory in the working directory.\n\
"

if [ $# -ne 1 ]; then
  echo -e $USAGE
  exit 1
fi

if [[ ! -f hairpin_hsa.fa ]]; then
 echo "hairpin_hsa.fa not found in current directory!"
 exit 1
fi

if [[ ! -f mature_hsa.fa ]]; then
 echo "mature_hsa.fa not found in current directory!"
 exit 1
fi

pre=$1
fp=$dir/${pre}_cutadapt-trimmed_cleaned.fastq.gz
if [[ ! -f $fp ]]; then
 echo "Error: $fp missing?"
 exit 1
fi

rf=reads_${pre}.collapsed.fa

set -e
if [[ ! -d $pre ]]; then
 mkdir $pre
fi

cd $pre

mapper.pl $fp -e -h -l 18 -v -m -s $rf -O m
quantifier.pl -d -p ../hairpin_hsa.fa -m ../mature_hsa.fa -r $rf -t hsa -y r


#!/bin/bash
dir='/mnt/compute_data/smallRNA/smallRNA_platesALL/cutadapt_trimmed_cleaned'
USAGE="Usage: \n\
 run_quant_mod.sh <sample_prefix>\n\
 \n\
 Using input base directory: $dir\n\
 Creates a <sample_prefix> directory in the working directory.\n\
"
bwt=./bwt

if [ $# -ne 1 ]; then
  echo -e $USAGE
  exit 1
fi

bwt=$(readlink -f $bwt)

if [[ ! -d $bwt ]]; then
 echo "$bwt bowtie index directory not found!"
 exit 1
fi

if [[ ! -f $bwt/miRNA_mature.1.ebwt ]]; then
 echo "$bwt/miRNA_mature bowtie index not found!"
 exit 1
fi

#if [[ ! -f $bwt/miRNA_precursor.1.ebwt ]]; then
# echo "$bwt/miRNA_precursor bowtie index not found!"
# exit 1
#fi

pre=$1

rf=reads_${pre}.collapsed.fa

set -e
if [[ ! -d $pre ]]; then
  echo "Error: $pre file not found!"
fi

cd $pre

if [[ ! -f $rf ]]; then
 echo "Error: $rf file not found!"
 exit 1
fi

#mapper.pl $fp -e -h -l 18 -v -m -s $rf -O m
echo "execute: quant_mod.pl -d -B $bwt -r $rf -y r"
quant_mod.pl -d -B $bwt -r $rf -y r > $pre.quant.log 2>&1

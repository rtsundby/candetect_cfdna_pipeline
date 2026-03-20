#!/usr/bin/env bash
input=$1
genome_fasta=$2
output=$3
threads=$4

module load picard
module load java

java -jar -Xmx128g -XX:ParallelGCThreads=32 $PICARDJARPATH/picard.jar CollectWgsMetrics\
       INPUT=$input \
       OUTPUT=$output \
       READ_LENGTH=150 \
       REFERENCE_SEQUENCE=$genome_fasta

#!/bin/sh
for file in 22010/22010_RawData/*fastq.gz*
do
  fastqc "$file" --outdir=22010/FASTQC
done

# run multiqc afterwards to compile all results
multiqc 22010/FASTQC --export
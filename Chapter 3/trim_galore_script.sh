#!/bin/sh

#Set paths
raw=~/Gentianella/whole_data/SRA_Sequence/
trim_out=~/Gentianella/whole_data/SRA_Sequence/trimmed_100BP/

#Do the trimming!
for samp in ${raw}*1.fastq.gz
do
base=$(basename ${samp} .1.fastq.gz)
echo "Running Trim_galore for ${base}..."
trim_galore --paired \
    --cores 8 \
    --length 100 \
    --output_dir ${trim_out}/ \
    --no_report_file \
    ${raw}${base}.1.fastq.gz \
    ${raw}${base}.2.fastq.gz
done

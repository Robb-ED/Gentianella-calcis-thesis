#!/bin/bash

# Path to working directory
work=WORKING DIRECTORY GOES HERE

# Parameter end points
ENDM=12
ENDi=9

# Loop over the M and m values
# Iterates over every value between 1-12 for M and 3, 6, 9 for m
for M in $(seq 1 $ENDM)
do
    for i in $(seq 3 3 $ENDi)
    do
    # This creates a new output directory per Stacks run
    out=$work/denovo.M${M}.m${i}
    mkdir $out
    # Move into the new directory
    cd $out
        # Stacks command
        denovo_map.pl --samples $work --popmap $work/popmap.tsv \
            --out-path $out \
            --paired \
            -M $M \
            -n $M \
            -T 20 \
            --min-samples-per-pop 0.8 \
            -X "ustacks: --force-diff-len" \
            -X "ustacks: -m ${i}"

    done
done

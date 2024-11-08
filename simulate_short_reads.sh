#!/bin/bash

# Configuration Parameters
GENOME_FA="test_data/genome.fa"
READ_LENGTH=100
COVERAGE=10
NUM_SAMPLES=6
OUTPUT_DIR="test_data/short_reads"

# Create output directory
mkdir -p ${OUTPUT_DIR}

for i in $(seq 1 ${NUM_SAMPLES}); do
    SAMPLE="Sample${i}"
    echo "Simulating short reads for ${SAMPLE}"
    art_illumina -ss HS25 -i ${GENOME_FA} -p -l ${READ_LENGTH} -f ${COVERAGE} -m 200 -s 20 -o ${OUTPUT_DIR}/${SAMPLE}
    # Compress the FASTQ files
    gzip ${OUTPUT_DIR}/${SAMPLE}1.fq
    gzip ${OUTPUT_DIR}/${SAMPLE}2.fq
done

echo "Short-read simulation completed."
